## ====================== packages ======================
pacman::p_load(data.table,impute,RiskScorescvd)

## ====================== paths & data ==================
ckm_path <- "D:/ckm2/data/phe/ckm.rds"
dat      <- readRDS(ckm_path)
dat      <- as.data.table(dat)

## ====================== 0) basic / categorical filter ==
# basic continuous predictors
basic_covs <- c("age", "sex", "ethnic.c", paste0("PC", 1:10))
basic_covs <- intersect(basic_covs, names(dat))

# key categorical variables
cat_vars <- c("sex", "ethnic.c", "smoke_status", "alcohol_status","drug.lipid","drug.htn","drug.dm")
cat_vars <- intersect(cat_vars, names(dat))

# N before and after dropping missing basic / categorical
  cat("N before basic/categorical filter:", nrow(dat), "\n")
 # N before basic/categorical filter: 48175 
keep_cov <- c(basic_covs, cat_vars)
dat <- dat[complete.cases(dat[, ..keep_cov])]
  cat("N after basic/categorical filter:", nrow(dat), "\n")
 # N after basic/categorical filter: 47372 

## ====================== 1) metabolite KNN imputation ===
met_cols <- grep("^met_", names(dat), value = TRUE)

if (length(met_cols) > 0) {
  # count before
    cat("metabolites before:", length(met_cols), "\n")
    cat("individuals before:", nrow(dat), "\n")
   # metabolites before: 251 
   # individuals before: 47372 
  # filter >10% missing metabolites & individuals
  row_miss <- rowMeans(is.na(dat[, ..met_cols]))
  col_miss <- colMeans(is.na(dat[, ..met_cols]))
  keep_rows <- which(row_miss <= 0.10)
  keep_cols <- met_cols[col_miss <= 0.10]
  dat       <- dat[keep_rows]
  met_cols  <- keep_cols
    cat("metabolites after:", length(met_cols), "\n")
    cat("individuals after:", nrow(dat), "\n")
   # metabolites after: 251 
   # individuals after: 46907
  # KNN imputation
  met_mat  <- as.matrix(dat[, ..met_cols])
  imp_met  <- impute.knn(met_mat, rowmax = 0.10, colmax = 0.10, k = 10)
  dat[, (met_cols) := as.data.table(imp_met$data)]
}
## ====================== 2) continuous covariate imputation ==
cont_vars <- c("bmi", "waist","sbp_auto_i0_a0", "dbp_auto_i0_a0","bb_TC", "bb_HDL", "bb_TG",
               "bb_GLU","bb_HBA1C", "bb_CRE", "bb_LDL","days_pa_mod", "tdi")
cont_vars <- intersect(cont_vars, names(dat))

impute_reg <- function(dt, y, xvars) {
  if (!y %in% names(dt)) return(dt)
  miss_ind <- paste0(y, "_missing")
  dt[[miss_ind]] <- as.integer(is.na(dt[[y]]))
  if (all(is.na(dt[[y]]))) return(dt)
  fit  <- lm(reformulate(xvars, y), data = dt)
  pred <- predict(fit, newdata = dt[, ..xvars])
  dt[[y]] <- ifelse(is.na(dt[[y]]), pred, dt[[y]])
  dt
}

for (v in cont_vars) {
  dat <- impute_reg(dat, v, basic_covs)
}

## enforce no missing in key covariates (basic + categorical + continuous + metabolites)
key_vars <- unique(c(basic_covs, cat_vars, cont_vars, met_cols))
dat <- dat[complete.cases(dat[, ..key_vars])]
  cat("N after all covariate / metabolite imputation & NA removal:", nrow(dat), "\n")
 # N after all covariate / metabolite imputation & NA removal: 46907  
## ====================== 3) eGFR (CKD-EPI 2021) =========
dat[, cre_mgdl := bb_CRE / 88.4]
dat[, male := as.integer(sex == 1)]
dat[, `:=`(
  kappa = ifelse(male == 1, 0.9, 0.7),
  alpha = ifelse(male == 1, -0.302, -0.241)
)]
dat[, egfr := 142 *
      (pmin(cre_mgdl / kappa, 1)^alpha) *
      (pmax(cre_mgdl / kappa, 1)^(-1.200)) *
      (0.9938^age) *
      ifelse(male == 1, 1, 1.012)]
dat[egfr < 0 | is.infinite(egfr), egfr := NA_real_]

## ====================== 4) dates & baseline diagnoses ===
date_cols <- c(
  "icd10Date_ihd", "icd10Date_hfail", "icd10Date_stroke",
  "icd10Date_pad", "icd10Date_afib",
  "icd10Date_mid_high_ckd", "icd10Date_extreme_ckd",
  "icd10Date_htn", "icd10Date_t2dm",
  "date_attend", "date_death"
)
date_cols <- intersect(date_cols, names(dat))
dat[, (date_cols) := lapply(.SD, as.IDate), .SDcols = date_cols]

follow_end <- as.IDate("2024-12-31")

limit_cols <- setdiff(date_cols, "date_attend")
for (v in limit_cols) {
  dat[get(v) > follow_end, (v) := as.IDate(NA)]
}

b4 <- function(dx, base) !is.na(dx) & (is.na(base) | dx <= base)

dat[, `:=`(
  has_cad = b4(icd10Date_ihd,    date_attend),
  has_as  = b4(icd10Date_stroke, date_attend),
  has_hf  = b4(icd10Date_hfail,  date_attend),
  has_af  = b4(icd10Date_afib,   date_attend),
  has_pad = b4(icd10Date_pad,    date_attend),
  has_ht  = b4(icd10Date_htn,    date_attend),
  has_t2d = b4(icd10Date_t2dm,   date_attend),
  has_ckd_mid_high_icd = b4(icd10Date_mid_high_ckd, date_attend),
  has_ckd_extreme_icd  = b4(icd10Date_extreme_ckd,  date_attend)
)]

## ====================== 5) metabolic status =============
dat[, `:=`(
  hyperTG_135 = bb_TG >= 1.53,
  lowHDL      = (male == 1 & bb_HDL < 1.036) | (male == 0 & bb_HDL < 1.29),
  ht_meas     = (sbp_auto_i0_a0 >= 140 | dbp_auto_i0_a0 >= 90)
)]
dat[, HT_any := (ht_meas | has_ht)]

dat[, `:=`(
  preDM = (bb_HBA1C >= 39 & bb_HBA1C < 48) |
    (is.na(bb_HBA1C) & bb_GLU >= 5.6 & bb_GLU < 7.0),
  T2D   = (has_t2d == TRUE) | (bb_HBA1C >= 48) | (bb_GLU >= 7.0)
)]
dat[, wc_high := (male == 1 & waist >= 102) | (male == 0 & waist >= 88)]
dat[, metab_n := as.integer(wc_high) + as.integer(lowHDL) +
      as.integer(bb_TG >= 1.69) + as.integer(HT_any) +
      as.integer(bb_GLU >= 5.55)]
dat[, MetS := metab_n >= 3]

## ====================== 6) CKD tiers ====================
dat[, CKD_mid_high := ((!is.na(egfr) & egfr < 60  & egfr >= 30) | has_ckd_mid_high_icd)]
dat[, CKD_extreme  := ((!is.na(egfr) & egfr < 30)               | has_ckd_extreme_icd)]

## ====================== 7) SCORE inputs =================
dat[, `:=`(
  prior_cvd = as.integer(has_cad | has_as | has_hf | has_af | has_pad),
  Gender    = ifelse(male == 1, "male", "female"),
  smoker    = as.integer(smoke_status == 2),
  diabetes  = as.integer(T2D == TRUE)
)]

cap <- function(x, lo, hi) pmin(pmax(x, lo), hi)

dat[, Age         := cap(age,            40, 89)]
dat[, systolic.bp := cap(sbp_auto_i0_a0, 90, 200)]
dat[, total.chol  := cap(bb_TC,           2, 10)]
dat[, total.hdl   := cap(bb_HDL,        0.5,  3)]
dat[, eGFR        := cap(egfr,           15, 120)]
dat[, HbA1c_mmol  := bb_HBA1C]

## ====================== 8) SCORE2 / SCORE2-DM / SCORE2-CKD ==
sub_s2 <- dat[prior_cvd == 0,
              .(eid, Age, Gender, smoker, systolic.bp, diabetes,
                total.chol, total.hdl)]
s2 <- as.data.table(
  SCORE2_scores(data = sub_s2, Risk.region = "Low", classify = TRUE)
)[, .(eid, SCORE2_score, SCORE2_strat)]
s2[, score2_any := SCORE2_score / 100]

sub_dm <- dat[prior_cvd == 0 & diabetes == 1 & Age >= 40 & Age <= 69,
              .(eid, Age, Gender, smoker,
                systolic.bp, total.chol, total.hdl,
                diabetes = 1L,
                diabetes.age = Age,
                HbA1c = HbA1c_mmol,
                eGFR = eGFR)]

s2_dm <- if (nrow(sub_dm) > 0) {
  tmp <- copy(sub_dm)
  tmp[, SCORE2_DM := SCORE2_Diabetes(
    Risk.region = "Low",
    Age         = Age,
    Gender      = Gender,
    smoker      = smoker,
    systolic.bp = systolic.bp,
    total.chol  = total.chol,
    total.hdl   = total.hdl,
    diabetes    = 1L,
    diabetes.age = diabetes.age,
    HbA1c       = HbA1c,
    eGFR        = eGFR,
    classify    = FALSE
  ), by = eid]
  tmp[, .(eid, SCORE2_DM)]
} else data.table(eid = integer(), SCORE2_DM = numeric())

sub_ckd <- dat[prior_cvd == 0,
               .(eid, Age, Gender,
                 smoker = smoker,
                 systolic.bp, total.chol, total.hdl,
                 diabetes, eGFR, ACR = NA_real_)]
s2_ckd <- as.data.table(
  SCORE2_CKD_scores(data = sub_ckd, Risk.region = "Low", classify = TRUE)
)[, .(eid, SCORE2_CKD_score, SCORE2_CKD_strat)]

dat <- merge(dat, s2,    by = "eid", all.x = TRUE)
dat <- merge(dat, s2_dm, by = "eid", all.x = TRUE)
dat <- merge(dat, s2_ckd, by = "eid", all.x = TRUE)

gen_high <- (!is.na(dat$SCORE2_strat)     & dat$SCORE2_strat     == "High risk")
dm_high  <- (!is.na(dat$SCORE2_DM)        & dat$SCORE2_DM / 100  >= 0.10)
ckd_high <- (!is.na(dat$SCORE2_CKD_strat) & dat$SCORE2_CKD_strat == "High risk")

dat[, highrisk := as.logical(
  ifelse(is.na(gen_high) & is.na(dm_high) & is.na(ckd_high),
         NA,
         (gen_high %in% TRUE) | (dm_high %in% TRUE) | (ckd_high %in% TRUE))
)]

## ====================== 9) CKM staging (0–4, no 4a/4b) ==
dat[, baseline_cvd := (has_cad | has_as | has_hf | has_af | has_pad)]

dat[, stage := NA_integer_]
dat[baseline_cvd == TRUE, stage := 4L]
dat[ is.na(stage) & (highrisk == TRUE | CKD_extreme == TRUE), stage := 3L]
dat[ is.na(stage) & (hyperTG_135 == TRUE | HT_any == TRUE |
                       MetS == TRUE | T2D == TRUE | CKD_mid_high == TRUE), stage := 2L]
dat[ is.na(stage) & ((bmi >= 25) | wc_high == TRUE | preDM == TRUE), stage := 1L]
dat[ is.na(stage), stage := 0L]

dat[, CKM_stage_baseline := factor(
  stage,
  levels = c(0L, 1L, 2L, 3L, 4L),
  labels = c("Stage 0", "Stage 1", "Stage 2", "Stage 3", "Stage 4")
)]
print(table(dat$CKM_stage_baseline, useNA = "ifany"))

## keep only baseline Stage 0–3
dat <- dat[CKM_stage_baseline %chin% c("Stage 0", "Stage 1", "Stage 2", "Stage 3")]
dat[, CKM_stage_baseline := factor(CKM_stage_baseline,
                                   levels = c("Stage 0", "Stage 1", "Stage 2", "Stage 3"))]
  cat("N after excluding baseline Stage 4:", nrow(dat), "\n")
 # N after excluding baseline Stage 4: 43734  
## ====================== 10) CKM4 and death outcomes =====
dat[, icdDate_cvd := do.call(
  pmin,
  c(.SD, list(na.rm = TRUE))
), .SDcols = c("icd10Date_ihd", "icd10Date_hfail", "icd10Date_stroke",
               "icd10Date_pad", "icd10Date_afib")]
dat[is.infinite(icdDate_cvd), icdDate_cvd := as.IDate(NA)]

# CKM4: first incident CVD after baseline
dat[, icd10Date_ckm4 := icdDate_cvd]
dat[!is.na(icd10Date_ckm4) & icd10Date_ckm4 <= date_attend, icd10Date_ckm4 := as.IDate(NA)]
dat[!is.na(icd10Date_ckm4) & icd10Date_ckm4 > follow_end,   icd10Date_ckm4 := as.IDate(NA)]

# death: from date_death, constrained to (baseline, follow_end]
dat[, icd10Date_death := date_death]
dat[!is.na(icd10Date_death) & icd10Date_death <= date_attend, icd10Date_death := as.IDate(NA)]
dat[!is.na(icd10Date_death) & icd10Date_death >  follow_end,  icd10Date_death := as.IDate(NA)]

dat[, prog_ckm4  := !is.na(icd10Date_ckm4)]
dat[, prog_death := !is.na(icd10Date_death)]

## ====================== 11) final checks & save =========
print(table(dat$CKM_stage_baseline, useNA = "ifany"))
 # Stage 0 Stage 1 Stage 2 Stage 3 
 #   5874    6138   28383    3339 
  cat("CKM4 events:", sum(dat$prog_ckm4,  na.rm = TRUE), "\n")
  cat("Death events:", sum(dat$prog_death, na.rm = TRUE), "\n")
 # CKM4 events: 7666 
 # Death events: 4272  
saveRDS(dat, ckm_path)
