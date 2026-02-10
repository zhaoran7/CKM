pacman::p_load(data.table,impute,ukbtools)

## paths
all_path      <- "D:/ukb/phe/Rdata/all.plus.rds"
prot_raw_path <- "D:/ukb/phe/rap/prot.tab.gz"
rel_path      <- "D:/ukb/raw/ukb_rel.dat"
withdraw_path <- "D:/ukb/raw/withdrawn_66137_20250818.csv"
out_dir       <- "D:/ckm2/data/phe"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## load data
dat       <- readRDS(all_path)
prot_raw  <- fread(prot_raw_path)
rel       <- fread(rel_path)
withdrawn <- fread(withdraw_path)
setnames(withdrawn, 1L, "eid")
setDT(dat)

## baseline: individuals with raw proteomics
base_ids <- intersect(prot_raw$eid, dat$eid)
dat      <- dat[eid %in% base_ids]
prot_all <- grep("^prot_", names(dat), value = TRUE)
pcols    <- prot_all
  cat("Baseline participants:", nrow(dat), "proteins =", length(pcols), "\n")
 # Baseline participants: 52995 proteins = 2923

## remove withdrawn and related samples
dat <- dat[!eid %in% withdrawn$eid]
  cat("After withdrawn removal:", nrow(dat), "proteins =", length(pcols), "\n")
 # After withdrawn removal: 52995 proteins = 2923

rel_rm <- ukbtools::ukb_gen_samples_to_remove(data = rel, ukb_with_data = dat$eid, cutoff = 0.0882)
dat <- dat[!eid %in% rel_rm]
  cat("After relatedness removal:", nrow(dat), "proteins =", length(pcols), "\n")
 # After relatedness removal: 52498 proteins = 2923 

## missingness filtering: proteins >10% and  individuals>10%
col_miss <- dat[, sapply(.SD, function(x) mean(is.na(x))), .SDcols = pcols]
keep_p   <- pcols[col_miss <= 0.10]
dat1     <- dat[, c("eid", keep_p), with = FALSE]
pcols    <- keep_p
  cat("After protein missingness filter:", nrow(dat1), "proteins =", length(pcols), "\n")
 # After protein missingness filter: 52498 proteins = 1459  

row_miss <- dat1[, rowMeans(is.na(.SD)), .SDcols = pcols]
keep_row <- row_miss <= 0.10
dat2     <- dat1[keep_row]
  cat("After individual missingness filter:", nrow(dat2), "proteins =", length(pcols), "\n")
 # After individual missingness filter: 48175 proteins = 1459

## KNN imputation (k = 10)
mat <- as.matrix(dat2[, -1, with = FALSE])
rownames(mat) <- dat2$eid
imp     <- impute::impute.knn(t(mat), k = 10)
mat_imp <- t(imp$data)
dt      <- data.table(eid = dat2$eid, mat_imp)
pnames  <- setdiff(names(dt), "eid")

## RINT + Z-standardization
rint <- function(x) {qnorm((rank(x, ties.method = "average") - 0.5) / length(x))}

dt[, (pnames) := lapply(.SD, rint), .SDcols = pnames]
dt[, (pnames) := lapply(.SD, scale), .SDcols = pnames]

## merge processed proteins back to phenotype data
dat_all     <- dat[eid %in% dt$eid]
dat_phe  <- dat_all[, setdiff(names(dat_all), prot_all), with = FALSE]
ckm          <- merge(dat_phe, dt, by = "eid")
  cat("Final CKM data:", nrow(ckm),"proteins =", sum(grepl("^prot_", names(ckm))), "\n")
 # Final CKM data: 48175 proteins = 1459 

## save
saveRDS(ckm, file.path(out_dir, "ckm.rds"))
