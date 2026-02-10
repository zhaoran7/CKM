## ============================================================
## CKM Clinical / Protein / Metabolite / Combined Models
## Elastic Net Cox, 100 resamples (7:3), 2 outcomes (CKM4 / Death)
## ============================================================
pacman::p_load(data.table, dplyr, lubridate, glmnet, survival, timeROC, readr, fs, pROC, readxl)

## ---------------------- paths ---------------------- ##
base_dir  <- "D:/ckm2"
phe_dir   <- file.path(base_dir, "data/phe")
model_dir <- file.path(base_dir, "res/model2")
dir_create(phe_dir)
dir_create(model_dir)

met_annot_path <- file.path(base_dir, "res/metabolite_annotation.xlsx")
met_annot <- as.data.table(readxl::read_xlsx(met_annot_path))
met_dict  <- setNames(as.character(met_annot$Metabolite), as.character(met_annot$`UKB field ID`))

## ------------------- helper funcs ------------------ ##
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))

map_feat <- function(x){
  x <- as.character(x)
  x <- ifelse(grepl("^prot_", x), toupper(sub("^prot_", "", x)), x)
  is_met <- grepl("^met_", x)
  if(any(is_met)){
    id <- sub("^met_", "", x[is_met])
    nm <- unname(met_dict[id])
    x[is_met] <- ifelse(is.na(nm) | nm=="", id, nm)
  }
  x
}

auc_ci <- function(y, s){
  r <- pROC::roc(response=y, predictor=s, levels=c(0,1), direction="<", quiet=TRUE)
  ci <- as.numeric(pROC::ci.auc(r, conf.level=0.95))
  list(roc=r, auc=as.numeric(pROC::auc(r)), l=ci[1], u=ci[3])
}

make_surv_outcome <- function(dat, outcome=c("ckm4","death"), follow_end=as.Date("2023-03-31")){
  outcome <- match.arg(outcome)
  ev_var <- if(outcome=="ckm4") "icd10Date_ckm4" else "icd10Date_death"
  dat %>%
    mutate(
      date_attend = as.Date(date_attend),
      date_death  = if ("date_death" %in% names(dat)) as.Date(date_death) else as.Date(NA),
      event_date  = as.Date(.data[[ev_var]])
    ) %>%
    mutate(
      follow_date = case_when(
        !is.na(event_date) & event_date >= date_attend ~ event_date,
        !is.na(date_death)                             ~ date_death,
        TRUE                                           ~ follow_end
      ),
      t2e  = as.numeric(follow_date - date_attend)/365.25,
      Yt2e = ifelse(!is.na(event_date) & event_date >= date_attend & event_date <= follow_end, 1L, 0L)
    ) %>%
    filter(is.finite(t2e) & t2e > 0) %>%
    {if(outcome=="ckm4") rename(., ckm4.t2e=t2e, ckm4.Yt2e=Yt2e) else rename(., death.t2e=t2e, death.Yt2e=Yt2e)}
}

split_train_test <- function(df, Y, frac=0.7, seed=1){
  set.seed(seed)
  y_col <- paste0(Y, ".Yt2e")
  cases <- df %>% filter(.data[[y_col]]==1L)
  ctrls <- df %>% filter(.data[[y_col]]==0L)
  n1 <- floor(nrow(cases)*frac); n0 <- floor(nrow(ctrls)*frac)
  tr1 <- sample(cases$eid, n1); tr0 <- sample(ctrls$eid, n0)
  train_ids <- c(tr1, tr0); test_ids <- setdiff(df$eid, train_ids)
  list(train_ids=train_ids, test_ids=test_ids)
}

run_enet_cox <- function(Y, X_label, X.covs, train_df, test_df, out_dir, seed, alpha=0.5, folds=10){
  y_t2e  <- paste0(Y, ".t2e"); y_Yt2e <- paste0(Y, ".Yt2e")
  X.train <- data.matrix(train_df[, X.covs, drop=FALSE])
  X.test  <- data.matrix(test_df[,  X.covs, drop=FALSE])
  Y.train <- Surv(train_df[[y_t2e]], train_df[[y_Yt2e]])
  Y.test  <- Surv(test_df[[y_t2e]],  test_df[[y_Yt2e]])
  
  cvfit <- cv.glmnet(X.train, Y.train, family="cox", alpha=alpha, nfolds=folds, type.measure="C")
  best_lambda <- if(length(X.covs) <= 20) cvfit$lambda.min else cvfit$lambda.1se
  
  beta_mat <- as.matrix(coef(cvfit, s=best_lambda))
  feat_raw <- rownames(beta_mat)
  feat_out <- map_feat(feat_raw)
  features <- data.frame(feature=feat_out, beta.net.cox=as.numeric(beta_mat), stringsAsFactors=FALSE)
  
  pred_test <- as.numeric(predict(cvfit, newx=X.test, s=best_lambda, type="link"))
  risk_col  <- sprintf("risk.%s.%s", "net", "cox")
  risk <- data.frame(eid=test_df$eid, t2e=test_df[[y_t2e]], Yt2e=test_df[[y_Yt2e]], stringsAsFactors=FALSE)
  risk[[risk_col]] <- pred_test
  
  c_ind <- survival::concordance(Y.test ~ pred_test, reverse=TRUE)$concordance
  roc0  <- auc_ci(risk$Yt2e, risk[[risk_col]])
  
  pars <- data.frame(
    seed=seed, outcome=Y, model=X_label, alpha=alpha, lambda=best_lambda,
    n_vars=sum(beta_mat!=0), c_index=c_ind, auc=roc0$auc, auc_l=roc0$l, auc_u=roc0$u,
    stringsAsFactors=FALSE
  )
  
  dir_create(out_dir)
  write_csv(pars,     file.path(out_dir, sprintf("%s.%s.seed%03d.0.params.csv",   Y, X_label, seed)))
  write_csv(features, file.path(out_dir, sprintf("%s.%s.seed%03d.1.features.csv", Y, X_label, seed)))
  write_csv(risk,     file.path(out_dir, sprintf("%s.%s.seed%03d.2.risk.csv",     Y, X_label, seed)))
  
  list(pars=pars, pred=pred_test, y=risk$Yt2e, feat_raw=feat_raw, beta=as.numeric(beta_mat))
}

## ============================================================
## 1. Read CKM data, define covariate sets
## ============================================================
dat <- readRDS(file.path(phe_dir, "ckm.rds"))
dat <- as.data.table(dat)

follow_end <- as.Date("2023-03-31")

## clinical variables
clin_vars <- c(
  "age","sex","bmi","waist","ethnic.c","tdi","bb_GLU",
  "sbp_auto_i0_a0","dbp_auto_i0_a0",
  "bb_LDL","bb_TC","bb_CRE","bb_TG","bb_HDL","bb_HBA1C",
  "smoke_status","drug.lipid","drug.htn","drug.dm",
  "days_pa_mod"
)
clin_vars <- intersect(clin_vars, names(dat))

prot_vars <- grep("^prot_", names(dat), value=TRUE)
met_vars  <- grep("^met_",  names(dat), value=TRUE)

clin_cat_vars  <- intersect(c("sex","ethnic.c","smoke_status","drug.lipid","drug.htn","drug.dm"), clin_vars)
clin_cont_vars <- setdiff(clin_vars, clin_cat_vars)

for(v in clin_cat_vars) if(v %in% names(dat) && !is.factor(dat[[v]])) dat[[v]] <- factor(dat[[v]])
lvl_list <- lapply(clin_cat_vars, function(v) levels(dat[[v]])); names(lvl_list) <- clin_cat_vars

build_clin_matrix <- function(df){
  Xc <- if(length(clin_cont_vars)>0) scale(as.data.frame(df[, ..clin_cont_vars])) else NULL
  Xd <- NULL
  if(length(clin_cat_vars)>0){
    for(v in clin_cat_vars) df[[v]] <- factor(df[[v]], levels=lvl_list[[v]])
    Xd <- model.matrix(~ . - 1, data=df[, ..clin_cat_vars])
  }
  if(!is.null(Xc) && !is.null(Xd)) cbind(Xc, Xd) else if(!is.null(Xc)) Xc else Xd
}

## ============================================================
## 2. Build survival datasets for CKM4 and Death
## ============================================================
## CKM4: restrict to baseline stage 0-3 if ckm.stage exists
dat_ckm4_src <-  dat
dat_ckm4  <- make_surv_outcome(as.data.frame(dat_ckm4_src), outcome="ckm4",  follow_end=follow_end)
dat_death <- make_surv_outcome(as.data.frame(dat),         outcome="death", follow_end=follow_end)

saveRDS(dat_ckm4,  file.path(phe_dir, "ckm_surv_ckm4_20230331.rds"))
saveRDS(dat_death, file.path(phe_dir, "ckm_surv_death_20230331.rds"))

cat("N (CKM4): ",  nrow(dat_ckm4),  " Events:", sum(dat_ckm4$ckm4.Yt2e),  "\n")
cat("N (Death): ", nrow(dat_death), " Events:", sum(dat_death$death.Yt2e), "\n")

## ============================================================
## 3. 100 resamples (seed 1-100), 4 models × 2 outcomes
## ============================================================
run_one_outcome <- function(dat_surv, Y, out_root){
  dir_create(out_root); dir_create(file.path(out_root,"splits"))
  auc_rows <- list()
  
  ## selection counters (raw names)
  prot_all <- prot_vars; met_all <- met_vars
  prot_cnt_prot <- setNames(integer(length(prot_all)), prot_all)
  prot_cnt_comb <- setNames(integer(length(prot_all)), prot_all)
  met_cnt_met   <- setNames(integer(length(met_all)),  met_all)
  met_cnt_comb  <- setNames(integer(length(met_all)),  met_all)
  
  for(seed in 1:100){
    split0 <- split_train_test(dat_surv, Y=Y, frac=0.7, seed=seed)
    
    write_lines(as.character(split0$train_ids), file.path(out_root,"splits", sprintf("train_eids_%s_seed%03d.txt", Y, seed)))
    write_lines(as.character(split0$test_ids),  file.path(out_root,"splits", sprintf("test_eids_%s_seed%03d.txt",  Y, seed)))
    
    train <- dat_surv %>% filter(eid %in% split0$train_ids)
    test  <- dat_surv %>% filter(eid %in% split0$test_ids)
    
    ## clinical
    X_clin_tr <- build_clin_matrix(as.data.table(train))
    X_clin_te <- build_clin_matrix(as.data.table(test))
    clin_cols <- colnames(X_clin_tr)
    
    df_clin_tr <- cbind(train[, c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_clin_tr)
    df_clin_te <- cbind(test[,  c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_clin_te)
    
    ## proteins / metabolites (keep your original: inormal + scale, train/test separately)
    X_prot_tr <- if(length(prot_vars)>0) scale(as.data.frame(lapply(as.data.table(train)[, prot_vars, with=FALSE], inormal))) else NULL
    X_prot_te <- if(length(prot_vars)>0) scale(as.data.frame(lapply(as.data.table(test)[,  prot_vars, with=FALSE], inormal))) else NULL
    X_met_tr  <- if(length(met_vars)>0)  scale(as.data.frame(lapply(as.data.table(train)[, met_vars,  with=FALSE], inormal))) else NULL
    X_met_te  <- if(length(met_vars)>0)  scale(as.data.frame(lapply(as.data.table(test)[,  met_vars,  with=FALSE], inormal))) else NULL
    
    df_prot_tr <- if(!is.null(X_prot_tr)) cbind(train[, c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_prot_tr) else NULL
    df_prot_te <- if(!is.null(X_prot_te)) cbind(test[,  c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_prot_te) else NULL
    df_met_tr  <- if(!is.null(X_met_tr))  cbind(train[, c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_met_tr)  else NULL
    df_met_te  <- if(!is.null(X_met_te))  cbind(test[,  c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_met_te)  else NULL
    
    df_comb_tr <- cbind(train[, c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_clin_tr, X_prot_tr, X_met_tr)
    df_comb_te <- cbind(test[,  c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e"))], X_clin_te, X_prot_te, X_met_te)
    
    ## run 4 models
    out_dir <- out_root
    
    r_clin <- run_enet_cox(Y, "clin", clin_cols,
                           as.data.frame(df_clin_tr), as.data.frame(df_clin_te),
                           out_dir, seed)
    
    r_prot <- if(!is.null(df_prot_tr)) run_enet_cox(Y, "prot", colnames(X_prot_tr),
                                                    as.data.frame(df_prot_tr), as.data.frame(df_prot_te),
                                                    out_dir, seed) else NULL
    
    r_met  <- if(!is.null(df_met_tr)) run_enet_cox(Y, "met",  colnames(X_met_tr),
                                                   as.data.frame(df_met_tr),  as.data.frame(df_met_te),
                                                   out_dir, seed) else NULL
    
    comb_cols <- colnames(df_comb_tr)[!(colnames(df_comb_tr) %in% c("eid", paste0(Y,".t2e"), paste0(Y,".Yt2e")))]
    r_comb <- run_enet_cox(Y, "comb", comb_cols,
                           as.data.frame(df_comb_tr), as.data.frame(df_comb_te),
                           out_dir, seed)
    
    ## selection count (raw names)
    if(!is.null(r_prot)){
      sel <- r_prot$feat_raw[r_prot$beta!=0]
      if(length(sel)) prot_cnt_prot[sel] <- prot_cnt_prot[sel] + 1L
    }
    if(!is.null(r_met)){
      sel <- r_met$feat_raw[r_met$beta!=0]
      if(length(sel)) met_cnt_met[sel] <- met_cnt_met[sel] + 1L
    }
    selc <- r_comb$feat_raw[r_comb$beta!=0]
    if(length(selc)){
      sp <- selc[grepl("^prot_", selc)]
      sm <- selc[grepl("^met_",  selc)]
      if(length(sp)) prot_cnt_comb[sp] <- prot_cnt_comb[sp] + 1L
      if(length(sm)) met_cnt_comb[sm]  <- met_cnt_comb[sm]  + 1L
    }
    
    ## AUC summary + DeLong vs clin (same test set)
    roc_clin <- auc_ci(r_clin$y, r_clin$pred)
    
    mk_row <- function(tag, r){
      if(is.null(r)) return(list(auc=NA_real_, l=NA_real_, u=NA_real_, d=NA_real_, p=NA_real_))
      rr <- auc_ci(r$y, r$pred)
      dt <- rr$auc - roc_clin$auc
      pv <- as.numeric(pROC::roc.test(roc_clin$roc, rr$roc, method="delong")$p.value)
      list(auc=rr$auc, l=rr$l, u=rr$u, d=dt, p=pv)
    }
    
    p_prot <- mk_row("prot", r_prot)
    p_met  <- mk_row("met",  r_met)
    p_comb <- mk_row("comb", r_comb)
    
    auc_rows[[seed]] <- data.frame(
      seed=seed,
      clin_auc=roc_clin$auc, clin_l=roc_clin$l, clin_u=roc_clin$u,
      prot_auc=p_prot$auc, prot_l=p_prot$l, prot_u=p_prot$u, prot_dAUC=p_prot$d, prot_pDeLong=p_prot$p,
      met_auc =p_met$auc,  met_l =p_met$l,  met_u =p_met$u,  met_dAUC =p_met$d,  met_pDeLong =p_met$p,
      comb_auc=p_comb$auc, comb_l=p_comb$l, comb_u=p_comb$u, comb_dAUC=p_comb$d, comb_pDeLong=p_comb$p,
      stringsAsFactors=FALSE
    )
    
    cat("DONE", Y, "seed", seed, "\n")
  }
  
  auc_df <- do.call(rbind, auc_rows)
  write_csv(auc_df, file.path(out_root, sprintf("%s.auc_100seeds.csv", Y)))
  
  ## selection tables (proteins / metabolites)
  prot_tab <- data.table(
    protein = map_feat(names(prot_cnt_prot)),
    sel_prot = as.integer(prot_cnt_prot),
    sel_comb = as.integer(prot_cnt_comb)
  )
  prot_tab[, sel_total := sel_prot + sel_comb][order(-sel_total, -sel_prot, -sel_comb)]
  fwrite(prot_tab, file.path(out_root, sprintf("%s.protein_selection_counts.csv", Y)))
  
  met_tab <- data.table(
    metabolite = map_feat(names(met_cnt_met)),
    sel_met = as.integer(met_cnt_met),
    sel_comb = as.integer(met_cnt_comb)
  )
  met_tab[, sel_total := sel_met + sel_comb][order(-sel_total, -sel_met, -sel_comb)]
  fwrite(met_tab, file.path(out_root, sprintf("%s.metabolite_selection_counts.csv", Y)))
  
  invisible(TRUE)
}

## ============================================================
## 4. Run (CKM4 / Death) separately
## ============================================================
out_dir_ckm4  <- file.path(model_dir, "ckm4")
out_dir_death <- file.path(model_dir, "death")

run_one_outcome(dat_ckm4,  Y="ckm4",  out_root=out_dir_ckm4)
run_one_outcome(dat_death, Y="death", out_root=out_dir_death)

count_met <- function(out_root, Y){
  mf <- list.files(out_root, pattern=paste0("^",Y,"\\.met\\.seed\\d+\\.1\\.features\\.csv$"),  full.names=TRUE)
  cf <- list.files(out_root, pattern=paste0("^",Y,"\\.comb\\.seed\\d+\\.1\\.features\\.csv$"), full.names=TRUE)
  met_all <- fread(mf[1])$feature
  sm <- table(unlist(lapply(mf, \(f) fread(f)[beta.net.cox!=0, feature])))
  sc <- table(unlist(lapply(cf, \(f) fread(f)[beta.net.cox!=0 & feature %in% met_all, feature])))
  out <- data.table(metabolite=met_all, sel_met=as.integer(sm[met_all]), sel_comb=as.integer(sc[met_all]))
  out[is.na(sel_met), sel_met:=0L][is.na(sel_comb), sel_comb:=0L][, sel_total:=sel_met+sel_comb][order(-sel_total,-sel_met,-sel_comb)]
}

out_root <- "D:/ckm2/res/model2/ckm4";  fwrite(count_met(out_root,"ckm4"),  file.path(out_root,"ckm4.metabolite_selection_counts.csv"))
out_root <- "D:/ckm2/res/model2/death"; fwrite(count_met(out_root,"death"), file.path(out_root,"death.metabolite_selection_counts.csv"))

