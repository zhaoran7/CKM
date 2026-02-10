pacman::p_load(data.table, cmprsk)

mk_ts <- function(dt, ev, cp, ed){
  da  <- as.Date(dt$date_attend)
  evd <- as.Date(dt[[ev]])
  cpd <- as.Date(dt[[cp]])
  f <- pmin(evd, cpd, ed, na.rm=TRUE)
  t <- as.numeric(f - da)/365.25
  ev_ok <- !is.na(evd) & evd>=da & evd<=ed
  cp_ok <- !is.na(cpd) & cpd>=da & cpd<=ed
  ev1 <- ev_ok & (!cp_ok | evd<=cpd)
  cp2 <- cp_ok & (!ev_ok | cpd<evd)
  list(t=t, s=ifelse(ev1,1L, ifelse(cp2,2L,0L)))
}

fit_fg <- function(dt, t, s, covs, feat, grp, out){
  ok <- is.finite(t) & t>0 & is.finite(dt$fx) & complete.cases(dt[, ..covs])
  if(!any(ok)) return(data.table(group=grp,outcome=out,protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_))
  X <- model.matrix(~., as.data.frame(dt[ok, ..covs])); X <- X[,colnames(X)!="(Intercept)",drop=FALSE]
  z <- as.numeric(scale(dt$fx[ok])); keep <- is.finite(z)
  ft <- t[ok][keep]; fs <- s[ok][keep]
  if(!length(ft)) return(data.table(group=grp,outcome=out,protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_))
  cov1 <- cbind(z=z[keep], X[keep,,drop=FALSE])
  m <- colMeans(cov1); v <- colMeans(cov1^2) - m^2
  cov1 <- cov1[, is.finite(v) & v>0, drop=FALSE]
  if(!ncol(cov1)) return(data.table(group=grp,outcome=out,protein=feat,N=length(ft),event1=sum(fs==1L),event2=sum(fs==2L),beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_))
  q <- qr(cov1)
  cov1 <- cov1[, q$pivot[seq_len(q$rank)], drop=FALSE]
  fit <- try(cmprsk::crr(ft, fs, cov1, failcode=1, cencode=0), silent=TRUE)
  if(inherits(fit,"try-error")) return(data.table(group=grp,outcome=out,protein=feat,N=length(ft),event1=sum(fs==1L),event2=sum(fs==2L),beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_))
  b <- unname(fit$coef[1]); v1 <- fit$var; se <- sqrt(if(is.null(dim(v1))) v1[1] else v1[1,1])
  p <- 2*pnorm(-abs(b/se))
  data.table(group=grp,outcome=out,protein=feat,N=length(ft),event1=sum(fs==1L),event2=sum(fs==2L),
             beta=b,se=se,sHR=exp(b),Lower=exp(b-1.96*se),Upper=exp(b+1.96*se),p=p)
}

run_fg_one <- function(feat, root, ed){
  base <- as.data.table(readRDS(file.path(root,"data/fg_base.rds")))
  fx_path <- if(grepl("^prot_", feat)) file.path(root,"data/fg_prot",paste0(feat,".rds")) else file.path(root,"data/fg_met",paste0(feat,".rds"))
  base[, fx:=as.numeric(readRDS(fx_path))]

  dc <- intersect(c("date_attend","icd10Date_ckm4","icd10Date_death","icd10Date_ihd","icd10Date_pad",
                    "icd10Date_hfail","icd10Date_afib","icd10Date_stroke","icd10Date_cvd_death"), names(base))
  if(length(dc)) base[,(dc):=lapply(.SD, as.Date), .SDcols=dc]

  covs <- intersect(c("age","sex","bmi","waist","ethnic.c","tdi","bb_GLU","sbp_auto_i0_a0","dbp_auto_i0_a0",
                      "bb_LDL","bb_TC","bb_CRE","bb_TG","bb_HDL","bb_HBA1C","smoke_status","drug.lipid",
                      "drug.htn","drug.dm","days_pa_mod",paste0("PC",1:10)), names(base))

  if("icd10Date_cvd_death" %in% names(base)){
    base[, noncvd_death := as.Date(NA)]
    base[!is.na(icd10Date_death) & is.na(icd10Date_cvd_death), noncvd_death := as.Date(icd10Date_death)]
  }

  tmp <- file.path(root,"res/tmp_fg")
  dir.create(file.path(tmp,"main_ckm4"), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(tmp,"main_death"), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(tmp,"sub"), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(tmp,"stage"), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(tmp,"sex"), showWarnings=FALSE, recursive=TRUE)

  ts1 <- mk_ts(base,"icd10Date_ckm4","icd10Date_death",ed)
  ts2 <- mk_ts(base,"icd10Date_death","icd10Date_ckm4",ed)
  fwrite(fit_fg(base, ts1$t, ts1$s, covs, feat, "All", "FG_ckm4_comp_death"), file.path(tmp,"main_ckm4", paste0(feat,".csv")))
  fwrite(fit_fg(base, ts2$t, ts2$s, covs, feat, "All", "FG_death_comp_ckm4"), file.path(tmp,"main_death", paste0(feat,".csv")))

  subs <- c(ihd="icd10Date_ihd", pad="icd10Date_pad", hfail="icd10Date_hfail", afib="icd10Date_afib", stroke="icd10Date_stroke")
  sub_res <- rbindlist(lapply(names(subs), function(nm){
    ev <- subs[[nm]]; if(!ev %in% names(base)) return(NULL)
    ts <- mk_ts(base, ev, "icd10Date_death", ed)
    fit_fg(base, ts$t, ts$s, covs, feat, "All", paste0("FG_",nm,"_comp_death"))
  }), fill=TRUE)
  if(all(c("icd10Date_cvd_death","noncvd_death") %in% names(base))){
    ts <- mk_ts(base,"icd10Date_cvd_death","noncvd_death",ed)
    sub_res <- rbind(sub_res, fit_fg(base, ts$t, ts$s, covs, feat, "All", "FG_cvd_death_comp_noncvd"), fill=TRUE)
  }
  if(!nrow(sub_res)) sub_res <- data.table(group="All",outcome=NA_character_,protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_)
  fwrite(sub_res, file.path(tmp,"sub", paste0(feat,".csv")))

  st_levels <- c("Stage 0","Stage 1","Stage 2","Stage 3")
  st_res <- rbindlist(lapply(st_levels, function(k){
    dtk <- base[CKM_stage_baseline==k]
    if(!nrow(dtk)) return(rbind(
      data.table(group=k,outcome="FG_ckm4_comp_death",protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_),
      data.table(group=k,outcome="FG_death_comp_ckm4",protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_)
    ))
    ts1 <- mk_ts(dtk,"icd10Date_ckm4","icd10Date_death",ed)
    ts2 <- mk_ts(dtk,"icd10Date_death","icd10Date_ckm4",ed)
    rbind(
      fit_fg(dtk, ts1$t, ts1$s, covs, feat, k, "FG_ckm4_comp_death"),
      fit_fg(dtk, ts2$t, ts2$s, covs, feat, k, "FG_death_comp_ckm4")
    )
  }), fill=TRUE)
  fwrite(st_res, file.path(tmp,"stage", paste0(feat,".csv")))

  sx <- as.character(base$sex)
  sx_res <- rbindlist(lapply(list(Male="1", Female="0"), function(v){
    dts <- base[sx %in% v]
    if(!nrow(dts)) return(rbind(
      data.table(group=names(v),outcome="FG_ckm4_comp_death",protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_),
      data.table(group=names(v),outcome="FG_death_comp_ckm4",protein=feat,N=0L,event1=0L,event2=0L,beta=NA_real_,se=NA_real_,sHR=NA_real_,Lower=NA_real_,Upper=NA_real_,p=NA_real_)
    ))
    ts1 <- mk_ts(dts,"icd10Date_ckm4","icd10Date_death",ed)
    ts2 <- mk_ts(dts,"icd10Date_death","icd10Date_ckm4",ed)
    rbind(
      fit_fg(dts, ts1$t, ts1$s, covs, feat, names(v), "FG_ckm4_comp_death"),
      fit_fg(dts, ts2$t, ts2$s, covs, feat, names(v), "FG_death_comp_ckm4")
    )
  }), fill=TRUE)
  fwrite(sx_res, file.path(tmp,"sex", paste0(feat,".csv")))
}
