pacman::p_load(data.table)

root <- "/work/sph-zhaor/analysis/ckm"
ckm_path <- file.path(root,"data/ckm.rds")

dir.create(file.path(root,"data/fg_prot"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root,"data/fg_met"),  showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root,"jobs/r"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root,"jobs/sh"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root,"jobs/log"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root,"res/tmp_fg"), showWarnings=FALSE, recursive=TRUE)

dat <- as.data.table(readRDS(ckm_path))

covs <- intersect(c("age","sex","bmi","waist","ethnic.c","tdi","bb_GLU","sbp_auto_i0_a0","dbp_auto_i0_a0",
                    "bb_LDL","bb_TC","bb_CRE","bb_TG","bb_HDL","bb_HBA1C","smoke_status","drug.lipid",
                    "drug.htn","drug.dm","days_pa_mod",paste0("PC",1:10)), names(dat))
p_cols <- grep("^prot_", names(dat), value=TRUE)
m_cols <- grep("^met_",  names(dat), value=TRUE)

cat_covs <- intersect(c("sex","ethnic.c","smoke_status","drug.lipid","drug.htn","drug.dm"), covs)
if(length(cat_covs)) dat[,(cat_covs):=lapply(.SD,as.factor), .SDcols=cat_covs]

need <- intersect(c("date_attend","icd10Date_ckm4","icd10Date_death","icd10Date_ihd","icd10Date_pad",
                    "icd10Date_hfail","icd10Date_afib","icd10Date_stroke","icd10Date_cvd_death",
                    "CKM_stage_baseline"), names(dat))

base <- dat[, c(need, covs), with=FALSE]
saveRDS(base, file.path(root,"data/fg_base.rds"))

writeLines(p_cols, file.path(root,"data/proteins.txt"))
writeLines(m_cols, file.path(root,"data/metabolites.txt"))
feat_cols <- c(p_cols, m_cols)
writeLines(feat_cols, file.path(root,"data/features.txt"))

for(p in p_cols) saveRDS(as.numeric(dat[[p]]), file.path(root,"data/fg_prot", paste0(p,".rds")))
for(m in m_cols) saveRDS(as.numeric(dat[[m]]), file.path(root,"data/fg_met",  paste0(m,".rds")))

worker <- file.path(root,"code/fg_worker.R")
mergeR <- file.path(root,"code/fg_merge.R")
follow_end <- "2023-03-31"

for(i in seq_along(feat_cols)){
  f <- feat_cols[i]
  rfile <- file.path(root,"jobs/r", sprintf("fg_%04d.R", i))
  shfile <- file.path(root,"jobs/sh", sprintf("fg_%04d.sh", i))
  logo <- file.path(root,"jobs/log", sprintf("fg_%04d.out", i))
  loge <- file.path(root,"jobs/log", sprintf("fg_%04d.err", i))

  writeLines(c(
    sprintf('source("%s")', worker),
    sprintf('run_fg_one("%s","%s",as.Date("%s"))', f, root, follow_end)
  ), rfile)

  writeLines(c(
    "#!/usr/bin/env bash",
    sprintf("#BSUB -J FG_%04d", i),
    "#BSUB -q short",
    "#BSUB -n 40",
    sprintf("#BSUB -o %s", logo),
    sprintf("#BSUB -e %s", loge),
    "source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh",
    "conda activate r",
    sprintf("Rscript %s", rfile)
  ), shfile)
}

submit <- file.path(root,"jobs/submit_all.sh")
writeLines(c(
  "#!/usr/bin/env bash",
  "source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh",
  "conda activate r",
  "set -euo pipefail",
  sprintf("ROOT=%s", root),
  "cd ${ROOT}",
  "for f in ${ROOT}/jobs/sh/fg_*.sh; do bsub < ${f}; done",
  sprintf('bsub -J FG_MERGE -q short -n 40 -o %s -e %s "source /work/sph-zhaor/miniconda3/etc/profile.d/conda.sh; conda activate r; Rscript %s"',
          file.path(root,"jobs/log","merge.out"), file.path(root,"jobs/log","merge.err"), mergeR)
), submit)
Sys.chmod(submit, "0755")
