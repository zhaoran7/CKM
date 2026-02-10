pacman::p_load(data.table)

root <- "/work/sph-zhaor/analysis/ckm"
feat <- readLines(file.path(root,"data/features.txt"))
tmp  <- file.path(root,"res/tmp_fg")
resd <- file.path(root,"res"); dir.create(resd, showWarnings=FALSE, recursive=TRUE)

need_wait <- function(dir, n){
  repeat{
    k <- length(list.files(dir, pattern="\\.csv$", full.names=TRUE))
    if(k >= n) break
    Sys.sleep(60)
  }
}

read_bind <- function(d){
  ff <- list.files(d, pattern="\\.csv$", full.names=TRUE)
  if(!length(ff)) return(data.table())
  rbindlist(lapply(ff, fread, showProgress=FALSE), fill=TRUE, use.names=TRUE)
}

split_write <- function(dt, tag){
  if(!nrow(dt)) return(invisible(NULL))
  dt[, p_bonf := p.adjust(p, "bonferroni")]
  prot <- dt[grepl("^prot_", protein)]
  met  <- dt[grepl("^met_",  protein)]
  if(nrow(prot)) fwrite(prot, file.path(resd, paste0(tag,"__prot.csv")))
  if(nrow(met))  fwrite(met,  file.path(resd, paste0(tag,"__met.csv")))
}

need_wait(file.path(tmp,"main_ckm4"), length(feat))
need_wait(file.path(tmp,"main_death"), length(feat))
need_wait(file.path(tmp,"sub"),       length(feat))
need_wait(file.path(tmp,"stage"),     length(feat))
need_wait(file.path(tmp,"sex"),       length(feat))

m1 <- read_bind(file.path(tmp,"main_ckm4"))
m2 <- read_bind(file.path(tmp,"main_death"))
s1 <- read_bind(file.path(tmp,"sub"))
g1 <- read_bind(file.path(tmp,"stage"))
x1 <- read_bind(file.path(tmp,"sex"))

keep <- c("group","outcome","protein","N","event1","event2","beta","se","sHR","Lower","Upper","p")
fix <- function(dt){
  miss <- setdiff(keep, names(dt))
  if(length(miss)) dt[,(miss):=NA]
  setcolorder(dt, keep)
  dt
}

split_write(fix(m1), "FG_main_ckm4")
split_write(fix(m2), "FG_main_death")
split_write(fix(s1), "FG_sub_outcomes")
split_write(fix(g1), "FG_stage_stratified")
split_write(fix(x1), "FG_sex_stratified")
