pacman::p_load(data.table,dplyr,stringr,clusterProfiler,org.Hs.eg.db,ReactomePA,ggplot2,patchwork,tibble,openxlsx)
options(download.file.method="wininet")
dir   <- "D:/ckm2/res/assoc/FG_sens"
typ   <- "prot"
outs  <- c("ckm4","death","ckm4_death")
p_thr <- 0.05/(15*1459)

ann <- fread("D:/ckm2/res/protein_annotations.csv")
feat_col <- if("feature"%chin%names(ann))"feature" else if("Feature"%chin%names(ann))"Feature" else stop("No feature col")
gene_col <- if("Gene name"%chin%names(ann))"Gene name" else if("Gene_name"%chin%names(ann))"Gene_name" else stop("No gene col")
abbr_col <- if("Protein (abbreviation)"%chin%names(ann))"Protein (abbreviation)" else if("Protein_abbreviation"%chin%names(ann))"Protein_abbreviation" else stop("No abbr col")

ann2 <- ann[,.(feature=tolower(as.character(get(feat_col))),
               gene=as.character(get(gene_col)),
               abbr=tolower(as.character(get(abbr_col))))]

map_tbl <- rbindlist(list(
  ann2[,.(key=feature,gene)],
  ann2[!is.na(abbr)&abbr!="",.(key=abbr,gene)],
  ann2[!is.na(gene)&gene!="",.(key=tolower(gene),gene)]
), use.names=TRUE, fill=TRUE)
map_tbl <- unique(map_tbl[!is.na(key)&key!=""&!is.na(gene)&gene!=""])

extract_top <- function(er, tag, n=5){
  if(is.null(er) || nrow(er@result)==0) return(NULL)
  as_tibble(er@result) |>
    arrange(p.adjust) |>
    slice_head(n=n) |>
    transmute(Term=Description,p.adjust,Count=as.numeric(Count),Ontology=tag,negLogFDR=-log10(p.adjust))
}

run_enrich_plot <- function(o, group_use="All"){
  z <- fread(file.path(dir, sprintf("FG_main_%s__%s.csv", o, typ)))
  z[, p := as.numeric(sub("^<","",as.character(p)))]
  z <- z[group==group_use & is.finite(p)]
  z[, key := tolower(as.character(protein))]
  z <- merge(z, map_tbl, by="key", all.x=TRUE)
  
  sig <- z[p < p_thr & !is.na(gene) & gene!="", unique(gene)]
  syms <- unique(unlist(str_split(paste(sig,collapse=";"),"[,;/\\s]+")))
  syms <- syms[syms!=""]
  if(!length(syms)) return(ggplot()+theme_void()+ggtitle(o))
  
  mp <- bitr(syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if(is.null(mp) || !nrow(mp)) return(ggplot()+theme_void()+ggtitle(o))
  
  eg_bp <- enrichGO(mp$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  eg_mf <- enrichGO(mp$ENTREZID, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  eg_cc <- enrichGO(mp$ENTREZID, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  ekg   <- enrichKEGG(mp$ENTREZID, organism="hsa", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1)
  erc   <- enrichPathway(gene=mp$ENTREZID, organism="human", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  
  df <- bind_rows(
    extract_top(eg_bp,"GO BP",5),
    extract_top(eg_mf,"GO MF",5),
    extract_top(eg_cc,"GO CC",5),
    {x<-extract_top(ekg,"KEGG",5); if(!is.null(x)) x$Term<-gsub(" - Homo sapiens \\(human\\)","",x$Term); x},
    extract_top(erc,"Reactome",5)
  )
  if(is.null(df) || !nrow(df)) return(ggplot()+theme_void()+ggtitle(o))
  
  df <- df |>
    mutate(Term=str_wrap(Term,40),
           Ontology=factor(Ontology,levels=c("GO BP","GO MF","GO CC","KEGG","Reactome"))) |>
    group_by(Ontology) |>
    mutate(k=row_number()) |>
    ungroup() |>
    mutate(x_axis=paste(Ontology,k,sep="_"))
  
  ggplot(df,aes(x=x_axis,y=negLogFDR,fill=Ontology))+
    geom_col(width=0.7,linewidth=0.1)+
    geom_text(aes(label=Count),vjust=-0.3,size=2.6)+
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="grey40",linewidth=0.1)+
    scale_x_discrete(labels=df$Term)+
    scale_fill_manual(values=c("GO BP"="#3c78ba","GO MF"="#64b5cd","GO CC"="#d05b69","KEGG"="#57a27d","Reactome"="#d8a74a"),drop=FALSE)+
    labs(title=o,x=NULL,y=expression(-log[10](italic(P)[FDR])))+
    theme_classic(base_size=8)+
    theme(axis.line=element_line(linewidth=0.1,color="black"),
          axis.ticks=element_line(linewidth=0.1,color="black"),
          axis.ticks.length=unit(2,"pt"),
          axis.text.x=element_text(angle=45,hjust=1,size=7,color="black"),
          axis.text.y=element_text(size=7,color="black"),
          legend.position="top",
          legend.title=element_blank(),
          plot.title=element_text(hjust=0.5,face="italic",size=8))
}

p1 <- run_enrich_plot("ckm4")
p2 <- run_enrich_plot("death")
p3 <- run_enrich_plot("ckm4_death")
(p1 / p2 / p3)
topptx(filename = "d:/ckm2/res/plot/enrich.pptx",height=10,width=10)

export_enrich_xlsx <- function(o, group_use="All", out_xlsx=NULL){
  z <- fread(file.path(dir, sprintf("FG_main_%s__%s.csv", o, typ)))
  z[, p := as.numeric(sub("^<","",as.character(p)))]
  z <- z[group==group_use & is.finite(p)]
  z[, key := tolower(as.character(protein))]
  z <- merge(z, map_tbl, by="key", all.x=TRUE)
  
  sig <- z[p < p_thr & !is.na(gene) & gene!="", unique(gene)]
  syms <- unique(unlist(str_split(paste(sig,collapse=";"),"[,;/\\s]+")))
  syms <- syms[syms!=""]
  if(!length(syms)) return(NULL)
  
  mp <- bitr(syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if(is.null(mp) || !nrow(mp)) return(NULL)
  
  eg_bp <- enrichGO(mp$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  eg_mf <- enrichGO(mp$ENTREZID, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  eg_cc <- enrichGO(mp$ENTREZID, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  ekg   <- enrichKEGG(mp$ENTREZID, organism="hsa", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1)
  erc   <- enrichPathway(gene=mp$ENTREZID, organism="human", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
  
  to_dt <- function(er, tag){
    if(is.null(er) || nrow(er@result)==0) return(NULL)
    as.data.table(er@result)[, Ontology := tag]
  }
  
  dt <- rbindlist(list(
    to_dt(eg_bp,"GO BP"),
    to_dt(eg_mf,"GO MF"),
    to_dt(eg_cc,"GO CC"),
    {x<-to_dt(ekg,"KEGG"); if(!is.null(x)) x[, Description := gsub(" - Homo sapiens \\(human\\)","",Description)]; x},
    to_dt(erc,"Reactome")
  ), use.names=TRUE, fill=TRUE)
  
  if(is.null(dt) || !nrow(dt)) return(NULL)
  setcolorder(dt, c("Ontology", setdiff(names(dt), "Ontology")))
  dt
}

out_xlsx <- file.path(dir, sprintf("enrichment_%s_%s.xlsx", typ, format(Sys.Date(), "%Y%m%d")))
wb <- createWorkbook()

for(o in outs){
  dt <- export_enrich_xlsx(o)
  addWorksheet(wb, sheetName = o)
  if(is.null(dt)){
    writeData(wb, o, data.table(note="No significant genes or no enrichment results."))
  } else {
    writeData(wb, o, dt)
  }
}

saveWorkbook(wb, out_xlsx, overwrite=TRUE)
cat("Saved:", out_xlsx, "\n")