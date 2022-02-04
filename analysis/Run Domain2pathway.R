library(HVSlimPred)
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
pfam_doms = pfam_doms %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")

InterPro_PFAM = read.delim("data/input/Domain2pathway/Interpro_pfam.tsv", na.strings = c("", NA))
InterPro_PFAM = InterPro_PFAM[!is.na(InterPro_PFAM$Integrated.Into),]

tmp = tempfile()
x = tempdir()
download.file(url = "https://netbiolab.org/wiki/images/6/68/Ps.zip", destfile = tmp)
unzip(zipfile = tmp, exdir = x)

real_dom_profile = read.delim(file.path(x, "ps/real_domainprofile.txt"), header = F)
real_Interpro_doms = as.character(real_dom_profile[1,-1])
length(intersect(real_Interpro_doms, InterPro_PFAM$Integrated.Into))

InterPro_PFAM = InterPro_PFAM[which(InterPro_PFAM$Integrated.Into %in% real_Interpro_doms),]
PS_KEGG_interpro = read.delim("data/input/Domain2pathway/PS_KEGG_InterPro.txt", header = F)
colnames(PS_KEGG_interpro) = c("Interpro_id","pathway", "PS_score")
PS_KEGG_interpro = PS_KEGG_interpro[which(PS_KEGG_interpro$Interpro_id %in% InterPro_PFAM$Integrated.Into),]
PS_KEGG_interpro = PS_KEGG_interpro[which(PS_KEGG_interpro$PS_score > 0.056),]
PS_KEGG_interpro = dplyr::left_join(PS_KEGG_interpro, InterPro_PFAM[,c(1,5)], by = c("Interpro_id" = "Integrated.Into"))
PS_KEGG_interpro = PS_KEGG_interpro[!duplicated(PS_KEGG_interpro),]

PS_GO_interpro = read.delim("data/input/Domain2pathway/PS_GOBP_real_interpro.txt", header = F)
colnames(PS_GO_interpro) = c("Interpro_id","pathway", "PS_score")
PS_GO_interpro = PS_GO_interpro[which(PS_GO_interpro$Interpro_id %in% InterPro_PFAM$Integrated.Into),]
PS_GO_interpro = PS_GO_interpro[which(PS_GO_interpro$PS_score > 0.056),]
PS_GO_interpro = dplyr::left_join(PS_GO_interpro, InterPro_PFAM[,c(1,5)], by = c("Interpro_id" = "Integrated.Into"))
PS_GO_interpro = PS_GO_interpro[!duplicated(PS_GO_interpro),]

All_filters_dcGOR = readRDS("data/output/plot_objects/HV_all_filters_dcGOR.rds")

for (i in 1:length(All_filters_dcGOR)){
  for (j in 1:2){
    x = PS_GO_interpro[which(PS_GO_interpro$Accession %in% All_filters_dcGOR[[i]][[j]][["query"]]),]
    y = PS_KEGG_interpro[which(PS_KEGG_interpro$Accession %in% All_filters_dcGOR[[i]][[j]][["query"]]),]
    All_filters_dcGOR[[i]][[j]][["PS"]] = list("GO" = x, "KEGG" = y)
  }
}

PS_info_all_filters_KEGG = list()
PS_info_all_filters_GO = list()
counter = 1
for (i in 1:length(All_filters_dcGOR)){
  for (j in 1:2){
    if ("PS" %in% names(All_filters_dcGOR[[i]][[j]])){
      if (nrow (All_filters_dcGOR[[i]][[j]][["PS"]][["GO"]]) != 0){
        res = All_filters_dcGOR[[i]][[j]][["PS"]][["GO"]][,-1]
        res$dom_prot_type = names(All_filters_dcGOR[[i]])[j]
        res$filter = names(All_filters_dcGOR)[i]
        PS_info_all_filters_GO[[counter]] = res
      }
      else{
        PS_info_all_filters_GO[[counter]] = NULL
      }
      if (nrow (All_filters_dcGOR[[i]][[j]][["PS"]][["KEGG"]]) != 0){
        res = All_filters_dcGOR[[i]][[j]][["PS"]][["KEGG"]][,-1]
        res$dom_prot_type = names(All_filters_dcGOR[[i]])[j]
        res$filter = names(All_filters_dcGOR)[i]
        PS_info_all_filters_KEGG[[counter]] = res
      }
      else{
        PS_info_all_filters_KEGG[[counter]] = NULL
      }
    }
    counter = counter + 1
  }
}
PS_info_all_filters_GO = dplyr::bind_rows(PS_info_all_filters_GO)
PS_info_all_filters_GO$PS_z.score = scale(PS_info_all_filters_GO$PS_score)
PS_info_all_filters_KEGG = dplyr::bind_rows(PS_info_all_filters_KEGG)
PS_info_all_filters_KEGG$PS_z.score = scale(PS_info_all_filters_KEGG$PS_score)

saveRDS(PS_info_all_filters_GO, "data/output/plot_objects/HV_PS_all_filters_GO.rds")
saveRDS(PS_info_all_filters_KEGG, "data/output/plot_objects/HV_PS_all_filters_KEGG.rds")
