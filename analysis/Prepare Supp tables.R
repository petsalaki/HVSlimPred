library(readxl)
library(openxlsx)
hs <- createStyle(
  textDecoration = "BOLD", fontColour = "#000000", fontSize = 12,
  fontName = "Arial", fgFill = "#B8CCE4", halign = "center", valign = "center", border = "LeftRight")

write_supp_table = function(data_list, filename){
  if (class(data_list) != "list"){
    data_list = list(x = data_list)
    names(data_list) = filename
  }
  out_file = file.path("../Manuscript/Supp_tables_after_review", paste0(filename,".xlsx"))
  if (file.exists(out_file)){
    if (length(excel_sheets(out_file)) == 2){
      Col_desc = readxl::read_xlsx(out_file, sheet = 1)
      data_list = c(list("Column Description" = Col_desc), data_list)
      write.xlsx(data_list, paste0("../Manuscript/Supp_tables_after_review/", filename, ".xlsx"),
                 keepNA = T, na.string = "NA", colWidths = "auto", headerStyle = hs, borders = "columns")
    }
    else{
      write.xlsx(data_list, paste0("../Manuscript/Supp_tables_after_review/", filename, ".xlsx"),
                 keepNA = T, na.string = "NA", colWidths = "auto", headerStyle = hs, borders = "columns")
    }
  }
  else{
    write.xlsx(data_list, paste0("../Manuscript/Supp_tables_after_review/", filename, ".xlsx"),
               keepNA = T, na.string = "NA", colWidths = "auto", headerStyle = hs, borders = "columns")
  }
}

out_dir = "../Manuscript/Supp_tables_after_review/"

x1 = readRDS("data/output/HV_final_hits.rds")
z1 = Classify_insts_by_score(x1)
y1 = readRDS("data/output/HV_classified_insts.rds")
y1$Dataset = paste0(y1$V1, "_", y1$H1)
y1 = y1 %>% dplyr::select(-V1, -H1)
y1 = y1[,c(11,1:10)]
y1 = y1[!duplicated(y1),]
x1 = left_join(x1, y1)
x1 = x1 %>% dplyr::select(-Known_ELM) %>% dplyr::left_join(z1)
write_supp_table(x1, filename = "Host_Viral_raw_hits")


x2 = readRDS("data/output/human_final_hits.rds")
z2 = Classify_insts_by_score(x2)
y2 = readRDS("data/output/human_classified_insts.rds")
y2 = y2[!duplicated(y2),]
y2 = y2[!duplicated(y2[,c(1:4)]),] %>% dplyr::select(-H1)
x2 = left_join(x2,y2)
x2 = x2 %>% dplyr::select(-Known_ELM) %>% dplyr::left_join(z2)
write_supp_table(x2, filename = "Human_only_raw_hits")

x1 = readRDS("data/output/HV_Clinvar_and_PTM.rds")
write_supp_table(list("HV_ClinVar_All" = x1$ClinVar_path$ClinVar_all), filename = "HV_ClinVar_All")
write_supp_table(list("HV_ClinVar_Pathogenic" = x1$ClinVar_path$ClinVar_path), filename = "HV_ClinVar_Pathogenic")

x1 = read.csv("data/input/ChEMBL/filtred_chembl_map.csv.gz")
x2 = read.csv("data/input/ChEMBL/Drug_indication_chembl_27.csv")
x2 = x2[which(x2$molregno %in% x1$molregno),]
x2 = dplyr::left_join(x2, x1[,c(1,2)])
x2 = x2 %>% dplyr::select(compound_chembl_id, max_phase_for_ind, efo_id, efo_term)
x2 = x2[!duplicated(x2),]
x1 = x1 %>% dplyr::select(compound_chembl_id, accession)
x1 = x1[!duplicated(x1),]
x3 = readRDS("data/output/HV_H1_chembl.rds")
x3 = x3$H1_and_Clinvar
write_supp_table(list("ChEMBL Drug-target raw" = x1), filename = "curated_Chembl_Drug_Target_pairs")
write_supp_table(list("Drug Indications" = x2), filename = "curated_Chembl_Drug_indications")
write_supp_table(list("HV_Drug_mapping" = x3), filename = "HV_Drug_mapping")


x1 = read.csv("data/input/SuppData1_Host-Vir-ints.csv")
x1 = x1[,-c(10,11)]

x2 = read.csv("data/input/human_human_intact.csv.gz")
x2$no_in_viral = as.character(x2$no_in_viral)
x2$no_in_viral = stringr::str_replace_all(x2$no_in_viral, c("1" = "A", "2" = "B"))

x3.1 = x1 %>% group_by(org_A) %>% dplyr::summarise(n=n())
x3.1 = x3.1[which(x3.1$org_A != "human"),]
x3.2 = x1 %>% group_by(org_B) %>% dplyr::summarise(n=n())
x3.2 = x3.2[which(x3.2$org_B != "human"),]
y = unique(c(x3.1$org_A, x3.2$org_B))
x3 = data.frame(viral_strain = y, No.of_interactions = 0)
for (i in 1:nrow(x3)){
  counter = 0
  if(x3$viral_strain[i] %in% x3.1$org_A){
    counter = counter + x3.1$n[which(x3.1$org_A == x3$viral_strain[i])]
  }
  if(x3$viral_strain[i] %in% x3.2$org_B){
    counter = counter + x3.2$n[which(x3.2$org_B == x3$viral_strain[i])]
  }
  x3$No.of_interactions[i] = counter
}

write_supp_table(list("HV_intact_PPI" = x1), filename = "HV_intact_PPI")
write_supp_table(list("Human_Human_intact_PPI" = x2), filename = "Human_Human_intact_PPI")
write_supp_table(list("Interactions per viral strain" = x3), filename = "Interactions per viral strain")


ELM = read.csv("data/input/ELM/elm_all_instances.csv")
colnames(ELM) = c("Accession","Class_type","Id","Entry_name","uniprot","other_uniprots","start","end","Pubmeds","Method","Logic","PDBs","Organism")
ELM = ELM[which(ELM$Organism == "Homo sapiens"),]
ELM = ELM[stringr::str_detect(ELM$Entry_name, "HUMAN"),]
ELM = ELM[which(ELM$Logic == "true positive"),]

pwm_prmdb_scan = read.csv("data/input/PRMDB/pwm_prmdb_scan.csv.gz")
pwm_ids = unique(pwm_prmdb_scan$Id)
domain_ids_pwm = list()
for (i in 1:length(pwm_ids)){
  if (nchar(pwm_ids[i]) == 8){
    domain_ids_pwm[[i]] = data.frame(Id = pwm_ids[i], domain_id = pwm_ids[i])
  }
  else{
    domain_ids_pwm[[i]] = data.frame(Id = pwm_ids[i],
                                     domain_id = paste(strsplit(pwm_ids[i], "_")[[1]][c(1,2)], collapse = "_"))
  }
}
domain_ids_pwm = dplyr::bind_rows(domain_ids_pwm)
pwm_prmdb_scan = dplyr::left_join(pwm_prmdb_scan, domain_ids_pwm)

PRM_master = read.csv("data/input/PRMDB/PRM_Master.csv")
PRM_master = PRM_master[which(PRM_master$Species == "Human"),]
PRM_master = PRM_master[,c(1:7)]
colnames(PRM_master)[1] = "domain_id"
pwm_prmdb_scan = pwm_prmdb_scan[which(pwm_prmdb_scan$domain_id %in% PRM_master$domain_id),]
names(PRM_master) = c("domain_id", "Domain_protein_name", "Domain_class", "Domain_class_number", "Domain_uniprot","Domain_start", "Domain_end")
PRMdb_data = dplyr::left_join(pwm_prmdb_scan, PRM_master)

write_supp_table(list("ELM" = ELM), filename = "ELM_bg")
write_supp_table(list("PRMdb" = PRMdb_data), filename = "PRMdb_bg")


pfam_doms = readRDS("data/input/pfam_doms_human.rds")
pfam_doms = pfam_doms %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")
dom_id_name = pfam_doms %>% dplyr::select(PFAM_ID, PFAM_name)
dom_id_name = dom_id_name[!duplicated(dom_id_name),]

dom_clusts = readRDS("data/output/GOSemSim_domain_clusts.rds")

dom_clust_annot = list()
for (i in 1:length(dom_clusts$domain_clusters)){
  res = data.frame(Clust = names(dom_clusts$domain_clusters)[i], PFAM_name = dom_clusts$domain_clusters[[i]]$members)
  res$Clust_center = ifelse(res$PFAM_name %in% dom_clusts$domain_clusters[[i]]$min_dist_center, "1","0")
  res = dplyr::left_join(res, dom_id_name)
  dom_clust_annot[[i]] = res
}
dom_clust_annot = dplyr::bind_rows(dom_clust_annot)
dom_clust_annot = dplyr::left_join(dom_clust_annot, dom_clusts$clust_GO_annotation)
dom_clust_annot = dom_clust_annot[,c(2,4,1,3,5,6)]

write_supp_table(list("Domain_enrichment_clusters" = dom_clust_annot), filename = "Domain_enrichment_clusters")


GOBP_clusts = readRDS("data/output/GOSemSim_GOBP_prot_clusts.rds")

GO_all = ontologyIndex::get_OBO("data/input/dcGOR/go-basic.obo")
GO_annot = data.frame(GO_all$name)
GO_annot$GO_ID = rownames(GO_annot)
colnames(GO_annot) = c("GO_name", "GO_ID")


GOBP_clusts_annot = list()
for (i in 1:length(GOBP_clusts$GOBP_clusters)){
  res = data.frame(Clust = names(GOBP_clusts$GOBP_clusters)[i], GO_ID = GOBP_clusts$GOBP_clusters[[i]]$members)
  res$Clust_center = ifelse(res$GO_ID %in% GOBP_clusts$GOBP_clusters[[i]]$min_dist_center, "1","0")
  res = dplyr::left_join(res, GO_annot)
  GOBP_clusts_annot[[i]] = res
}
GOBP_clusts_annot = dplyr::bind_rows(GOBP_clusts_annot)
GOBP_clusts_annot = GOBP_clusts_annot[,c(4,2,1,3)]

write_supp_table(list("GOBP_prot_enrichment_clusters" = GOBP_clusts_annot), filename = "GOBP_prot_enrichment_clusters")

#Count classified insts per filter
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
HV_class_insts = Classify_insts_by_score(HV_final_hits)
all_data = readRDS("data/output/HV_all_filters_data.rds")

final = list()
for (i in 1:nrow(all_data$all_hits$Metadata)){
  res = plyr::match_df(HV_class_insts, all_data$all_hits[[all_data$all_hits$Metadata$sym[i]]])
  if (nrow(res) == 0){next()}
  res = as.data.frame(table(res$inst_type))
  res = t(res)
  colnames(res) = res[1,]
  res = res[-1,]
  res = t(res)
  res = as.data.frame(res)
  res$Filter = all_data$all_hits$Metadata$Term[i]
  res = res[,c(ncol(res),1:(ncol(res) - 1))]
  final[[i]] = res
}
final = dplyr::bind_rows(final)
final = as.matrix(final)
final[is.na(final)] = 0
final = as.data.frame(final)
final$Filter = stringr::str_replace_all(final$Filter, c("ClinVar_Path" = "ClinVar"))
final$Filter = stringr::str_replace_all(final$Filter, c("_" = " ", "&" = ","))

for (i in 1:ncol(final)){
  final[,i] = trimws(final[,i])
}
write.csv(final, "../Manuscript/Supp_tables_after_review/Classification_per_filter.csv", row.names = F)
