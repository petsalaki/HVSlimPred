library(HVSlimPred)

uniprot_PTM = read.table("data/input/UP000005640_9606_mod_res.bed", sep = "\t", quote = "")

#In the time this study was done, ClinVar was downloaded on March 2021, so the same data will be used for reproducibility. However, the user can choose
# to download the latest version by using the "use_latest_clinvar" function.
tmp = tempfile()
download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2021-03.txt.gz", destfile = tmp)
Clinvar = read.delim(tmp, header = T, sep = "\t")
ready_noAA = readRDS("data/input/qslim_insts_clinvar_no_AA.rds")


# Host viral ClinVar and PTM ----------------------------------------------

HV_final_hits = readRDS("data/output/HV_final_hits.rds")
HV_enriched_doms = readRDS("data/output/HV_enriched_doms_with_emppval.rds")
HV_all_pfam_doms = prepare_pfam_doms_enrich(HV_enriched_doms)

#PTM
HV_mapped_PTM = Map_PTM(slim_hits = HV_final_hits, all_pfam_dom = HV_all_pfam_doms, uniprot_mod_res = uniprot_PTM)
#Clinvar
HV_final_clinvar_path = Map_ClinVar(slim_hits = HV_final_hits, Clinvar_variant_summary = Clinvar, pathogenic_only = T, qslim_noAA_gncoords = ready_noAA)
HV_final_clinvar_path$ClinVar_all = HV_final_clinvar_path$ClinVar_all[!is.na(HV_final_clinvar_path$ClinVar_all$Type),]

remove_cols = intersect(colnames(HV_final_clinvar_path$ClinVar_all), c("RS._dbSNP.", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "MatchPos", "NormIC", "new_score"))
for (i in 1:length(HV_final_clinvar_path)){
  HV_final_clinvar_path[[i]] = HV_final_clinvar_path[[i]][,which(colnames(HV_final_clinvar_path[[i]]) %nin% remove_cols)]
  HV_final_clinvar_path[[i]] = HV_final_clinvar_path[[i]][!duplicated(HV_final_clinvar_path[[i]]),]
}


saveRDS(list("ClinVar_path" = HV_final_clinvar_path, "PTM" = HV_mapped_PTM), "data/output/HV_Clinvar_and_PTM.rds")

# Human only ClinVar and PTM ----------------------------------------------

human_final_hits = readRDS("data/output/human_final_hits.rds")
human_enriched_doms = readRDS("data/output/human_only_enriched_doms.rds")
human_all_pfam_doms = prepare_pfam_doms_enrich(HV_enriched_doms)

#PTM
human_mapped_PTM = Map_PTM(slim_hits = human_final_hits, all_pfam_dom = human_all_pfam_doms, uniprot_mod_res = uniprot_PTM)
#Clinvar
human_final_clinvar_path = Map_ClinVar(slim_hits = human_final_hits, Clinvar_variant_summary = Clinvar, pathogenic_only = T, qslim_noAA_gncoords = ready_noAA)
human_final_clinvar_path$ClinVar_all = human_final_clinvar_path$ClinVar_all[!is.na(human_final_clinvar_path$ClinVar_all$Type),]

remove_cols = intersect(colnames(human_final_clinvar_path$ClinVar_all), c("RS._dbSNP.", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "MatchPos", "NormIC", "new_score"))
for (i in 1:length(human_final_clinvar_path)){
  human_final_clinvar_path[[i]] = human_final_clinvar_path[[i]][,which(colnames(human_final_clinvar_path[[i]]) %nin% remove_cols)]
  human_final_clinvar_path[[i]] = human_final_clinvar_path[[i]][!duplicated(human_final_clinvar_path[[i]]),]
}


saveRDS(list("ClinVar_path" = human_final_clinvar_path, "PTM" = human_mapped_PTM), "data/output/human_Clinvar_and_PTM.rds")
