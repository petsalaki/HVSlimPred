library(HVSlimPred)
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
human_fasta_seqs = readRDS("data/input/Human_prots_fasta.rds")
human_human_intact = read.csv("data/input/human_human_intact.csv.gz")

HV_final_hits = readRDS("data/output/HV_final_hits.rds")
human_final_hits = readRDS("data/output/human_final_hits.rds")

HV_motif_prots = unique(HV_final_hits$uniprot)
human_motif_prots = unique(human_final_hits$uniprot)

# Domain enrichment will be calculated on the domain-carrying proteins that directly interact with identified
# motif carrying proteins. CD-Hit was used to remove redundancy of sequences > 70% similar. For the sake of performance
# CD-hit output is already provided. However, users should be able to call CD-hit from within the function if the executable
# CD-hit binary path is provided.

HV_cdhit_out = readRDS("data/input/HV_cdhit_output.rds")
human_cdhit_out = readRDS("data/input/human_only_cdhit_output.rds")

HV_enriched_doms = Domain_Enrichment(human_fasta = human_fasta_seqs, human_human_int = human_human_intact,
                                     mot_carprots = HV_motif_prots, cdhit_output = HV_cdhit_out,
                                     pfam_doms_list = pfam_doms)
human_enriched_doms = Domain_Enrichment(human_fasta = human_fasta_seqs, human_human_int = human_human_intact,
                                        mot_carprots = human_motif_prots, cdhit_output = human_cdhit_out,
                                        pfam_doms_list = pfam_doms)

saveRDS(HV_enriched_doms, "data/output/HV_enriched_doms.rds")
saveRDS(human_enriched_doms, "data/output/human_only_enriched_doms.rds")
