#!/bin/bash
for file in /home/bwadie/Hmmer_H1/pfam_hmms/*.hmm ; do
	filename=$(basename -- "$file")
        file_name="${filename%.*}"
	hmmsearch -E1 -o /home/bwadie/Hmmer_H1/human_slimfinder_HMM_outputs/all_outs/$file_name.out --domtblout /home/bwadie/Hmmer_H1/human_slimfinder_HMM_outputs/pfam_H1_domtblouts/$file_name.txt pfam_hmms/$filename human_only_H1_Prot_fasta.fasta
done

for file in /home/bwadie/Hmmer_H1/domain_all/*.hmm ; do
        filename=$(basename -- "$file")
        file_name="${filename%.*}"
        hmmsearch -E1 -o /home/bwadie/Hmmer_H1/human_slimfinder_HMM_outputs/all_outs/$file_name.out --domtblout /home/bwadie/Hmmer_H1/human_slimfinder_HMM_outputs/domain_all_H1_domtblouts/$file_name.txt domain_all/$filename human_only_H1_Prot_fasta.fasta
done
