
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HVSlimPred

<!-- badges: start -->

<!-- badges: end -->

This R package complements the performance evaluation analysis for the
manuscript entitled “Identifying novel functional linear motifs using
the host-viral protein interaction network and the principle of
convergent evolution”.

## Installation

You can install the released version of HVSlimPred from gitlab with:

``` r
devtools::install_gitlab("petsalakilab/HVSlimPred", host = "gitlab.ebi.ac.uk")
```

## Get Protein-level evaluation metrics

For protein-level enrichment, we measured the enrichment of
true-positives in our predicted dataset using a one-tailed fisher-exact
test, where the odds ratio represents the magnitude of the enrichment.
True-positives are the number of motif-carrying proteins present in both
the predicted dataset and the ELM dataset regardless of whether the
predicted protein has the right motif or found in the right location

The output of the following command is a data frame containing all the
relevant protein-level performance metrics for each domain enrichment
filter in addition to the non-filtered qslim output.

``` r
library(HVSlimPred)
prot_eval_metrics = HV_prot_level_eval()
```

## Get Motif-level evaluation metrics

For motif-level enrichment, we simply cannot use a binary classification
as we did for protein-level evaluation because in reality, predicted
motifs are partially correct to some extent as they might contain
true-positive residues in a given sequence stretch, and therefore we
used a re-implemented version of the evaluation protocol proposed in
[Prytuliak et
al. 2017](https://academic.oup.com/nar/article/45/W1/W470/3782606)
instead of binary classification, where we computed the common
performance metrics (Recall, precision F1, etc .. ) both residue-wise
and site-wise given that the motif-carrying proteins are also found in
the ELM benchmarking dataset. So this analysis was not performed on
proteins not reported in the ELM dataset.

The output of the following command is a data frame containing all the
relevant motif-level performance metrics for each domain enrichment
filter in addition to the non-filtered qslim output.

``` r
library(HVSlimPred)
motif_eval_metrics = HV_motif_level_eval()
```

## Get Protein-domain interactions evaluation metrics

For evaluating protein-domain interactions we measured the enrichment of
true-positive interactions between a given motif-carrying protein and
its associated domains as reported in the ELM interaction dataset. As in
the motif-level evaluation, this analysis was performed only on the
motif-carrying proteins reported in the ELM interactions dataset, where
true-positives represents the number of correctly associated domains for
a given motif-carrying protein and then summed over all motif-carrying
proteins in the predicted dataset.

The output of the following command is a data frame containing all the
relevant protein-domain interactions’ performance metrics for each
domain enrichment filter.

``` r
library(HVSlimPred)
ProtDom_int_eval_metrics = HV_prot_dom_int_eval()
```
