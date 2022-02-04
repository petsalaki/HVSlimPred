library(HVSlimPred)
require(dcGOR)

# Prepare custom data for GOBP -dcGOR -------------------------------------------------

#Prepare Domain.RData
Pfam <- dcRDataLoader('Pfam')
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
pfam_doms = pfam_doms %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")
Pfam_human_bg = pfam_doms[,c(2,3)]
Pfam_human_bg = Pfam_human_bg[!duplicated(Pfam_human_bg),]
Pfam_human_bg$level = "pfam"
Pfam_human_bg = Pfam_human_bg[,c(1,3,2)]
colnames(Pfam_human_bg) = c("PFAM_ID", "level", "name")
write.table(Pfam_human_bg, "data/input/dcGOR/pfam_human_bg.txt", row.names = F, quote = F, sep = "\t")
dcBuildInfoDataFrame(input.file = "data/input/dcGOR/pfam_human_bg.txt", output.file = "data/input/dcGOR/custom_Domain.RData")

#Prepare Onto.RData
GO_obo = ontologyIndex::get_OBO("data/input/dcGOR/go-basic.obo")
GO_names = data.frame(ID = GO_obo$id, Name = GO_obo$name)
all_BP_terms = ontologyIndex::get_descendants(GO_obo, roots = "GO:0008150")
BP_relations = list()
for (i in 1:length(all_BP_terms)){
  x = GO_obo$children[[all_BP_terms[i]]]
  if (length(x) == 0){
    BP_relations[[i]] = NULL
  }
  else{
    BP_relations[[i]] = data.frame(from = all_BP_terms[i], to = x)
  }
}
BP_relations = dplyr::bind_rows(BP_relations)
BP_ig <- igraph::graph.data.frame(d = BP_relations, directed = T,
                                  vertices = all_BP_terms)

trial = igraph::distances(BP_ig, to = "GO:0008150")
distances_to_root = data.frame(ID = rownames(trial), distance = trial)
colnames(distances_to_root) = c("ID", "Distance")
GOBP_term_info = data.frame(ID = all_BP_terms, Namespace = "biological_process")
GOBP_term_info = dplyr::left_join(GOBP_term_info, GO_names)
GOBP_term_info = dplyr::left_join(GOBP_term_info, distances_to_root)
GOBP_term_info = GOBP_term_info[,c(1,3,2,4)]
write.table(GOBP_term_info, "data/input/dcGOR/GOBP_term_info.txt", row.names = F, quote = F, sep = "\t")

#Prepare association PFAM2GO file
pfam2GO = read.delim("data/input/dcGOR/pfam2go", header = F)
pfam2GO = pfam2GO %>% tidyr::separate(V1, into = c("pfam_info", "GO_info"), sep = " > ") %>%
  tidyr::separate(GO_info, into = c("rm1", "ID"), sep = " ; ") %>% dplyr::select(-rm1)
pfam2GO$pfam_info = sub(" .*", "", sub(".*Pfam:","",pfam2GO$pfam_info))
colnames(pfam2GO)[1] = "PFAM_ID"
pfam2GO = pfam2GO[which(pfam2GO$PFAM_ID %in% Pfam_human_bg$PFAM_ID),]
pfam2GO = pfam2GO[which(pfam2GO$ID %in% GO_obo$id),]
write.table(pfam2GO, "data/input/dcGOR/pfam2GO_assoc.txt", row.names = F, quote = F, sep = "\t")

#Build dcAnno object
#This will be added after line 69 in dcBuildAnnot : f_col = f_col[!is.na(f_col)]
dcBuildAnnot = function (domain_info.file, term_info.file, association.file,
          output.file = "Anno.RData")
{
  if (is.null(domain_info.file) | is.na(domain_info.file)) {
    stop("The file 'domain_info.file' must be provided!\n")
  }
  if (is.null(term_info.file) | is.na(term_info.file)) {
    stop("The file 'term_info.file' must be provided!\n")
  }
  if (is.null(association.file) | is.na(association.file)) {
    stop("The file 'association.file' must be provided!\n")
  }
  if (is.null(output.file)) {
    warnings("Since the output file is not provided, the function will use the default output file 'Anno.RData'!\n")
    output.file <- "Anno.RData"
  }
  f <- function(x) {
    if (is.numeric(x)) {
      x
    }
    else {
      iconv(x, "latin1", "ASCII", sub = "")
    }
  }
  tab <- utils::read.delim(association.file, header = F, sep = "\t",
                           nrows = 50, skip = 1)
  association <- utils::read.table(association.file, header = F,
                                   sep = "\t", skip = 1, colClasses = sapply(tab, class))
  x <- association[, 1]
  y <- association[, 2]
  z <- rep(1, nrow(association))
  domains_in_asso <- sort(unique(x))
  terms_in_asso <- sort(unique(y))
  data <- Matrix::sparseMatrix(i = match(x, domains_in_asso),
                               j = match(y, terms_in_asso), x = z)
  rownames(data) <- domains_in_asso
  colnames(data) <- terms_in_asso
  term_info <- utils::read.delim(term_info.file, header = T)
  term_info <- as.data.frame(apply(term_info, 1:2, f))
  colnames(term_info) <- c("ID", "Name", "Namespace", "Distance")
  rownames(term_info) <- term_info[, 1]
  term_info <- term_info[order(term_info[, 4]), ]
  termID <- intersect(terms_in_asso, rownames(term_info))
  if (sum(is.na(suppressWarnings(as.numeric(termID)))) >=
      1) {
    termID <- sort(termID)
  }
  else {
    termID <- sort(as.numeric(termID))
  }
  flag <- match(termID, rownames(term_info))
  term_info <- term_info[flag, ]
  term_info <- term_info[order(term_info[, 4]), ]
  domain_info <- utils::read.delim(domain_info.file, header = T)
  rownames(domain_info) <- domain_info[, 1]
  domain_info <- domain_info[order(domain_info[, 1]), ]
  domainID <- intersect(domains_in_asso, rownames(domain_info))
  if (sum(is.na(suppressWarnings(as.numeric(domainID)))) >=
      1) {
    domainID <- sort(domainID)
  }
  else {
    domainID <- sort(as.numeric(domainID))
  }
  flag <- match(domainID, rownames(domain_info))
  domain_info <- domain_info[flag, ]
  domain_info <- domain_info[order(domain_info[, 1]), ]
  f_row <- match(rownames(data), rownames(domain_info))
  f_col <- match(colnames(data), rownames(term_info))
  f_col[!is.na(f_col)]
  annoData <- data[f_row, f_col]
  x <- new("Anno", annoData = annoData, termData = as(term_info,
                                                      "InfoDataFrame"), domainData = as(domain_info, "InfoDataFrame"))
  output.var <- gsub(".RData$", "", output.file, ignore.case = T,
                     perl = T)
  output.var <- gsub(".RDat$", "", output.var, ignore.case = T,
                     perl = T)
  output.var <- gsub(".RDa$", "", output.var, ignore.case = T,
                     perl = T)
  do.call(assign, list(output.var, x))
  save(list = output.var, file = output.file)
  if (file.exists(output.file)) {
    message(sprintf("An object of S4 class 'Anno' has been built and saved into '%s'.",
                    file.path(getwd(), output.file)), appendLF = T)
  }
  invisible(x)
}
environment(dcBuildAnnot) <- asNamespace('dcGOR')
assignInNamespace("dcBuildAnnot", dcBuildAnnot, ns = "dcGOR")

dcBuildAnno(domain_info.file = "data/input/dcGOR/pfam_human_bg.txt",
            term_info.file = "data/input/dcGOR/GOBP_term_info.txt",
            association.file = "data/input/dcGOR/pfam2GO_assoc.txt",
            output.file = "data/input/dcGOR/custom_Pfam2GOBP.RData")

#Build dcOnto object
nodes_file = GOBP_term_info
nodes_file$name = nodes_file$ID
nodes_file = nodes_file[,c(5,1:4)]
colnames(nodes_file) = c("name", "term_id", "term_name", "term_namespace", "term_distance")
write.table(nodes_file, "data/input/dcGOR/GOBP_nodes.txt", row.names = F, quote = F, sep = "\t")
write.table(BP_relations, "data/input/dcGOR/GOBP_relations.txt", row.names = F, quote = F, sep = "\t")

dcBuildOnto(relations.file = "data/input/dcGOR/GOBP_relations.txt",
            nodes.file = "data/input/dcGOR/GOBP_nodes.txt",
            output.file = "data/input/dcGOR/custom_onto_GOBP.RData")

# Perform Domain Enrichment all filters -dcGOR  ---------------------------
new_hits_ints = readRDS("data/output/HV_final_hits.rds")
pepsite_requests = readRDS("data/output/All_HV_pepsite_reqs_with_AF.rds")
colnames(pepsite_requests)[which(names(pepsite_requests) == "motif_protein")] = "uniprot"
iELM = readRDS("data/output/HV_HMM_iELM.rds")
EM = readRDS("data/output/HV_enriched_doms_with_emppval.rds")
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
pfam_doms = pfam_doms %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")
human_human_intact = read.csv("data/input/human_human_intact.csv.gz")
all_data = readRDS("data/output/HV_all_filters_data.rds")

select_bg_dcGOR_filter = function(req_data, filter_sym, restrict_iELM = F, H1_dom_grp = F){
  if (stringr::str_detect(filter_sym, "A")){
    x = dplyr::bind_rows(EM$Enriched_domains)
    if ("emp_pval" %in% names(x)){
      x = x[which(x$emp_pval < 0.05),]
    }
    x = x %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")
    A_bg = unique(x$PFAM_ID)
  }
  else{
    A_bg = unique(pfam_doms$PFAM_ID)
  }
  B_bg = unique(pfam_doms$PFAM_ID)

  if (stringr::str_detect(filter_sym, "C")){
    pep_reqs = pepsite_requests
    if (H1_dom_grp){
      H1_prots = unique(sub(".*_", "", new_hits_ints$Dataset))
      pep_reqs = pep_reqs[which(pep_reqs$domain_protein %in% H1_prots),]
    }
    x = req_data[[filter_sym]]
    y = plyr::match_df(pep_reqs, x, on = c("uniprot", "Start_Pos", "End_Pos"))
    C_bg = unique(y$pfam_id)
  }
  else{
    C_bg = unique(pfam_doms$PFAM_ID)
  }

  if (stringr::str_detect(filter_sym, "D") & restrict_iELM){
    x = req_data[[filter_sym]]
    y = plyr::match_df(new_hits_ints, x, on = c("uniprot", "Start_Pos", "End_Pos"))
    y = y %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
    iELM_doms = iELM[which(iELM$H1 %in% y$H1),]
    iELM_doms = iELM_doms[which(iELM_doms$HMM_ELM_type == "pfam"),]
    D_bg = unique(sub("[.].*", "", iELM_doms$query_accession))
  }
  else{
    D_bg = unique(pfam_doms$PFAM_ID)
  }
  final_bg = Reduce(intersect, list(A_bg, B_bg, C_bg, D_bg))
  return(final_bg)
}
create_input_dcGOR_filter = function(req_data, filter_sym, restrict_iELM = F, H1_dom_grp = F){
  if (H1_dom_grp){
    H1_prots = unique(sub(".*_", "", new_hits_ints$Dataset))
    H1_doms = unique(pfam_doms$PFAM_ID[which(pfam_doms$uniprot %in% H1_prots)])
  }
  main = req_data[[filter_sym]]
  if (stringr::str_detect(filter_sym, "A")){
    x = plyr::match_df(req_data$A, main, on = c("uniprot", "Start_Pos", "End_Pos"))
    A_input = unique(x$pfam_id)
    if (H1_dom_grp){
      A_input = intersect(A_input, H1_doms)
    }
  }
  else{
    A_input = unique(pfam_doms$PFAM_ID)
  }
  if (stringr::str_detect(filter_sym, "B")){
    x = plyr::match_df(new_hits_ints, main, on = c("uniprot", "Start_Pos", "End_Pos"))
    x = x %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_") %>% dplyr::select(H1)
    if (H1_dom_grp){
      B_input = unique(pfam_doms$PFAM_ID[which(pfam_doms$uniprot %in% x$H1)])
    }
    else{
      motif_prots = unique(main$uniprot)
      all_human_ints = list()
      for (i in motif_prots){
        all_human_ints[[i]] = get_human_interactors(human_human_intact = human_human_intact, uniprot_id = i, return_fasta = F)
      }
      all_human_ints = unique(unlist(all_human_ints))
      # B_input = unique(pfam_doms$PFAM_ID)
      B_input = unique(pfam_doms$PFAM_ID[which(pfam_doms$uniprot %in% all_human_ints)])
    }
  }
  else{
    B_input = unique(pfam_doms$PFAM_ID)
  }

  if (stringr::str_detect(filter_sym, "C")){
    new_hits_pep = plyr::match_df(new_hits_ints, main, on = c("uniprot", "Start_Pos", "End_Pos"))
    if (H1_dom_grp){
      pep_H1_prots = unique(sub(".*_", "", new_hits_pep$Dataset))
      pep_hits = req_data$C[which(req_data$C$domain_protein %in% pep_H1_prots),]
      pep_hits = plyr::match_df(pep_hits, main, on = c("uniprot", "Start_Pos", "End_Pos"))
      C_input = unique(pep_hits$pfam_id)
    }
    else{
      pep_hits = plyr::match_df(req_data$C, main, on = c("uniprot", "Start_Pos", "End_Pos"))
      C_input = unique(pep_hits$pfam_id)
    }
  }
  else{
    C_input = unique(pfam_doms$PFAM_ID)
  }

  if (stringr::str_detect(filter_sym, "D")){
    if (H1_dom_grp){
      if (restrict_iELM){
        D_hits = plyr::match_df(req_data$D, main, on = c("uniprot", "Start_Pos", "End_Pos"))
        iELM_doms = iELM[which(iELM$H1 %in% D_hits$domain_protein),]
        iELM_doms = iELM_doms[which(iELM_doms$HMM_ELM_type == "pfam"),]
        D_input = unique(sub("[.].*", "", iELM_doms$query_accession))
      }
      else{
        x = plyr::match_df(new_hits_ints, main, on = c("uniprot", "Start_Pos", "End_Pos"))
        x = x %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_") %>% dplyr::select(H1)
        D_input = unique(pfam_doms$PFAM_ID[which(pfam_doms$uniprot %in% x$H1)])
      }
    }
    else{
      motif_prots = unique(main$uniprot)
      all_human_ints = list()
      for (i in motif_prots){
        all_human_ints[[i]] = get_human_interactors(human_human_intact = human_human_intact, uniprot_id = i, return_fasta = F)
      }
      all_human_ints = unique(unlist(all_human_ints))
      # B_input = unique(pfam_doms$PFAM_ID)
      D_input = unique(pfam_doms$PFAM_ID[which(pfam_doms$uniprot %in% all_human_ints)])
      # D_input = unique(pfam_doms$PFAM_ID)
    }
  }
  else{
    D_input = unique(pfam_doms$PFAM_ID)
  }
  final_input = Reduce(intersect, list(A_input, B_input, C_input, D_input))
  return(final_input)
}
prepare_data_dcGOR_enrich = function(req_data, new_hits_ints, filter_sym, restrict_iELM = F, domain_grp = c("H1", "pepsite", "EM"),min_set_size = 5){
  main = req_data[[filter_sym]]
  if (domain_grp == "H1"){
    H1_prots = unique(sub(".*_", "", new_hits_ints$Dataset))
    H1_doms = unique(pfam_doms$PFAM_ID[which(pfam_doms$uniprot %in% H1_prots)])
    bg = select_bg_dcGOR_filter(req_data = req_data, filter_sym = filter_sym, restrict_iELM = restrict_iELM, H1_dom_grp = T)
    bg = intersect(bg, H1_doms)
    if (length(bg) < min_set_size){
      return(NULL)
    }
    input = create_input_dcGOR_filter(req_data = req_data, filter_sym = filter_sym, restrict_iELM = restrict_iELM, H1_dom_grp = T)
    if (length(input) < min_set_size){
      return(NULL)
    }
    return(list("query" = input, "bg" = bg))
  }
  else{
    bg = select_bg_dcGOR_filter(req_data = req_data, filter_sym = filter_sym, restrict_iELM = F, H1_dom_grp = F)
    if (length(bg) < min_set_size){
      return(NULL)
    }
    input = create_input_dcGOR_filter(req_data = req_data, filter_sym = filter_sym, restrict_iELM = F, H1_dom_grp = F)
    if (length(input) < min_set_size){
      return(NULL)
    }
    # if (filter_sym %in% c("B", "D", "BD")){
    #   return(NULL)
    # }
    else{
      return(list("query" = input, "bg" = bg))
    }
  }

}

run_dcGOR_enrich = function(new_hits_ints,req_data, cons = F, restrict_iELM = F){
  if (cons){
    req_data = req_data[["conserved_only"]]
    All_filters_dcGOR = list()
    for (i in 1:nrow(req_data$Metadata)){
      x = prepare_data_dcGOR_enrich(req_data = req_data, new_hits_ints = new_hits_ints, filter_sym = req_data$Metadata$sym[i],domain_grp = "H1")
      y = prepare_data_dcGOR_enrich(req_data = req_data, new_hits_ints = new_hits_ints, filter_sym = req_data$Metadata$sym[i],domain_grp = "pepsite")
      All_filters_dcGOR[[req_data$Metadata$Term[i]]] = list("H1" = x, "Motif_prots" = y)
    }
  }
  else{
    req_data = req_data[["all_hits"]]
    All_filters_dcGOR = list()
    for (i in 1:nrow(req_data$Metadata)){
      x = prepare_data_dcGOR_enrich(req_data = req_data, new_hits_ints = new_hits_ints, filter_sym = req_data$Metadata$sym[i],domain_grp = "H1")
      y = prepare_data_dcGOR_enrich(req_data = req_data, new_hits_ints = new_hits_ints, filter_sym = req_data$Metadata$sym[i],domain_grp = "pepsite")
      All_filters_dcGOR[[req_data$Metadata$Term[i]]] = list("H1" = x, "Motif_prots" = y)
    }
  }

  for (i in 1:length(All_filters_dcGOR)){
    for (j in 1:2){
      if (!is.null(All_filters_dcGOR[[i]][[j]])){
        x = dcEnrichment(data = All_filters_dcGOR[[i]][[j]]$query, background = All_filters_dcGOR[[i]][[j]]$bg,
                         p.adjust.method = "BH", test = "HypergeoTest", ontology.algorithm = "lea",
                         domain.RData = "data/input/dcGOR/custom_Domain.RData",
                         ontology.RData = "data/input/dcGOR/custom_onto_GOBP.RData",
                         annotations.RData = "data/input/dcGOR/custom_Pfam2GOBP.RData")
        if (is.null(x)){
          next()
        }
        All_filters_dcGOR[[i]][[j]][["eoutput"]] = list("lea" = x)
      }
    }
  }
  return(All_filters_dcGOR)
}



dcDAGannotate = function (g, annotations, path.mode = c("all_paths", "shortest_paths",
                                                        "all_shortest_paths"), verbose = TRUE)
{
  path.mode <- match.arg(path.mode)
  if (class(g) == "Onto") {
    ig <- dcConverter(g, from = "Onto", to = "igraph", verbose = F)
  }
  else {
    ig <- g
  }
  if (class(ig) != "igraph") {
    stop("The function must apply to either 'igraph' or 'Onto' object.\n")
  }
  if (class(annotations) != "Anno") {
    stop("The function must apply to 'Anno' object.\n")
  }
  D <- annoData(annotations)
  originAnnos <- sapply(1:ncol(D), function(j) {
    names(which(D[, j] != 0))
  })
  names(originAnnos) <- colnames(D)
  if (is.list(originAnnos)) {
    originNodes <- names(originAnnos)
    ind <- match(originNodes, V(ig)$name)
    nodes_mapped <- originNodes[!is.na(ind)]
    if (length(nodes_mapped) == 0) {
      stop("The input annotations do not contain terms matched to the nodes/terms in the input graph.\n")
    }
  }
  dag <- dnet::dDAGinduce(ig, originNodes, path.mode = path.mode)
  allNodes <- V(dag)$name
  node2domain.HoH <- new.env(hash = T, parent = emptyenv())
  lapply(allNodes, function(node) {
    e <- new.env(hash = T, parent = emptyenv())
    if (node %in% originNodes) {
      sapply(originAnnos[[node]], function(domain) {
        assign(as.character(domain), "origin", envir = e)
      })
    }
    assign(node, e, envir = node2domain.HoH)
  })
  level2node <- dnet::dDAGlevel(dag, level.mode = "longest_path",
                                return.mode = "level2node")
  level2node.Hash <- list2env(level2node)
  nLevels <- length(level2node)
  for (i in nLevels:1) {
    currNodes <- get(as.character(i), envir = level2node.Hash,
                     mode = "character")
    adjNodesList <- lapply(currNodes, function(node) {
      neighs.in <- igraph::neighborhood(dag, order = 1,
                                        nodes = node, mode = "in")
      setdiff(V(dag)[unlist(neighs.in)]$name, node)
    })
    names(adjNodesList) <- currNodes
    lapply(currNodes, function(node) {
      domainsID <- ls(get(node, envir = node2domain.HoH,
                          mode = "environment"))
      lapply(adjNodesList[[node]], function(adjNode) {
        adjEnv <- get(adjNode, envir = node2domain.HoH,
                      mode = "environment")
        sapply(domainsID, function(domainID) {
          assign(domainID, "inherit", envir = adjEnv)
        })
      })
    })
    if (verbose) {
      message(sprintf("\tAt level %d, there are %d nodes, and %d incoming neighbors.",
                      i, length(currNodes), length(unique(unlist(adjNodesList)))),
              appendLF = T)
    }
  }
  node2domains <- as.list(node2domain.HoH)[allNodes]
  domain_annotations <- sapply(node2domains, function(node) {
    names(unlist(as.list(node)))
  })
  V(dag)$annotations <- domain_annotations
  counts <- sapply(domain_annotations, length)
  IC <- -1 * log10(counts/max(counts))
  V(dag)$IC <- IC
  if (class(g) == "Onto") {
    dag <- dcConverter(dag, from = "igraph", to = "Onto",
                       verbose = F)
  }
  return(dag)
}

environment(dcDAGannotate) <- asNamespace('dcGOR')
assignInNamespace("dcDAGannotate", dcDAGannotate, ns = "dcGOR")

dcEnrichment = function (data, background = NULL, domain = c(NA, "SCOP.sf",
                                                             "SCOP.fa", "Pfam", "InterPro", "Rfam"), ontology = c(NA,
                                                                                                                  "GOBP", "GOMF", "GOCC", "DO", "HPPA", "HPMI", "HPON", "MP",
                                                                                                                  "EC", "KW", "UP"), sizeRange = c(10, 1000), min.overlap = 3,
                         which_distance = NULL, test = c("HypergeoTest", "FisherTest",
                                                         "BinomialTest"), p.adjust.method = c("BH", "BY", "bonferroni",
                                                                                              "holm", "hochberg", "hommel"), ontology.algorithm = c("none",
                                                                                                                                                    "pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2,
                         verbose = T, domain.RData = NULL, ontology.RData = NULL,
                         annotations.RData = NULL, RData.location = "http://dcgor.r-forge.r-project.org/data")
{
  startT <- Sys.time()
  message(paste(c("Start at ", as.character(startT)), collapse = ""),
          appendLF = T)
  message("", appendLF = T)
  domain <- match.arg(domain)
  ontology <- match.arg(ontology)
  test <- match.arg(test)
  p.adjust.method <- match.arg(p.adjust.method)
  ontology.algorithm <- match.arg(ontology.algorithm)
  if (is.vector(data)) {
    data <- unique(data)
    data <- data[!is.null(data)]
    data <- data[!is.na(data)]
  }
  else {
    stop("The input data must be a vector.\n")
  }
  if (!is.na(domain) & !is.na(ontology)) {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("First, load the ontology '%s', the domain '%s', and their associations (%s) ...",
                      ontology, domain, as.character(now)), appendLF = T)
    }
    g <- dcRDataLoader(paste("onto.", ontology, sep = ""),
                       RData.location = RData.location)
    if (class(g) == "Onto") {
      g <- dcConverter(g, from = "Onto", to = "igraph",
                       verbose = F)
    }
    Domain <- dcRDataLoader(domain, RData.location = RData.location)
    Anno <- dcRDataLoader(domain = domain, ontology = ontology,
                          RData.location = RData.location)
  }
  else if (file.exists(domain.RData) & file.exists(ontology.RData) &
           file.exists(annotations.RData)) {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("First, load customised ontology '%s', the domain '%s', and their associations '%s' (%s)...",
                      ontology.RData, domain.RData, annotations.RData,
                      as.character(now)), appendLF = T)
    }
    g <- ""
    eval(parse(text = paste("g <- get(load('", ontology.RData,
                            "'))", sep = "")))
    if (class(g) == "Onto") {
      g <- dcConverter(g, from = "Onto", to = "igraph",
                       verbose = F)
    }
    ontology <- ontology.RData
    Domain <- ""
    eval(parse(text = paste("Domain <- get(load('", domain.RData,
                            "'))", sep = "")))
    domain <- domain.RData
    Anno <- ""
    eval(parse(text = paste("Anno <- get(load('", annotations.RData,
                            "'))", sep = "")))
  }
  else {
    stop("There is no input for domains and/or ontology and/or annotation.\n")
  }
  if (1) {
    if (is.vector(background)) {
      background <- unique(background)
      background <- background[!is.null(background)]
      background <- background[!is.na(background)]
    }
    if (length(background) > 0) {
      ind <- match(background, rowNames(Domain))
      background <- background[!is.na(ind)]
      if (length(background) > 0) {
        background <- union(background, data)
        ind <- match(domainNames(Anno), background)
        Anno <- Anno[!is.na(ind), ]
      }
    }
  }
  dag <- dcDAGannotate(g, annotations = Anno, path.mode = "all_paths",
                       verbose = F)
  ind <- match(data, rowNames(Domain))
  domains.group <- data[!is.na(ind)]
  distance <- V(dag)$term_distance
  if (!is.null(which_distance) & sum(is.na(distance)) == 0) {
    set_filtered <- sapply(which_distance, function(x) {
      V(dag)$term_id[(distance == as.integer(x))]
    })
    set_filtered <- unlist(set_filtered)
  }
  else {
    set_filtered <- V(dag)$term_id
  }
  dag <- dnet::dDAGinduce(dag, nodes_query = set_filtered,
                          path.mode = "all_paths")
  gs.length <- sapply(V(dag)$annotations, length)
  ind.length <- which(gs.length >= sizeRange[1] & gs.length <=
                        sizeRange[2])
  dag <- dnet::dDAGinduce(dag, nodes_query = V(dag)$term_id[ind.length],
                          path.mode = "all_paths")
  if (length(V(dag)) == 0) {
    stop("There is no term being used.\n")
  }
  doFisherTest <- function(domains.group, domains.term, domains.universe) {
    domains.hit <- intersect(domains.group, domains.term)
    X <- length(domains.hit)
    K <- length(domains.group)
    M <- length(domains.term)
    N <- length(domains.universe)
    cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2,
                   dimnames = list(c("anno", "notAnno"), c("group",
                                                           "notGroup")))
    p.value <- ifelse(all(cTab == 0), 1, stats::fisher.test(cTab,
                                                            alternative = "greater")$p.value)
    return(p.value)
  }
  doHypergeoTest <- function(domains.group, domains.term,
                             domains.universe) {
    domains.hit <- intersect(domains.group, domains.term)
    X <- length(domains.hit)
    K <- length(domains.group)
    M <- length(domains.term)
    N <- length(domains.universe)
    x <- X
    m <- M
    n <- N - M
    k <- K
    p.value <- ifelse(m == 0 || k == 0, 1, stats::phyper(x,
                                                         m, n, k, lower.tail = F, log.p = F))
    return(p.value)
  }
  doBinomialTest <- function(domains.group, domains.term,
                             domains.universe) {
    domains.hit <- intersect(domains.group, domains.term)
    X <- length(domains.hit)
    K <- length(domains.group)
    M <- length(domains.term)
    N <- length(domains.universe)
    p.value <- ifelse(K == 0 || M == 0 || N == 0, 1, stats::pbinom(X,
                                                                   K, M/N, lower.tail = F, log.p = F))
    return(p.value)
  }
  zscoreHyper <- function(domains.group, domains.term, domains.universe) {
    domains.hit <- intersect(domains.group, domains.term)
    X <- length(domains.hit)
    K <- length(domains.group)
    M <- length(domains.term)
    N <- length(domains.universe)
    if (1) {
      x.exp <- K * M/N
      var.exp <- K * M/N * (N - M)/N * (N - K)/(N - 1)
      if (var.exp == 0) {
        z <- NA
      }
      else {
        suppressWarnings(z <- (X - x.exp)/sqrt(var.exp))
      }
    }
    else {
      x <- X
      m <- M
      n <- N - M
      k <- K
      suppressWarnings(d <- stats::dhyper(x, m, n, k,
                                          log = TRUE) - log(2))
      suppressWarnings(pupper <- stats::phyper(x, m, n,
                                               k, lower.tail = FALSE, log.p = TRUE))
      suppressWarnings(plower <- stats::phyper(x - 1,
                                               m, n, k, lower.tail = TRUE, log.p = TRUE))
      d[is.na(d)] <- -Inf
      pupper[is.na(pupper)] <- -Inf
      plower[is.na(plower)] <- -Inf
      a <- pupper
      b <- d - pupper
      a[b > 0] <- d[b > 0]
      b <- -abs(b)
      pmidupper <- a + log1p(exp(b))
      pmidupper[is.infinite(a)] <- a[is.infinite(a)]
      a <- plower
      b <- d - plower
      a[b > 0] <- d[b > 0]
      b <- -abs(b)
      pmidlower <- a + log1p(exp(b))
      pmidlower[is.infinite(a)] <- a[is.infinite(a)]
      up <- pmidupper < pmidlower
      if (any(up))
        z <- stats::qnorm(pmidupper, lower.tail = FALSE,
                          log.p = TRUE)
      if (any(!up))
        z <- stats::qnorm(pmidlower, lower.tail = TRUE,
                          log.p = TRUE)
    }
    return(z)
  }
  terms <- V(dag)$term_id
  gs <- V(dag)$annotations
  names(gs) <- terms
  domains.universe <- unique(unlist(V(dag)$annotations))
  domains.group <- intersect(domains.universe, domains.group)
  if (length(domains.group) == 0) {
    warnings("There is no domain being used.\n")
    return(F)
  }
  if (ontology.algorithm == "none") {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("Second, perform enrichment analysis using %s (%s) ...",
                      test, as.character(now)), appendLF = T)
      if (is.null(which_distance)) {
        message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations",
                        length(terms), paste(sizeRange, collapse = ",")),
                appendLF = T)
      }
      else {
        message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations and [%s] distance",
                        length(terms), paste(sizeRange, collapse = ","),
                        paste(which_distance, collapse = ",")), appendLF = T)
      }
    }
    pvals <- sapply(terms, function(term) {
      domains.term <- unique(unlist(gs[term]))
      p.value <- switch(test, FisherTest = doFisherTest(domains.group,
                                                        domains.term, domains.universe), HypergeoTest = doHypergeoTest(domains.group,
                                                                                                                       domains.term, domains.universe), BinomialTest = doBinomialTest(domains.group,
                                                                                                                                                                                      domains.term, domains.universe))
    })
    zscores <- sapply(terms, function(term) {
      domains.term <- unique(unlist(gs[term]))
      zscoreHyper(domains.group, domains.term, domains.universe)
    })
  }
  else if (ontology.algorithm == "pc" || ontology.algorithm ==
           "elim" || ontology.algorithm == "lea") {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("Third, perform enrichment analysis using %s based on %s algorithm to respect ontology structure (%s) ...",
                      test, ontology.algorithm, as.character(now)),
              appendLF = T)
    }
    subg <- dag
    if (verbose) {
      message(sprintf("\tThere are %d terms being used",
                      length(V(subg))), appendLF = T)
    }
    level2node <- dnet::dDAGlevel(subg, level.mode = "longest_path",
                                  return.mode = "level2node")
    level2node.Hash <- list2env(level2node)
    nLevels <- length(level2node)
    node2pval.Hash <- new.env(hash = T, parent = emptyenv())
    node2zscore.Hash <- new.env(hash = T, parent = emptyenv())
    if (ontology.algorithm == "pc") {
      for (i in nLevels:2) {
        currNodes <- get(as.character(i), envir = level2node.Hash,
                         mode = "character")
        for (currNode in currNodes) {
          domains.term <- unique(unlist(gs[currNode]))
          pvalue_whole <- switch(test, FisherTest = doFisherTest(domains.group,
                                                                 domains.term, domains.universe), HypergeoTest = doHypergeoTest(domains.group,
                                                                                                                                domains.term, domains.universe), BinomialTest = doBinomialTest(domains.group,
                                                                                                                                                                                               domains.term, domains.universe))
          zscore_whole <- zscoreHyper(domains.group,
                                      domains.term, domains.universe)
          neighs.in <- igraph::neighborhood(subg, order = 1,
                                            nodes = currNode, mode = "in")
          adjNodes <- setdiff(V(subg)[unlist(neighs.in)]$name,
                              currNode)
          domains.parent <- unique(unlist(gs[adjNodes]))
          domains.group.parent <- intersect(domains.group,
                                            domains.parent)
          domains.term.parent <- intersect(domains.term,
                                           domains.parent)
          pvalue_relative <- switch(test, FisherTest = doFisherTest(domains.group.parent,
                                                                    domains.term.parent, domains.parent), HypergeoTest = doHypergeoTest(domains.group.parent,
                                                                                                                                        domains.term.parent, domains.parent), BinomialTest = doBinomialTest(domains.group.parent,
                                                                                                                                                                                                            domains.term.parent, domains.parent))
          zscore_relative <- zscoreHyper(domains.group.parent,
                                         domains.term.parent, domains.parent)
          pvalue <- max(pvalue_whole, pvalue_relative)
          assign(currNode, pvalue, envir = node2pval.Hash)
          zscore <- ifelse(pvalue_whole > pvalue_relative,
                           zscore_whole, zscore_relative)
          assign(currNode, zscore, envir = node2zscore.Hash)
        }
        if (verbose) {
          message(sprintf("\tAt level %d, there are %d nodes/terms",
                          i, length(currNodes), appendLF = T))
        }
      }
      root <- dnet::dDAGroot(subg)
      assign(root, 1, envir = node2pval.Hash)
      assign(root, 0, envir = node2zscore.Hash)
    }
    else if (ontology.algorithm == "elim") {
      sigNode2pval.Hash <- new.env(hash = T, parent = emptyenv())
      ancNode2domain.Hash <- new.env(hash = T, parent = emptyenv())
      if (is.null(elim.pvalue) || is.na(elim.pvalue) ||
          elim.pvalue > 1 || elim.pvalue < 0) {
        elim.pvalue <- 0.01
      }
      pval.cutoff <- elim.pvalue
      for (i in nLevels:1) {
        currNodes <- get(as.character(i), envir = level2node.Hash,
                         mode = "character")
        currAnno <- gs[currNodes]
        for (currNode in currNodes) {
          domains.term <- unique(unlist(gs[currNode]))
          if (exists(currNode, envir = ancNode2domain.Hash,
                     mode = "numeric")) {
            domains.elim <- get(currNode, envir = ancNode2domain.Hash,
                                mode = "numeric")
            domains.term <- setdiff(domains.term, domains.elim)
          }
          pvalue <- switch(test, FisherTest = doFisherTest(domains.group,
                                                           domains.term, domains.universe), HypergeoTest = doHypergeoTest(domains.group,
                                                                                                                          domains.term, domains.universe), BinomialTest = doBinomialTest(domains.group,
                                                                                                                                                                                         domains.term, domains.universe))
          zscore <- zscoreHyper(domains.group, domains.term,
                                domains.universe)
          assign(currNode, pvalue, envir = node2pval.Hash)
          assign(currNode, zscore, envir = node2zscore.Hash)
          if (pvalue < pval.cutoff) {
            assign(currNode, pvalue, envir = sigNode2pval.Hash)
            elimGenesID <- currAnno[[currNode]]
            dag.ancestors <- dnet::dDAGinduce(subg,
                                              currNode, path.mode = "all_paths")
            ancestors <- setdiff(V(dag.ancestors)$name,
                                 currNode)
            oldAncestors2GenesID <- sapply(ancestors,
                                           function(ancestor) {
                                             if (exists(ancestor, envir = ancNode2domain.Hash,
                                                        mode = "numeric")) {
                                               get(ancestor, envir = ancNode2domain.Hash,
                                                   mode = "numeric")
                                             }
                                           })
            newAncestors2GenesID <- lapply(oldAncestors2GenesID,
                                           function(oldGenes) {
                                             union(oldGenes, elimGenesID)
                                           })
            if (length(newAncestors2GenesID) > 0) {
              sapply(names(newAncestors2GenesID), function(ancestor) {
                assign(ancestor, newAncestors2GenesID[[ancestor]],
                       envir = ancNode2domain.Hash)
              })
            }
          }
        }
        if (verbose) {
          num.signodes <- length(ls(sigNode2pval.Hash))
          num.ancnodes <- length(ls(ancNode2domain.Hash))
          num.elimdomains <- length(unique(unlist(as.list(ancNode2domain.Hash))))
          message(sprintf("\tAt level %d, there are %d nodes/terms: up to %d significant nodes, %d ancestral nodes changed (%d domains eliminated)",
                          i, length(currNodes), num.signodes, num.ancnodes,
                          num.elimdomains), appendLF = T)
        }
      }
    }
    else if (ontology.algorithm == "lea") {
      node2pvalo.Hash <- new.env(hash = T, parent = emptyenv())
      if (is.null(lea.depth) || is.na(lea.depth) || lea.depth <
          0) {
        lea.depth <- 2
      }
      depth.cutoff <- as.integer(lea.depth)
      for (i in nLevels:1) {
        currNodes <- get(as.character(i), envir = level2node.Hash,
                         mode = "character")
        currAnno <- gs[currNodes]
        num.recalculate <- 0
        for (currNode in currNodes) {
          domains.term <- unique(unlist(gs[currNode]))
          pvalue.old <- switch(test, FisherTest = doFisherTest(domains.group,
                                                               domains.term, domains.universe), HypergeoTest = doHypergeoTest(domains.group,
                                                                                                                              domains.term, domains.universe), BinomialTest = doBinomialTest(domains.group,
                                                                                                                                                                                             domains.term, domains.universe))
          zscore.old <- zscoreHyper(domains.group, domains.term,
                                    domains.universe)
          assign(currNode, pvalue.old, envir = node2pvalo.Hash)
          neighs.out <- igraph::neighborhood(subg, order = depth.cutoff,
                                             nodes = currNode, mode = "out")
          adjNodes <- setdiff(V(subg)[unlist(neighs.out)]$name,
                              currNode)
          if (length(adjNodes) != 0) {
            if (1) {
              pvalue.children <- sapply(adjNodes, function(child) {
                if (exists(child, envir = node2pvalo.Hash,
                           mode = "numeric")) {
                  get(child, envir = node2pvalo.Hash,
                      mode = "numeric")
                }
              })
            }
            else {
              pvalue.children <- sapply(adjNodes, function(child) {
                if (exists(child, envir = node2pval.Hash,
                           mode = "numeric")) {
                  get(child, envir = node2pval.Hash,
                      mode = "numeric")
                }
              })
            }
            chNodes <- names(pvalue.children[pvalue.children <
                                               pvalue.old])
            if (length(chNodes) > 0) {
              num.recalculate <- num.recalculate + 1
              domains.elim <- unique(unlist(gs[chNodes]))
              domains.term.new <- setdiff(domains.term,
                                          domains.elim)
              pvalue.new <- switch(test, FisherTest = doFisherTest(domains.group,
                                                                   domains.term.new, domains.universe),
                                   HypergeoTest = doHypergeoTest(domains.group,
                                                                 domains.term.new, domains.universe),
                                   BinomialTest = doBinomialTest(domains.group,
                                                                 domains.term.new, domains.universe))
              zscore.new <- zscoreHyper(domains.group,
                                        domains.term.new, domains.universe)
              pvalue <- max(pvalue.new, pvalue.old)
              zscore <- ifelse(pvalue.new > pvalue.old,
                               zscore.new, zscore.old)
            }
            else {
              pvalue <- pvalue.old
              zscore <- zscore.old
            }
          }
          else {
            pvalue <- pvalue.old
            zscore <- zscore.old
          }
          assign(currNode, pvalue, envir = node2pval.Hash)
          assign(currNode, zscore, envir = node2zscore.Hash)
        }
        if (verbose) {
          message(sprintf("\tAt level %d, there are %d nodes/terms and %d being recalculated",
                          i, length(currNodes), num.recalculate),
                  appendLF = T)
        }
      }
    }
    pvals <- unlist(as.list(node2pval.Hash))
    zscores <- unlist(as.list(node2zscore.Hash))
  }
  if (verbose) {
    now <- Sys.time()
    message(sprintf("Last, adjust the p-values using the %s method (%s) ...",
                    p.adjust.method, as.character(now)), appendLF = T)
  }
  overlaps <- sapply(names(gs), function(term) {
    domains.term <- unique(unlist(gs[term]))
    intersect(domains.group, domains.term)
  })
  flag_filter <- sapply(overlaps, function(x) ifelse(length(x) >=
                                                       min.overlap, T, F))
  if (sum(flag_filter) == 0) {
    warnings("It seems there are no terms meeting the specified 'sizeRange' and 'min.overlap'.\n")
    return(F)
  }
  gs <- gs[flag_filter]
  overlaps <- overlaps[flag_filter]
  common <- intersect(names(gs), names(zscores))
  ind_gs <- match(common, names(gs))
  ind_zscores <- match(common, names(zscores))
  gs <- gs[ind_gs[!is.na(ind_gs)]]
  overlaps <- overlaps[ind_gs[!is.na(ind_gs)]]
  zscores <- zscores[ind_zscores[!is.na(ind_zscores)]]
  pvals <- pvals[ind_zscores[!is.na(ind_zscores)]]
  flag <- !is.na(zscores)
  gs <- gs[flag]
  overlaps <- overlaps[flag]
  zscores <- zscores[flag]
  pvals <- pvals[flag]

  if (length(zscores) == 0 | length(pvals) == 0){
    message("No GO terms are enriched")
    return(NULL)
  }

  zscores <- signif(zscores, digits = 3)
  pvals <- sapply(pvals, function(x) min(x, 1))
  adjpvals <- stats::p.adjust(pvals, method = p.adjust.method)
  pvals <- signif(pvals, digits = 2)
  adjpvals <- sapply(adjpvals, function(x) min(x, 1))
  adjpvals <- signif(adjpvals, digits = 2)
  endT <- Sys.time()
  message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
          appendLF = T)
  runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"),
                                 strptime(startT, "%Y-%m-%d %H:%M:%S"), units = "secs"))
  message(paste(c("Runtime in total is: ", runTime, " secs\n"),
                collapse = ""), appendLF = T)
  set_info <- get.data.frame(dag, what = "vertices")[names(gs),
                                                     c(2:5, 7)]
  annotations <- V(dag)$annotations
  names(annotations) <- V(dag)$term_id
  annotations <- annotations[names(gs)]
  if (1) {
    overlaps <- lapply(overlaps, function(x) {
      ind <- match(x, rowNames(Domain))
      names(x) <- as.character(Domain@data$description[ind])
      x
    })
  }
  eoutput <- new("Eoutput", domain = domain, ontology = ontology,
                 term_info = set_info, anno = annotations, data = domains.group,
                 background = domains.universe, overlap = overlaps, zscore = zscores,
                 pvalue = pvals, adjp = adjpvals)
  invisible(eoutput)
}

environment(dcEnrichment) <- asNamespace('dcGOR')
assignInNamespace("dcEnrichment", dcEnrichment, ns = "dcGOR")

dDAGlevel = function (g, level.mode = c("longest_path", "shortest_path"),
                      return.mode = c("node2level", "level2node"))
{
  level.mode <- match.arg(level.mode)
  return.mode <- match.arg(return.mode)
  if (class(g) == "graphNEL") {
    ig <- igraph.from.graphNEL(g)
  }
  else {
    ig <- g
  }
  if (class(ig) != "igraph") {
    stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
  }
  root <- dDAGroot(ig)
  if (length(root) > 1) {
    warning("The input DAG has multiple roots; recalculate the root after being reversed.\n")
    ig <- dDAGreverse(ig)
    root <- dDAGroot(ig)
  }
  if (is.null(root)) {
    stop("The function must have the root; check the eligibility of the input DAG.\n")
  }
  else if (length(root) > 1) {
    stop("Even after automatic reversing, the input DAG still has multiple roots; check the eligibility of the input DAG.\n")
  }
  if (level.mode == "longest_path") {
    edgelist <- get.data.frame(ig, what = "edges")
    nodeLookDown <- new.env(hash = T, parent = emptyenv())
    parents <- root
    level <- 1
    while (length(parents) > 0) {
      sapply(parents, function(parent) {
        assign(parent, level, envir = nodeLookDown)
      })
      children <- sapply(parents, function(parent) {
        edgelist[edgelist[, 1] == parent, 2]
      })
      level <- level + 1
      parents <- unique(unlist(children, use.names = F))
    }
    node2level <- unlist(as.list(nodeLookDown))
    node2level <- node2level[V(ig)$name]
  }
  else if (level.mode == "shortest_path") {
    vpaths <- get.shortest.paths(ig, from = root, to = V(ig),
                                 output = "vpath")
    if (length(vpaths) != length(V(ig)$name)) {
      vpaths <- vpaths$vpath
    }
    node2level <- sapply(1:length(vpaths), function(i) length(vpaths[[i]]))
    names(node2level) <- V(ig)$name
  }
  if (return.mode == "node2level") {
    return(node2level)
  }
  else if (return.mode == "level2node") {
    lvs <- sort(unique(node2level))
    level2node <- lapply(lvs, function(x) names(node2level[node2level ==
                                                             x]))
    names(level2node) <- lvs
    return(level2node)
  }
}

environment(dDAGlevel) <- asNamespace('dnet')
assignInNamespace("dDAGlevel", dDAGlevel, ns = "dnet")


All_filters_dcGOR = run_dcGOR_enrich(new_hits_ints = new_hits_ints, req_data = all_data, cons = F)
All_filters_dcGOR_cons = run_dcGOR_enrich(new_hits_ints = new_hits_ints, req_data = all_data, cons = T)


saveRDS(All_filters_dcGOR, "data/output/plot_objects/HV_all_filters_dcGOR.rds")
saveRDS(All_filters_dcGOR_cons, "data/output/plot_objects/All_filters_dcGOR_conserved.rds")


