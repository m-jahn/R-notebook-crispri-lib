
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PRE-DEPLETION ANALYSIS OF E.COLI/SYNECHOCYSTIS INITIAL SGRNA LIBRARIES
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# LOAD PACKAGES
library(lattice)
library(latticeExtra)
library(topGO)
library(tidyverse)
library(Rtools)
setwd("/home/michael/Documents/SciLifeLab/Experiments/20190122_Synechocystis_CRISPRi_library/")


# LOAD PROCESSED DATA FILE  ++++++++++++++++++++++++++++++++++++++++++++++++++++
df <- read_tsv("rawdata/190702_library_predepletion_deseq.tab")
# turn conditon into a factor for ordered plotting, input should be ungroup()ed!
df <- mutate(df, Condition = factor(Condition, c("LL", "HL", "DN")))

# add annotation to sgRNAs from existing
# map trivial names to LocusTags using manually curated list from Kiyan
# using left join function, filter for duplicates before
df <- df %>% mutate(sgRNA_short = 
    gsub("(\\_| )[0-9]+(\\_MULTI)?$", "", sgRNA)) %>% 
  
  left_join(read_tsv("rawdata/gene_protein_table.txt")[1:2] %>% 
      filter(!duplicated(gene)), 
    by = c("sgRNA_short" = "gene")) %>%
  
  # one more join, this time with the annotation db
  # containing uniprot and categorial (KEGG Brite) data
  left_join(read_csv("../../Resources/MS/databases/Synechocystis/Synechocystis_PCC6803_genome_annotation_20190614.csv")[c(1,4:6,11,20,22,27:32)], 
    by = c("locus" = "GeneID"))


# GENE SET ENRICHMENT USING TOPGO  +++++++++++++++++++++++++++++++++++++++++++++
#
# first we have a look at log2 FC disitribution over conditions
densityplot( ~ log2FoldChange | Reference * Condition, df,
  groups = padj > 0.01,
  as.table = TRUE, par.settings = custom.lattice, #xlim = c(-7.5,-0.5),
  scales = list(alternating = FALSE),
  panel = function(x, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.densityplot(x, pch = ".", lwd = 2, ...)
    #panel.key(labels = c("LL","HL","DN"), points = FALSE, lines = TRUE, lwd = 2)
  }
)


# Many are sgRNAs are changed. Let's filter by log2FC and p-value to
# construct an outlier group for GO enrichment
df <- group_by(df, Reference, Condition) %>%
  
  # we call it 'depleted' group because that's what these sgRNAs mostly are
  mutate(depleted = (abs(log2FoldChange) > 1 & padj <= 0.05) %>%
    as.numeric(.))
  


# Using GO enrichemnt, TopGO package by Alexa et al.

# the GetTopGO function takes as input either
#   1. a data frame with one row per gene as only presence in cluster is 
#      important here. Three columns are obligatory, "cluster", "GeneID", and 
#      "Gene.ontology.IDs" with GO terms separated by '; '
#   2. alternatively three vectors or lists corresponding to
#      the three columns can be provided

# reduce table to one sgRNA per row and call GetTopGO function
# execute GetTopGO function and collect list of results
TopGoResult <- mapply(SIMPLIFY = FALSE, function(cond, ref) {
    with(filter(df, Condition == cond, Reference == ref) %>% as.data.frame,
      GetTopGO(df = NULL,
        cluster = depleted,
        GeneID = sgRNA,
        Gene.ontology.IDs,
        topNodes = 50, 
        selected.cluster = 1
      )
    )
  }, cond = rep(c("LL", "HL", "DN"), 2), 
     ref = rep(c("genome", "plasmid"), each = 3)
  ) %>% setNames(c("G-LL", "G-HL", "G-DN", "P-LL", "P-HL", "P-DN")) %>%
    
  # turn into data frame with cluster as ID columns
  plyr::ldply(., .id = "Condition")
    
  # save unfiltered TopGO result
  #write_csv(TopGoResult, "TopGoResult_predepletion.csv")


# ADDITIONAL FILTERING
# optional filtering by REVIGO that excludes redundant GO terms
# To this end, the resulting terms have to be uploaded to http://revigo.irb.hr/
# and the trimmed list be downloaded. 

# merge TopGo result with REVIGO dispensibility analysis
TopGoResultFiltered <- left_join(as_tibble(TopGoResult), read_csv("REVIGO_predepletion.csv"), 
  by = c("GO.ID" = "term_ID")) %>%
  
  # trim data frame by p.value and GO term node size
  filter(dispensability <= 1.0,
    elimFisher < 0.025,
    Annotated >= 5,
    Annotated <= 200) %>%
  
  # optionally trim lengthy GO terms
  mutate(
    TermAbbr = gsub("(biosynthetic|metabolic|catabolic) proce?s?s?|\\.\\.\\.|metabolic.?.?.?", 
    replacement = "", Term) %>% paste0(., "...", Condition)
  )

# Plot results of GO enrichment for each cluster
# (p-value is probability that differential expression of interesting genes
# is identical to background)
custom_lattice = custom.lattice()
plot_TopGO_table <- xyplot(factor(TermAbbr, unique(TermAbbr)) ~ as.numeric(elimFisher), 
  TopGoResultFiltered %>% arrange(desc(Condition), desc(elimFisher)), 
  groups = factor(Condition), main = "GO enrichment, GO terms by p-value (p<0.025)",
  par.settings = custom_lattice, xlab = "p-value", ylab = "GO term",
  scales = list(alternating=FALSE, x = list(rot = 45)), 
  panel = function(x, y, groups = groups, subscripts = subscripts, ...) {
    panel.abline(h = 1:length(y), col = grey(0.9), lty = 2)
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.key(labels = as.character(TopGoResultFiltered$Condition %>% unique), 
      pch = 19, corner = c(0.1, 0.95))
    panel.superpose(x, y, groups, subscripts, ...)
  },
  panel.groups = function(x, y, groups = groups, subscripts = subscripts, ...) {
    panel.xyplot(x, y, cex = 1.5-(35*x), pch = 19,
      col = custom_lattice$superpose.polygon$col[subscripts])
  }
)

# Plot results of GO enrichment as a map of semantic similarity
plot_TopGO_map <- xyplot(plot_Y ~ plot_X, TopGoResultFiltered,
  Condition = as.numeric(TopGoResultFiltered$Condition),
  bubblesize = 3+0.08 * sapply(TopGoResultFiltered$SigGenes, function(x) {
    strsplit(x, split = ",") %>% unlist %>% length}), 
  xlim = c(-10, 10), ylim = c(-10, 10),
  main = "GO enrichment, GO terms by semantic similarity (p<0.025)",
  par.settings = custom_lattice, xlab = "semantic space X", ylab = "semantic space Y",
  panel = function(x, y, bubblesize, Condition, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.key(labels = as.character(TopGoResultFiltered$Condition %>% unique), 
      pch = 19, corner = c(0.9, 0.95))
    panel.key(labels = paste(c(10, 20, 40), "sgRNAs\n"), points = FALSE, col = grey(0.6), corner = c(0.01, 0.95))
    panel.xyplot(x = -5, y = c(9.1, 8.2, 7.3), col = grey(0.6), pch = 19, cex = 3+0.08*c(10, 20, 40))
    panel.xyplot(x, y, alpha = 0.5, cex = bubblesize, 
      col = custom_lattice$superpose.polygon$col[Condition], pch = 19)
    selected = TopGoResultFiltered$dispensability <= 0.7
    panel.text(x[selected], y[selected], cex = 0.8, col = grey(0.4),
      labels = TopGoResultFiltered$TermAbbr[selected])
  }
)

png("CRISPRi_library_TopGoResult_predepletion.png", width = 1400, height = 900, res = 100)
print(plot_TopGO_table, position = c(0, 0, 0.5, 1), more = TRUE)
print(plot_TopGO_map, position = c(0.5, 0.03, 1, 1))
dev.off()

