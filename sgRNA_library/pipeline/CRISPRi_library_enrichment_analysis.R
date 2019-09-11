
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CLUSTERING AND ENRICHMENT / DEPLETION ANALYSIS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# LOAD PACKAGES
library(lattice)
library(latticeExtra)
library(cluster)
library(dendextend)
library(topGO)
library(tidyverse)
library(Rtools)
library(grid)
setwd("/home/michael/Documents/SciLifeLab/Experiments/20190122_Synechocystis_CRISPRi_library/")


# LOAD PROCESSED DATA FILE  ++++++++++++++++++++++++++++++++++++++++++++++++++++
load(file = "processed_data/CRISPRi_library_df_annotated.Rdata")
# turn conditon into a factor for ordered plotting, input should be ungroup()ed!
df <- ungroup(df) %>% mutate(condition = factor(condition, c("LL", "HL", "DN", "LAC", "NACL")))


# BASIC STATISTICS  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Basic key numbers: number of unique sgRNAs and genes
paste(df$sgRNA %>% unique %>% length, "unique sgRNAs in all conditions")
paste(df$locus %>% unique %>% length, "unique genes covered")
# missing sgRNA quantification per condition, high for LAC
df %>% group_by(condition, induction, sgRNA) %>%
  summarise(NAs = sum(is.na(read.fraction.mean))) %>% summarise(sum(NAs > 0))


# distribution of raw intensities (read fractions)
# shows that the distribution is skewed towards the lower end 
# meaning that sgRNAs are being depleted
# there's a difference for DN, HL and LL, LL shows less changes over time and
# bimodal distribution (sgRNAs are not normal-distributed from start)

png("CRISPRi_library_read_density.png", width = 800, height = 1000, res = 110)
densityplot( ~ log10(read.fraction.mean) | 
    plyr::mapvalues(induction, c("i", "u"), c("induced", "uninduced")) * condition, 
  groups = timepoint, df,
  as.table = TRUE, par.settings = custom.lattice, xlim = c(-7.5,-0.5),
  scales = list(alternating = FALSE),
  panel = function(x, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.densityplot(x, pch = ".", lwd = 2, ...)
    panel.key(labels = c("0d","1d","2d","4d","8d","16d","32d"), points = FALSE, lines = TRUE, lwd = 2)
  }
)
dev.off()


# correlation between sgRNA 1 and 2 per gene, and condition
png("CRSIPRi_library_2sgRNA_correlation.png", width=1100, height=1000, res=110)
xyplot(`sgRNA-2` ~ `sgRNA-1` | factor(timepoint) * condition,
  # filter out genes that don't have two sgRNAs
  df[c(2,3,4,9,14,15,16)] %>% group_by(sgRNA_short) %>%
    filter(timepoint > 0) %>%
    spread(key = sgRNA_index, value = log2FoldChange) %>%
    rename_at(vars(6:7), function(x) paste0("sgRNA-", x)),
  groups = induction, as.table = TRUE, #layout = c(6, 4), 
  par.settings = custom.lattice, cex = 0.2, pch=19,
  scales = list(alternating = FALSE),
  xlim = c(-9, 4), ylim = c(-9, 4),
  xlab = "log2 FC sgRNA1", ylab = "log2 FC sgRNA2",
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.xyplot(x, y, ... )
    if (panel.number() %in% c(1:5, 7:11, 17:18, 23:24, 29:30)) {
      panel.quadrants(x, y, h = -1, v = -1, lwd = 2, ...)
    }
    panel.key(which.panel = 6, c("induced", "uninduced"), pch = 19, corner = c(0.1, 0.9))
  }
)
dev.off()


# ANALYSE ENRICHED SGRNAS IN LACTATE SAMPLES  ++++++++++++++++++++++++++++++++++
#
# volcanoplot with potential boundaries for enriched sgRNAs
plot_lactate_volcano <- df %>% filter(condition %in%  c("LAC", "NACL") & timepoint > 8) %>% 
  filter(!is.na(log2FoldChange)) %>%
  
  xyplot(-log10(padj) ~ log2FoldChange | 
      paste0(condition, ", ", timepoint, " d") %>% factor(., unique(.)[c(2,3,1,4)]), ., 
    as.table = TRUE, groups = induction, layout = c(4, 1),
    # groups = read.fraction.mean < 1e-05,
    par.settings = custom.lattice, 
    xlab = "log2 FC", ylab = "-log10 p-value",
    col = c("#D33F6A", "#3D8900"),
    cex = 0.2, pch=19, between = list(x = 0.6, y = 0.6),
    scales = list(alternating = FALSE),
    ylim = c(- 20, 200), xlim = c(-10, 10),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.xyplot(x, y, ...)
      panel.quadrants(x, y, h = 20, v = 2, labels = "events")
      panel.key(c("induced", "uninduced"), col = c("#D33F6A", "#3D8900"), pch = 19, corner = c(0.9, 0.7))
    }
  )

# filter sgRNAs according to pvalue and FC boundaries
candidates <- filter(df, condition %in% "LAC", induction == "i",
  log2FoldChange > 2, log10(padj) <= -20) %>%
  pull(sgRNA) %>% unique
# check which ones fall into same criteria for NACL control ('false positives')
falsepos <- filter(df, condition %in% "NACL", 
  log2FoldChange > 2 & log10(padj) <= -20) %>% 
  pull(sgRNA) %>% unique
# no overlap between Lactate and NaCl enriched sgRNAs
candidates %in% falsepos


# # alternatively filter for sgRNAs that have low fitness score in LAC but not
# # in NACL, meaning genes that could be specific transporters/tolerance genes
# # for LAC
# candidates <- filter(df, condition %in% c("NACL", "LAC"), induction == "i", timepoint == 0) %>%
#   
#   # optional filter for sgRNAs that have both sgRNAs detected in NACL and LAC
#   group_by(sgRNA_short) %>% filter(length(sgRNA_index) == 4) %>%
#   
#   # filter by highest fitness diff between sgRNA of two conditions (NACL-LAC)
#   # higher for stronger LAC depletion
#   group_by(sgRNA) %>% summarise(diff_fitness_score = diff(fitness_score)) %>% 
#   
#   # sort by differential fitness score and pull out candidates
#   arrange(desc(diff_fitness_score)) %>% slice(1:200) %>% pull(sgRNA)


# generate subset df with candidate sgRNAs only
df_lac <- filter(df, sgRNA %in% candidates & condition %in% c("NACL", "LAC"))

# plot induced versus uninduced lactate
plot_lactate_induction <- xyplot(log2FoldChange ~ generations | Process.abbr, 
  filter(df_lac, induction == "i", condition == "LAC"), 
  as.table = TRUE, groups = sgRNA,
  main = paste0("lactate, induced versus uninduced control, n = ", length(candidates)),
  par.settings = custom.lattice, type = "l",
  ylim = c(-1, 15), xlim = c(0, 30.5), #c(-7, 1)
  scales = list(alternating = FALSE),
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.xyarea(x, y, col.line = paste0(custom.lattice()$superpose.polygon$col[1], "50"), 
      border = NA, origin = 0, ...)
    panel.text(0, 12, cex = 0.7, pos = 4, col = grey(0.4), # y = -6
      paste0("n = ", length(x)/3, " sgRNAs"))
    panel.key(c("induced", "uninduced"), points = FALSE, lines = TRUE, cex = 0.7, corner = c(0.1, 0.6))
  }
) + 
as.layer(
  xyplot(log2FoldChange ~ generations | Process.abbr, 
    filter(df_lac, induction == "u", condition == "LAC"),
    type = "l", groups = sgRNA, 
    panel = function(x, y, ...) {
      panel.xyarea(x, y, col.line = paste0(custom.lattice()$superpose.polygon$col[2], "50"), 
        border = NA, origin = 0, ...)
    }
  )
)

# plot induced lactate versus induced NaCl control
plot_lactate_vs_nacl <- xyplot(log2FoldChange ~ generations | Process.abbr, 
  filter(df_lac, induction == "i", condition == "LAC"), 
  as.table = TRUE, groups = sgRNA,
  main = paste0("lactate versus NaCl control (both induced), n = ", length(candidates)),
  par.settings = custom.lattice, type = "l",
  ylim = c(-1, 15), xlim = c(0, 30.5), # c(-7, 1)
  scales = list(alternating = FALSE),
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.xyarea(x, y, col.line = paste0(custom.lattice()$superpose.polygon$col[1], "50"), 
      border = NA, origin = 0, ...)
    panel.text(0, 12, cex = 0.7, pos = 4, col = grey(0.4), #y = -6
      paste0("n = ", length(x)/3, " sgRNAs"))
    panel.key(c("lactate", "NaCl control"), points = FALSE, lines = TRUE, cex = 0.7, corner = c(0.1, 0.6))
  }
) + 
as.layer(
  xyplot(log2FoldChange ~ generations | Process.abbr, 
    filter(df_lac, induction == "i", condition == "NACL"),
    type = "l", groups = sgRNA, 
    panel = function(x, y, ...) {
      panel.xyarea(x, y, col.line = paste0(custom.lattice()$superpose.polygon$col[2], "50"), 
        border = NA, origin = 0, ...)
    }
  )
)

png("CRISPRi_library_LAC_enrichment_process.png", width = 900, height = 1400, res = 120)
print(plot_lactate_induction, position = c(0,0.5,1,1), more = TRUE)
print(plot_lactate_vs_nacl, position = c(0,0,1,0.5))
dev.off()


# subset data according to 1 or 2 sgRNA per gene enriched
selected_fitness <- df_lac %>% group_by(sgRNA_short) %>% 
  filter(condition == "LAC", induction == "i") %>% 
  summarise(
    average_fitness = mean(unique(fitness_score), na.rm =TRUE),
    Process = unique(Process.abbr),
    sgRNAs = length(unique(sgRNA_index))
  ) %>% arrange(Process)

plot_1sgRNAs <- xyplot(log2FoldChange ~ generations | 
    factor(sgRNA_short, arrange(df_lac, Process.abbr) %>% 
      pull(sgRNA_short) %>% unique),
  df_lac %>% group_by(condition, sgRNA_short) %>%
  filter(condition == "LAC", induction == "i", length(sgRNA_index) == 6), 
  auto.key = list(columns = 5, points = FALSE, lines = TRUE),
  main = "genes with 1 enriched sgRNA",
  as.table = TRUE, groups = Process.abbr,
  par.settings = custom.lattice, type = "p",
  scales = list(alternating = FALSE, cex = 0.8),
  ylim = c(-1, 15), xlim = c(0, 30.5), # c(-7, 1)
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.xyarea(x[1:3], y[1:3], border = NA, alpha = 0.7, origin = 0, ...)
    panel.text(0, 10, cex = 0.6, pos = 4, col = grey(0.4), # -4
      paste0("F = ", round(filter(selected_fitness, sgRNAs == 1)[panel.number(), 2], 2)))
  }
)


plot_2sgRNAs <- xyplot(log2FoldChange ~ generations | 
  factor(sgRNA_short, arrange(df_lac, Process.abbr) %>% 
    pull(sgRNA_short) %>% unique),
  df_lac %>% group_by(condition, sgRNA_short) %>%
    filter(condition == "LAC", induction == "i", length(sgRNA_index) == 12), 
  auto.key = list(columns = 3, points = FALSE, lines = TRUE), layout = c(4, 3),
  #main = "genes with 2 enriched sgRNAs",
  as.table = TRUE, groups = Process.abbr, ylab = "log2 FC",
  par.settings = custom.lattice, 
  between = list(x = 0.6, y = 0.6),
  scales = list(alternating = FALSE, cex = 0.8),
  ylim = c(-10, 10), xlim = c(0, 30.5), # c(-7, 1)
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.xyarea(x[1:3], y[1:3], border = NA, alpha = 0.5, origin = 0, ...)
    panel.xyarea(x[4:6], y[4:6], border = NA, alpha = 0.5, origin = 0, ...)
    panel.text(0, 8, cex = 0.7, pos = 4, col = grey(0.4),  # -4
      paste0("F = ", round(filter(selected_fitness, sgRNAs == 2)[panel.number(), 2], 2)))
  }
) +
as.layer(
  xyplot(log2FoldChange ~ generations | 
    factor(sgRNA_short, arrange(df_lac, Process.abbr)
      %>% pull(sgRNA_short) %>% unique),
    df_lac %>% group_by(condition, sgRNA_short) %>%
      filter(condition == "NACL", induction == "i", length(sgRNA_index) == 6),
    type = "l", groups = sgRNA, 
    panel = function(x, y, ...) {
      panel.xyarea(x, y, col.line = grey(0.5, 0.5), border = NA, origin = 0, ...)
    }
  )
)


png("CRISPRi_library_LAC_enrichment_process2.png", width = 1600, height = 700, res = 110)
print(plot_2sgRNAs, position = c(0.6, 0, 1, 1), more = TRUE)
print(plot_1sgRNAs, position = c(0, 0, 0.6, 1))
dev.off()


# FIGURE 4
#png("CRISPRi_library_volcanoplot.png", width = 1000, height = 1400, res = 110)
svg("CRISPRi_library_volcanoplot.svg", width = 6.5, height = 7.7)
print(plot_lactate_volcano, position = c(0, 0.64, 1, 1), more = TRUE)
print(plot_2sgRNAs, position = c(0.01, 0, 1, 0.66))
grid.text(label = c("A", "B"), x = c(0.03, 0.03), y = c(0.96, 0.64), gp = gpar(cex = 1.2))
dev.off()


# UNSUPERVISED CLUSTERING OF GENES BY SIMILARITY  ++++++++++++++++++++++++++++++
#
# for heatmaps and clustering,
# a matrix must be generated, using wide format and relative fold change
# first we paste conditions together (will be columns) 
mat <- unite(df, condition, condition, induction, timepoint) %>%
  
  # and sgRNA indices (will be rows) to help spreading
  unite(sgRNA, sgRNA_short, sgRNA_index) %>%
  
  # select only the required columns
  select(sgRNA, condition, log2FoldChange) %>%
  
  # replace NA in log2FoldChange with 1
  mutate(log2FoldChange = coalesce(log2FoldChange, 0)) %>%
  
  # and spread with condition as key and log2foldchange as value
  spread(key = condition, value = log2FoldChange) %>%
  
  # rename columns clumsily
  rename_at(vars(ends_with("_0")), funs(sub("_0$", "_00", .))) %>%
  rename_at(vars(ends_with("_1")), funs(sub("_1$", "_01", .))) %>%
  rename_at(vars(ends_with("_2")), funs(sub("_2$", "_02", .))) %>%
  rename_at(vars(ends_with("_4")), funs(sub("_4$", "_04", .))) %>%
  rename_at(vars(ends_with("_8")), funs(sub("_8$", "_08", .)))

# change row names
rownames(mat) <- mat$sgRNA

# reorder a bit and coerce to matrix
mat <- ungroup(mat) %>% select(-sgRNA) %>%
  select(contains("00"), contains("01"), contains("02"), contains("04"), 
         contains("08"), contains("16"), contains("32")) %>%
  select(contains("LL"), contains("HL"), contains("DN"), contains("LAC"), contains("NACL")) %>%
  select(contains("_u_"), contains("_i_")) %>% as.matrix
  
# optionally select only a subset of the data, e.g. no LAC/NACL
mat <- mat[, grepl("HL|LL|DN", colnames(mat))]

# compute dissimilarity matrix for genes and cluster using hclust
# for the clustering algorithm, see hclust manual
cluster <- hclust(dist(mat), method = "ward.D2")


# CHOOSE AND PLOT MAJOR CLUSTERS  ++++++++++++++++++++++++++++++++++++++++++++++
#
# Perform iterative silhouette analysis for a range of cluster numbers
# Prerequisite is only a cluster object (such as obtained from hclust())
silhouetteResult <- silhouetteAnalysis(mat, cluster, 2:20)


# plot the results from silhouette analysis
svg("CRISPRi_library_average_silhouette.svg", 9, 4.5)
print(silhouetteResult$plot.clusters, position = c(0, 0, 0.6, 1), more = TRUE)
print(silhouetteResult$plot.summary, position = c(0.6, 0, 1, 1))
dev.off()


# We can select the k upmost clusters using cutree. 
clusters <- vegan::cutreeord(cluster, k = 6); clusters %>% table
# combine clusters 4 and 5 in 5 as it is artificially separated; rename cluster 6 to 4.
clusters <- clusters %>% recode(`4` = 5L, `6` = 4L) %>% setNames(., names(clusters))
# The cluster order by cutree does not correspond to the one of the dendrogram
# therefore change order of clusters manually using factor levels
clusters <- factor(clusters, unique(clusters) %>% sort)
# make the silhouette for plotting
sil <- silhouette(as.numeric(clusters), dist(mat))
# merge predicted clusters with df
df <- group_by(df, condition, timepoint, induction) %>% mutate(cluster = clusters)
# custom color vector for clusters, without first color that equals last
#colorspace::hclwizard()
custom_palette <- colorspace::qualitative_hcl(n = length(unique(clusters)), h = c(30, 360), c = 100, l = 50)
custom_lattice <- custom.lattice()
custom_lattice$superpose.polygon$col = custom_palette
custom_lattice$superpose.symbol$col = custom_palette
custom_lattice$superpose.line$col = custom_palette


# re-order the dendrogram a bit for plotting, in order of clusters
cluster_reordered <- as.dendrogram(cluster) %>% 
  reorder(wts = clusters[cluster$labels], agglo.FUN = mean)


# plot colored dendrogram, heatmap and clusters all on one page
library(grid)
library(gridBase)

#svg("CRISPRi_library_heatmap.svg", width = 12, height = 7)
png("CRISPRi_library_heatmap.png", width = 1750, height = 800, res = 165)
plot.new()
viewport(x = 0.01, y = 0.21, width = 0.98, height = 0.79, just = c("left", "bottom")) %>%
  pushViewport

# plot dendrogram first
par(new = TRUE, fig = gridFIG())
plot(color_branches(
  cluster_reordered,
  #plyr::mutate(cluster_reordered, labels = rep("", length(labels))),
  #k = length(unique(clusters)),
  groupLabels = as.numeric(levels(clusters)), lwd = 1,
  col = custom_palette[as.numeric(levels(clusters))],
  clusters = sort(as.numeric(clusters))
))

upViewport()
viewport(x = 0, y = 0, width = 1, height = 0.35, just = c("left", "bottom")) %>%
  pushViewport

# plot heatmap according to clustered protein groups
lp <- levelplot(
  # limit plotted range to log2 FC of -5 to 5 (32 fold)
  mat[as.hclust(cluster_reordered)$order, ] %>% replace(., . > 5, 5) %>% replace(., . < -5, -5),
  col.regions = colorRampPalette(c(custom_palette[4], grey(0.95), custom_palette[1])),
  as.table = TRUE, aspect = "fill", 
  xlab = "sgRNA", ylab = "time [d]",
  scales = list(alternating = FALSE, x = list(draw = FALSE))
)
print(lp, more = TRUE, position=c(0.017, 0, 1, 1))
dev.off()


# plot colored dendrogram and silhouette width
png("CRISPRi_library_silhouette.png")
plot(sil, 
  border = NA, 
  col = custom_palette[sort(unique(clusters))])
dev.off()


# REGULATION FOR DIFFERENT FUNCTIONAL GROUPS   +++++++++++++++++++++++++++++++++
#
# general overview about clusters
plot_library_depletion <- xyplot(log2FoldChange ~ generations | paste(condition, " - ", cluster),
  # we have to add one more fake data point, at max of 50 generations, with NA
  # so that lines are broken after last time point
  bind_rows(
    df %>% ungroup %>% filter(induction == "i", condition %in% c("DN", "HL", "LL")),
    df %>% ungroup %>% filter(induction == "i", condition %in% c("DN", "HL", "LL"), timepoint == 0) %>%
      mutate(generations = 50, log2FoldChange = NA)
  ) %>% arrange(condition, cluster, sgRNA, generations),
  groups = cluster, type = "l", between = list(x = 0.6, y = 0.6), 
  as.table = TRUE, layout = c(length(unique(clusters)), 3),
  col = paste0(custom_lattice$superpose.polygon$col, "20"),
  xlab = "generations", ylab = "log2 FC", 
  ylim = c(-8, 4), xlim = c(0, 40),
  scales = list(alternating = FALSE), pch = 16, cex = 0.7, 
  par.settings = custom_lattice,
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.xyplot(x, y, ...)
    #panel.stripplot(x, y, jitter.data = TRUE, amount = 0.2, horizontal = FALSE, ...)
    medianFC <- tapply(y, x, function(x) median(x, na.rm=TRUE))
    panel.xyplot(names(medianFC) %>% as.numeric, medianFC, col = grey(0.3), type = "l", lwd = 3)
    panel.text(10, 2.5, paste(length(x)/length(unique(x)), "sgRNAs"), cex = 0.6, col = grey(0.5))
  }
)


# GENE SET ENRICHMENT USING TOPGO  +++++++++++++++++++++++++++++++++++++++++++++

# Using GO enrichemnt, TopGO package by Alexa et al.
library(topGO)

# the GetTopGO function takes as input either
#   1. a data frame with one row per gene as only presence in cluster is 
#      important here. Three columns are obligatory, "cluster", "GeneID", and 
#      "Gene.ontology.IDs" with GO terms separated by '; '
#   2. alternatively three vectors or lists corresponding to
#      the three columns can be provided

# reduce table to one sgRNA per row and call GetTopGO function
# execute GetTopGO function and collect list of results
TopGoResult <- lapply(1:length(unique(clusters)), function(i) {
    with(df[!duplicated(df$sgRNA), ] %>% as.data.frame,
      GetTopGO(df = NULL,
        cluster = as.numeric(cluster),
        GeneID = sgRNA,
        Gene.ontology.IDs,
        topNodes = 50, 
        selected.cluster = i
      )
    )
  }) %>%

  # regarding p-value adjustment: the authors discourage from multiple 
  # hypothesis testing and indeed, it turns most p-values insignificant
  setNames(1:length(unique(clusters))) %>%
    
  # turn into data frame with cluster as ID columns
  plyr::ldply(., .id = "cluster")
    
  # save unfiltered TopGO result
  write_csv(TopGoResult, "TopGoResult.csv")
  # or load existing one
  #TopGoResult <- read_csv("TopGoResult.csv")
  

# ADDITIONAL FILTERING
# optional filtering by REVIGO that excludes redundant GO terms
# To this end, the resulting terms have to be uploaded to http://revigo.irb.hr/
# and the trimmed list be downloaded. 

# merge TopGo result with REVIGO dispensibility analysis
TopGoResultFiltered <- left_join(as_tibble(TopGoResult), read_csv("REVIGO.csv"), 
  by = c("GO.ID" = "term_ID")) %>%
  
  # trim data frame by p.value and GO term node size
  filter(dispensability <= 0.5,
    elimFisher < 0.03,
    Annotated >= 5,
    Annotated <= 200) %>%
  
  # optionally trim lengthy GO terms and fuse with cluster (some show up several times)
  mutate(
    TermAbbr = gsub("(biosynthetic|metabolic|catabolic) proce?s?s?|\\.\\.\\.|metabolic.?.?.?", 
    replacement = "", Term) %>% paste0(., "...", cluster)
  ) %>% filter(!duplicated(TermAbbr))


# Plot results of GO enrichment for each cluster
# (p-value is probability that differential expression of interesting genes
# is identical to background)
plot_TopGO_table <- xyplot(factor(TermAbbr, unique(TermAbbr)) ~ as.numeric(elimFisher), 
  TopGoResultFiltered %>% filter(cluster != 5) %>% arrange(desc(cluster), desc(elimFisher)), 
  groups = factor(cluster), main = "", #"GO enrichment, GO terms by p-value",
  par.settings = custom_lattice, xlab = "p-value", ylab = "GO term",
  scales = list(alternating=FALSE, x = list(rot = 45)), 
  panel = function(x, y, groups = groups, subscripts = subscripts, ...) {
    panel.abline(h = 1:length(y), col = grey(0.9), lty = 2)
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.key(labels = paste("cluster", 1:4), pch = 19, corner = c(0.9, 0.95))
    panel.superpose(x, y, groups, subscripts, ...)
  },
  panel.groups = function(x, y, groups = groups, subscripts = subscripts, ...) {
    panel.xyplot(x, y, cex = 2.5-(50*x), pch = 19, alpha = 0.7,
      col = custom_lattice$superpose.polygon$col[subscripts])
  }
)

# Plot results of GO enrichment as a map of semantic similarity
plot_TopGO_map <- xyplot(plot_Y ~ plot_X, 
  TopGoResultFiltered %>%filter(cluster != 5),
  cluster = TopGoResultFiltered$cluster,
  bubblesize = 3+0.08 * sapply(TopGoResultFiltered$SigGenes, function(x) {
    strsplit(x, split = ",") %>% unlist %>% length}), 
  xlim = c(-7.5, 9.5), ylim = c(-7.5, 8.5),
  main = "GO enrichment, GO terms by semantic similarity (p<0.025)",
  par.settings = custom_lattice, xlab = "semantic space X", ylab = "semantic space Y",
  panel = function(x, y, bubblesize, cluster, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.key(labels = paste("cluster", 1:4), pch = 19, corner = c(0.9, 0.95))
    panel.key(labels = paste(c(10, 20, 40), "sgRNAs\n"), points = FALSE, col = grey(0.6), corner = c(0.01, 0.95))
    panel.xyplot(x = -4.1, y = c(7.6, 6.7, 5.8), col = grey(0.6), pch = 19, cex = 3+0.08*c(10, 20, 40))
    panel.xyplot(x, y, alpha = 0.5, cex = bubblesize, col = custom_palette[cluster], pch = 19)
    selected = TopGoResultFiltered$dispensability <= 0.5
    panel.text(x[selected], y[selected], cex = 0.8, col = grey(0.4),
      labels = TopGoResultFiltered$TermAbbr[selected])
  }
)


# PLOT ENRICHMENT / DEPLETION PER FUNCTIONAL GROUP  ++++++++++++++++++++++++++++
#
# investigating different functional groups that were enriched for each cluster
# generalized plotting function for selected genes/sgRNAs

plot.sgRNAs <- function(data, comment = NULL) {
  
  xyplot(log2FoldChange ~ factor(timepoint), 
    data, as.table=TRUE, groups = cluster,
    xlab = "time [d]", ylab = "log2 FC", 
    between = list(x = 0.5, y = 0.5),
    layout = c(5, 4),
    ylim = c(-9, 3), pch = 19, cex = 0.2, lwd = 2, 
    par.settings = custom_lattice,
    scales = list(alternating = FALSE),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.superpose(x, y, ...)
      panel.text(3, 2, labels = comment, cex = 0.6, col = grey(0.5))
    },
    panel.groups = function(x, y, ...) {
      panel.stripplot(x, y, jitter.data = TRUE, amount = 0.2, horizontal = FALSE, ...)
      panel.smoother(x, y, ...)
    }
  )
}

# make plots with genes for the most enriched GO terms
plot.list <- filter(TopGoResultFiltered, cluster != 5) %>%
  
  # select top 4 terms per cluster by p-value
  group_by(cluster) %>% slice(1:4) %>%
  
  # pass only the abbreviated GO term/cluster combination
  pull(TermAbbr) %>% lapply(., function(x) {
  
    # extract sgRNA or gene IDs
    sgRNAs = filter(TopGoResultFiltered, TermAbbr %in% x)$SigGenes %>%
      strsplit(split = ",") %>%
      unlist
  
    # generate plot of subsetted data
    plot.sgRNAs(
      data = filter(df, sgRNA %in% sgRNAs, induction == "i") %>%
        mutate(log2FoldChange = coalesce(log2FoldChange, 0)),
      comment = paste(length(sgRNAs), "sgRNAs")
    )
  }) %>%
  
  setNames(filter(TopGoResultFiltered, cluster != 5) %>% 
    group_by(cluster) %>% slice(1:4) %>% 
    pull(TermAbbr) %>% substr(1, 15))


#png("CRISPRi_library_TopGoResult.png", width = 1400, height = 800, res = 110)
svg("CRISPRi_library_TopGoResult.svg", width = 10, height = 11)
print(plot_library_depletion, position = c(0, 0.54, 1, 1), more = TRUE)
print(plot_TopGO_table, position = c(0, 0, 0.5, 0.56), more = TRUE)
print(do.call(c, plot.list), position = c(0.5, 0, 1, 0.56))
dev.off()


# EXPLORING DIFFERENCES BETWEEN CONDITIONS ON GENE LEVEL +++++++++++++++++++++++
#
# here the fitness score can be the major determinant to find differences between
# two (or more?) selected conditions, and clusters
# first have a look at global fitness score distribution
plot_global_fitness <- df %>% 
  
  # select conditions for comparison
  filter(condition %in% c("LL", "HL"), induction == "i", 
    timepoint == 0, cluster != 5) %>% 
  
  # create boxplot or violinplot for fitness score, direct comparison
  xyplot(fitness_score ~ factor(condition, c("LL", "HL")) | factor(cluster), .,
    as.table = TRUE, par.settings = custom_lattice,
    groups = cluster, between = list(x = 0.5, y = 0.5),
    scales = list(alternating = FALSE), 
    #main = "global fitness score distribution, LL vs HL",
    xlab = "", ylab = "fitness score F",
    ylim = c(-8, 3), 
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.stripplot(as.numeric(x)-0.25, y, jitter.data = TRUE, amount = 0.14, 
        pch = ".", cex = 2, horizontal = FALSE, ...)
      panel.violin(as.numeric(x)+0.25, y, horizontal = FALSE, box.ratio = 0.5, 
        border = grey(0.5), col = custom_palette[panel.number()], ...)
      panel.pvalue(x, y, std = "LL", offset = 0, fixed.pos = 0.5, 
        col = custom_palette[panel.number()], ...) 
    }
  )

# second step is to have a look at strongest difference between sgRNAs 
# of two conditions. First create a list of candidates, select conditions for comparison
df_selected_fitness <- df %>% 
  
  # select conditions for comparison
  filter(condition %in% c("LL", "HL"), induction == "i", 
    timepoint == 0, cluster != 5) %>%
  
  # optional filter for sgRNAs that have 2 members per cluster
  group_by(cluster, sgRNA_short) %>% filter(length(sgRNA_index) == 4) %>% 
  
  # filter by highest fitness diff between sgRNA of two conditions (HL-LL)
  # positive for stronger HL depletion, negative for stronger LL depletion
  group_by(cluster, sgRNA_short) %>% summarise(
    diff_fitness_score = mean(fitness_score[c(1,3)]) - mean(fitness_score[c(2,4)]),
    Process.abbr = unique(Process.abbr)
  ) %>% 
  
  # group by cluster and sort by fitness score for top n filtering
  group_by(cluster) %>% arrange(cluster, desc(diff_fitness_score)) 


# plot differential fitness score per gene, per cluster and Process
plot_diff_fitness <- xyplot(diff_fitness_score ~ factor(Process.abbr) | factor(cluster), 
    df_selected_fitness, group = between(diff_fitness_score, -3, 3),
    #main = "differential fitness score of genes, HL vs DN",
    xlab = "", ylab = "differential fitness score (dF = F_HL - F_LL)",
    par.settings = custom_lattice, type = "p", pch = 19, as.table = TRUE, 
    col = c("#D33F6A", "#3D8900"),
    scales = list(alternating = FALSE, x = list(rot = 45)), 
    ylim = c(-6, 6), between = list(x = 0.5, y = 0.5),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.abline(h = c(-3, 3), lty = 2, lwd = 1, col = grey(0.5))
      panel.stripplot(x, y, jitter.data = TRUE, amount = 0.2, 
        alpha = 0.6, horizontal = FALSE, ...)
      panel.text(1, 4, pos = 4, col = grey(0.5),
        labels = paste0(sum(y <= -3 | y >= 3, na.rm = TRUE), " out of ", length(y), " genes"))
    }
  )


# fish out most different ones, pick genes above certain threshold
# plot the log2FC over generations for the top differential genes
custom_lattice$superpose.polygon$col = c("#3D8900", "#3D8900", "#D33F6A", "#D33F6A")
plot_top_diff_fitness <- xyplot(log2FoldChange ~ generations | paste0(cluster, "   ", sgRNA_short), 
    df %>% filter(
      sgRNA_short %in% filter(df_selected_fitness, 
        !between(diff_fitness_score, -3, 3))$sgRNA_short, 
      induction == "i", condition %in% c("LL", "HL")), 
    #main = "genes with highest differential fitness score, HL vs DN",
    as.table = TRUE, groups = paste(condition, sgRNA_index) %>% factor(.), 
    par.settings = custom_lattice, type = "l", 
    scales = list(alternating = FALSE, cex = 0.8),
    ylim = c(-10, 10), xlim = c(0, 17),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.key(c("HL", "LL"), points = FALSE, col = c("#3D8900", "#D33F6A"),
        lines = TRUE, cex = 0.7, corner = c(0.05, 0.95), ...)
      panel.abline(h = 0, lty = 2, lwd = 1.5, col = grey(0.5))
      #panel.xyplot(x, y, ...)
      panel.xyarea(x, y, origin = 0, alpha = 0.6, ...)
    }
  )


# plot diff fitness score investigation on one page
#png("CRISPRi_library_fitness_score_HL_LL.png", width = 2000, height = 1300, res = 130)
svg("CRISPRi_library_fitness_score_HL_LL.svg", width = 12, height = 7)
print(plot_global_fitness, position = c(0, 0.5, 0.45, 1), more = TRUE)
print(plot_diff_fitness, position = c(0, -0.02, 0.45, 0.55), more = TRUE)
print(plot_top_diff_fitness, position = c(0.44, 0.03, 1, 1))
dev.off()


# SIMPLE FITNESS SCORE ANALYSIS LIKE IN S. GOLDEN PAPER
# 
# plot histogram of fitness scores for all genes
# select one timepoint is enough
df %>% filter(timepoint == 0, !is.na(fitness_score)) %>%
  
  histogram( ~ fitness_score | induction * condition, .,
    main = "distribution of average fitness score per gene",
    as.table = TRUE, nint = 50,
    par.settings = custom.lattice,
    scales = list(alternating = FALSE, cex = 0.8),
    xlim = c(-5, 5), 
    panel = function(x, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.histogram(x, alpha = 0.5, lwd = 1.5, ...)
      panel.abline(v = -1.5, col = grey(0.5), lty = 2)
      panel.text(2.5, 50, col = grey(0.5),
        labels = paste0(round(sum(x <= -1.5)/length(x)*100, 1), "% sgRNAs F <= -1.5"))
    }
  )


# CORRELATION OF SGRNA DEPLETION AND PROTEIN MASS ++++++++++++++++++++++++++++++
#
# first we need to think about mapping conditions to each other.
# the MS dataset does not contain DN cycle but HL and LL conditions

# We can map HL and LL to 300 and 100 µE respectively (this is what was used by 
# Lun Yao who did the experiments)
df <- read_csv("../20171204_MS_CO2/CO2_light/analysis/diffacto_weightedsum/df_long_massfrac.csv") %>%
  
  # load mass fraction data from Cell Reports paper
  filter(sample == "light", light %in% c(300, 100)) %>%
  
  # change sample annotation from [100, 300] to [LL, HL]
  mutate(condition = recode(light, `300` = "HL", `100` = "LL")) %>%
  
  # rename a column
  rename(locus = protein) %>%
  
  # select only the interesting columns
  select(locus, condition, mean.mass.fraction.norm) %>%
  
  # merge with df
  left_join(df, .) 
  
  # test if mass fraction sums up to _around 1_ (also depends on intersection
  # between proteomics and CRISPRi library gene sets)
  df %>% group_by(condition, timepoint, induction, sgRNA_index) %>% 
    summarise(sum(mean.mass.fraction.norm, na.rm=TRUE)) %>%
    print(n = Inf)

  
# Next step: Check how many sgRNA pairs end up in same cluster
# for this, all we need is one particular condition
sgRNA_plot1 <- df %>% filter(condition == "HL", timepoint == 16, induction == "i") %>%
  
  # now we need to group by cluster and locusTag
  group_by(cluster) %>% 
  
  # and summarise the number of unique sgRNAs for cluster and locusTag
  summarise(
    n_sgRNA_1 = sum(table(sgRNA_short) == 1), 
    n_sgRNA_2 = sum(table(sgRNA_short) == 2),
    n_sgRNA_total = length(unique(sgRNA_short))
  ) %>%
  
  # plot a barchart out of that. There is no interesting correlation, it's clearly related
  # to cluster size, meaning the probability for 2 sgRNAs to land in a bigger cluster is 
  # higher (poisson distribution)
  barchart(n_sgRNA_1/n_sgRNA_total * 100 + n_sgRNA_2/n_sgRNA_total * 100 ~ cluster, .,
    as.table = TRUE, par.settings = custom.lattice,
    main = "co-appearance of sgRNA1 and 2 per cluster",
    xlab = "cluster", ylab = "% 1 or 2 sgRNAs in cluster",
    ylim = c(-2, 122),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.barchart(x, y, lwd = 0, ...)
      panel.key(labels = c("1 sgRNA", "2 sgRNAs"), pch = 15)
    }
  )

# We can also check if sgRNA 1 or 2 have induce stronger depletion
sgRNA_plot2 <- df %>% filter(induction == "i", timepoint == 16, 
  condition %in% c("HL", "LL", "DN")) %>%
  
  # plot log2FoldChange per cluster, sgRNA_index, and conditions
  densityplot( ~ log2FoldChange | cluster * condition, .,
    as.table=TRUE, par.settings = custom.lattice,
    groups = sgRNA_index, scales = list(alternating = FALSE),
    main = "comparison of sgRNA1 and 2 at 16d, per cluster",
    xlab = "log2 FC", ylab = "frequency",
    xlim = c(-8, 3), pch = ".", lwd = 2,
    panel = function(x, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.abline(v = 0, col.line = grey(0.7), lty = 2, lwd = 1.5)
      panel.densityplot(x, ...)
      panel.key(labels = c("sgRNA 1", "sgRNA 2"), 
        points = FALSE, lines = TRUE, lwd = 2)
    }
  )


# Investigating the mass fraction of genes broken down by cluster
# we can plot the cumulative mass fraction of proteins
sgRNA_plot3 <- df %>% filter(condition %in% c("LL", "HL"), induction == "i", timepoint == 16) %>%
  
  # plot log2FoldChange per cluster, sgRNA_index, and conditions
  xyplot(mean.mass.fraction.norm ~ 1:length(sgRNA) | cluster * condition, .,
    as.table = TRUE, par.settings = custom.lattice,
    groups = sgRNA_index, type = "p", pch = 19, cex = 0.8,
    scales = list(alternating = FALSE), main = "cumulative mass fraction per cluster",
    xlab = "gene number", ylab = "cumulative protein mass fraction",
    ylim = c(0, 0.65), xlim = c(-100, 1800),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      y = cumsum(sort(na.omit(y), decreasing = TRUE))
      panel.xyplot(1:length(y), y, ...)
      #panel.abline(v = 0, col.line = grey(0.7), lty = 2, lwd = 1.5)
      panel.key(labels = c("sgRNA 1", "sgRNA 2"), pch = 19)
      #panel.text(-6, -1.3, labels = paste0(length(y), " sgRNAs"), col = grey(0.5))
    }
  )

png("CRISPRi_library_sgRNA_1vs2.png", width = 1100, height = 1000, res = 110)
print(sgRNA_plot1, position = c(0, 0.48, 0.4, 1), more = TRUE)
print(sgRNA_plot2, position = c(0.4, 0.48, 1, 1), more = TRUE)
print(sgRNA_plot3, position = c(0, 0, 1, 0.5))
#grid.text(label = c("A", "B", "C"), x = c(0.03, 0.43, 0.03), y = c(0.96, 0.96, 0.48), gp = gpar(cex = 1.3))
dev.off()


# STRING-DB TO EXPLORE GENE NETWORKS WITHIN CLUSTERS +++++++++++++++++++++++++++

# Exploring the correlation between sgRNA log2FC and network connections
# we can try with acyl carrier protein, GeneID = ssl2084
#read_tsv("https://string-db.org/api/tsv/network?identifier=ssl2084&species=1148")

# This function retrieves the network for an arbitrary list of genes 
# from our data, e.g. for all clusters
get_String_db <- function(geneList) {
  
  # trim to only unique terms
  unique(geneList) %>%
  # paste together for one single query
  paste(collapse = "%0D") %>%
  # channel further to STRINGdb API
  paste0("https://string-db.org/api/tsv/network?identifiers=", ., "&species=1148&required_score=900") %>%
  # retrieve the query result
  read_tsv(.)

}

# retrieve networks for the first 4 clusters
df_string <- lapply(c(3,4,5,6), function(i) {
  ungroup(df) %>% filter(cluster == i) %>% 
  pull(locus) %>% unique %>% na.omit %>%
  get_String_db
}) %>% bind_rows(.id = "cluster") %>% 
  mutate(cluster = as.numeric(cluster) + 2)


# retrieve networks for the large clusters that have to be processed in several 
# chunks, as one limitation with the API is that the URI can not exceed a certain length
df_string <- lapply(1:4, function(i) {
  ungroup(df) %>% filter(cluster == 1, !duplicated(locus)) %>% 
  mutate(chunk = rep(1:4, length.out = nrow(.)) %>% sort) %>% 
  filter(chunk == i) %>%
  select(locus) %>% unlist %>%
  get_String_db
}) %>% bind_rows %>% mutate(cluster = 1) %>% bind_rows(df_string, .)




# Analyze the network retrieved from STRINGdb. Gather the two interaction partners
# into one column and add annotation
df_string_summary <- df_string %>% gather(key = interaction_partner, 
  value = GeneID, preferredName_A:preferredName_B) %>%
  
  # Group by cluster and interaction partner
  group_by(cluster, GeneID) %>%
  
  # summarise the number of significant interactions
  summarise(interactions = length(score), mean_score = mean(score)) %>%
  
  # sort by interactions
  arrange(cluster, desc(interactions)) %>%
  
  # add a rank argument
  mutate(rank = 1:length(cluster)) %>% 
  
  # add pathway annotation from original df
  # first convert STRINGdb IDs back to locusTags
  mutate(Gene.names = sapply(GeneID, function(x) {
    grep(x, value = TRUE,
      ungroup(df) %>% filter(!duplicated(Gene.names)) %>%
      select(Gene.names) %>% unlist
    )[1]
  })) %>%
  
  left_join(., 
    ungroup(df) %>% filter(!duplicated(Gene.names)) %>% 
    select(Gene.names, locus, Process, Pathway, Process.abbr, Pathway.abbr)
  )


# also add number of interactions to each gene in main data frame
df <- df_string_summary %>% 
  
  # only select important columns from STRING summary
  select(cluster, Gene.names, interactions, mean_score) %>%
  
  # deal with duplicated entries: map interactions for the same gene in 2 different
  # cluster contexts (2 sgRNAs!) using cluster grouping
  filter(!is.na(Gene.names)) %>% left_join(df, .) %>%
  
  # turn character into factor
  mutate(cluster = as.numeric(cluster))


# plot interactions per cluster and gene
inter_plot1 <- xyplot(interactions ~ rank | paste("cluster", cluster), df_string_summary,
  as.table=TRUE, par.settings = custom.lattice, 
  layout = c(length(unique(df_string_summary$cluster)), 1),
  ylim = c(0, 40), ylab = "number of interactions",
  groups = cluster, scales = list(alternating = FALSE),
  main = "Number of network interaction partners (STRINGdb, score > 0.9)",
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.horizonplot(x, y, col.regions = custom.lattice$superpose.polygon$col[panel.number()], ...)
    panel.ablineq(h = mean(y), label = paste("mean interactions:", round(mean(y), 1)), 
      pos = 3, fontfamily = "FreeSans", col = grey(0.3), lty =2)
  }
)

# plot top 12 interactions per cluster, exclusively ribosomal proteins
inter_plot2 <- xyplot(interactions ~ factor(GeneID) | paste("cluster", cluster), 
  df_string_summary %>% filter(rank <= 12),
  as.table=TRUE, par.settings = custom.lattice,
  ylim = c(0, 50), xlim = c(0, 13), 
  layout = c(length(unique(df_string_summary$cluster)), 1),
  ylab = "number of interactions", xlab = "rank",
  groups = cluster, scales = list(alternating = FALSE),
  main = "Top 12 interaction partners with ribosomal proteins",
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.barplot(1:12, y, ewidth = 0.3, border = NA, custom.lattice$superpose.polygon$col[panel.number()])
    panel.text((1:12) + 0.5, y+3, pos = 3, labels = x, srt = 45, col =grey (0.4), cex = 0.8)
  }
)

# plot top 10 interactions per cluster, WITHOUT ribosomal proteins
inter_plot3 <- xyplot(interactions ~ factor(GeneID) | paste("cluster", cluster), 
  #filter ribosomal proteins and re-calculate ranks
  df_string_summary %>% 
    filter(!grepl("^rp[msl]", GeneID)) %>% 
    mutate(rank = 1:length(cluster)) %>% 
    filter(rank <= 12),
  as.table=TRUE, par.settings = custom.lattice,
  ylim = c(0, 50), xlim = c(0, 13), 
  layout = c(length(unique(df_string_summary$cluster)), 1),
  ylab = "number of interactions", xlab = "rank",
  groups = cluster, scales = list(alternating = FALSE),
  main = "Top 12 interaction partners without ribosomal proteins",
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.barplot(1:12, y, ewidth = 0.3, border = NA, custom.lattice$superpose.polygon$col[panel.number()])
    panel.text((1:12) + 0.5, y+3, pos = 3, labels = x, srt = 45, cex = 0.8, 
      col = ifelse(grepl("ndh|psa|psb|mur|mra", x), "#D13636", grey(0.4)))
  }
)


png("CRISPRi_library_network_interactions.png", width = 1000, height = 900, res = 100)
print(inter_plot1, position = c(0, 0.66, 1, 1), more = TRUE)
print(inter_plot2, position = c(0, 0.33, 1, 0.66), more = TRUE)
print(inter_plot3, position = c(0, 0.00, 1, 0.33))
grid.text(label = c("A", "B", "C"), x = c(0.03, 0.03, 0.03), y = c(0.96, 0.64, 0.31), gp = gpar(cex = 1.3))
dev.off()


# plotting a network from a list of edges and nodes using igraph/tidygraph packages
# input are two lists or data frames, 'nodes' with unique name and ID of genes
# and 'edges' indicating all recorded connections between nodes; However,
# the information in 'edges' already contains all possible nodes, so 'edges' alone can be sufficient
library(tidygraph)
library(ggplot2)
library(ggraph)


# construct main graph object containing data for all clusters
network <- lapply(c(3:6), function(i) {df_string %>% filter(cluster == i) %>%
  
  # select columns with Node A, Node B, and optional score or weight
  select(preferredName_A, preferredName_B, score) %>%

  # consctruct igraph object (tbl_graph is just an igraph wrapper)
  tbl_graph(
    edges = .,
    nodes = df_string_summary %>% ungroup %>% mutate(Process.abbr = factor(Process.abbr)) %>% filter(cluster == i) %>%
      select(GeneID, interactions, cluster, rank, mean_score, Process.abbr)
  ) %>%
  
  # filter out connections with just a few genes
  filter(interactions >= 4) %>%
  
  # plot the graph
  # the layout argument takes one of: 'star', 'circle', 'gem', 'dh', 
  # 'randomly', 'fr', 'kk', 'drl', 'lgl', 'graphopt', 'grid', 'mds'
  ggraph(., layout = 'nicely') +
  geom_edge_link(colour = grey(0.5)) + 
  geom_node_point(aes(colour = Process.abbr, size = interactions), show.legend = FALSE) +
  geom_node_text(aes(label = GeneID), repel = TRUE, col = grey(0.5), cex = 3) +
  theme_graph(background = grey(0.95), foreground = grey(0.5),
    plot_margin = margin(10, 10, 10, 10)) +
  facet_grid(~ cluster) +
  scale_colour_manual(values = colorRampPalette(custom_palette)(17), 
    name = "Process", drop = FALSE)
  
})

# make a custom legend and 
legend <- ggplot(df_string_summary, aes(interactions, GeneID)) +
    geom_point(aes(colour = Process.abbr)); legend = cowplot::get_legend(legend)


png("CRISPRi_library_network_topology.png", width = 2500, height = 500, res = 120)
do.call(gridExtra::grid.arrange, list(grobs = {network[[5]] = legend; network}, ncol = 5))
dev.off()



# SAVE PROCESSED DATA FOR RE-IMPORT  +++++++++++++++++++++++++++++++++++++++++++
# save data frame including clusters
save(df, file = "processed_data/CRISPRi_library_df_clusters.Rdata")
write_csv(df, "processed_data/CRISPRi_library_df_clusters.csv")
write_csv(ungroup(df) %>% filter(!duplicated(sgRNA)), "processed_data/CRISPRi_library_df_clusters_only.csv")
save(df, file = "/media/ProteomicsCyano/Resources/treemaps/Synechocystis_CRISPRi_library/CRISPRi_library_df_clusters.Rdata")
