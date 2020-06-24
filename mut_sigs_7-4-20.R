# load('~/Desktop/script_test/deconstructSigs/data/signatures.genome.cosmic.v3.may2019.rda')


excluded_samples <- c("B241R41-2",  "A373R7", "A512R17")
excluded_samples <- c(excluded_samples, "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R09", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
excluded_samples <- c(excluded_samples, "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
excluded_samples <- c(excluded_samples, "D106R1", "D106R3", "D106R5", "D106R7", "D106R9", "D106R11", "D106R13", "D106R15", "D106R17", "D106R19", "D106R21", "D106R23", "D106R25", "D106R27", "D106R29", "D106R31", "D106R33"  )
excluded_samples <- c(excluded_samples, "D197R09", "D197R11", "D197R13", "D197R15")
excluded_samples <- c(excluded_samples, "D265R01", "D265R03", "D265R05", "D265R07", "D265R09", "D265R11", "D265R13")


####### mutationalPatterns #####
library(BSgenome)
library("gridExtra")


ref_genome <- ("BSgenome.Dmelanogaster.UCSC.dm6")
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene", version = "3.8")

library(ref_genome, character.only = TRUE)
library(MutationalPatterns)

vcf_files <- list.files(pattern = "_consensus_filt.vcf", full.names = TRUE)
# Remove excluded samples
vcf_files <- grep(vcf_files, pattern = paste(excluded_samples, collapse = '|'), value = T, invert = T)

sample_names <- basename(vcf_files)
sample_names <- gsub("_.*", "", sample_names)

vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

# vcfs <- read_vcfs_as_granges(vcf_files, 'combined', ref_genome)

tissue <- rep("Male", length(sample_names))

muts <- mutations_from_vcf(vcfs[[1]])
types <- mut_type(vcfs[[1]])
context <- mut_context(vcfs[[1]], ref_genome)
type_context <- type_context(vcfs[[1]], ref_genome)


type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

## Rainfall
# chromosomes <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "X", "Y")
# plot_rainfall(vcfs[[1]], chromosomes,
#               title = "", cex = 2.5,
#               cex_text = 3, ylim = 1e+08
# )

# # plot spectrum
# p1 = plot_spectrum(type_occurrences)
# p1
#

# # Plot the mutation spectrum with distinction between C>T at CpG sites:
p2 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)
p2
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)
p4
# # Plot spectrum without legend:
# p3 = plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE)
# grid.arrange(p1, p2, p3, ncol=3, widths=c(2,2,3))

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
# plot_96_profile(mut_mat[, c(1, 5)], condensed = TRUE)
# plot_96_profile(mut_mat, condensed = TRUE)

# Download signatures from pan-cancer study (Alexandrov et al., 2013):

# sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")

cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)

new_order <- match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cancer_signatures <- cancer_signatures[as.vector(new_order), ]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) <- cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures <- as.matrix(cancer_signatures[, 4:33])

# Plot some signatures
# plot_96_profile(cancer_signatures[, c("Signature.15", "Signature.12", "Signature.25")], condensed = TRUE, ymax = 0.3)

# Do some similarity clustering between signatures
hclust_cosmic <- cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order <- colnames(cancer_signatures)[hclust_cosmic$order]
# plot(hclust_cosmic)


# Similarity between mutational profiles and COSMIC signatures
# The cosine similarity re- flects how well each mutational profile can be explained by each signature individually.

# for (n in colnames(mut_mat)) {
#   print(paste(n, round(cos_sim(mut_mat[, n], cancer_signatures[, 14]), 2)))
# }
# plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = FALSE, mode = "absolute")

# Calculate pairwise cosine similarity between mutational profiles and COSMIC signatures
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, cancer_signatures)
plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = TRUE)

fit_res <- fit_to_signatures(mut_mat, cancer_signatures)

# Select signatures with some contribution
select <- which(rowSums(fit_res$contribution) > 30)
# Plot contribution barplot
plot_contribution(fit_res$contribution[select, ],
                  cancer_signatures[, select],
                  coord_flip = FALSE,
                  mode = "absolute"
)

# Plot relative contribution of the cancer signatures in each sample as a heatmap with sample clustering:
plot_contribution_heatmap(fit_res$contribution[select, ],
                          cluster_samples = TRUE,
                          method = "complete"
)

plot_compare_profiles(mut_mat[, 2], fit_res$reconstructed[, 1], profile_names = c("Original", "Reconstructed"), condensed = TRUE)


# calculate all pairwise cosine similarities
cos_sim_ori_rec <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
# extract cosine similarities per sample between original and reconstructed
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))


colnames(cos_sim_ori_rec) <- "cos_sim"
cos_sim_ori_rec$sample <- row.names(cos_sim_ori_rec)

ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = sample)) +
  geom_bar(stat = "identity", fill = "skyblue4") +
  coord_cartesian(ylim = c(0.8, 1)) +
  # coord_flip(ylim=c(0.8,1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  xlab("") +
  # Reverse order of the samples such that first is up
  xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
  theme_bw() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90)
  ) +
  # Add cut.off line
  geom_hline(aes(yintercept = .90)) +
  coord_flip()









# snv_data <- getData(infile = "data/combined_annotated_snvs_91019.txt", exclude = F)
snv_data <- getData(!sample %in% excluded_samples, infile = "data/annotated_snvs.txt")

# signatures.genome.cosmic.v3.may2019
# signatures.nature2013
signature_classifcation <- 'signatures.nature2013'

contributions <- sigTypes(snv_data = snv_data, min_contribution = 0.1)

save(contributions, file = "data/sig_contributions_cosmic_nature_no_cutoff_7420")

load("data/sig_contributions_cosmic_nature_no_cutoff_7420")
# load("data/sig_contributions_cosmin_old_no_cutoff")

library(pheatmap)
library(tidyverse)

contributions <- as.data.frame(contributions)

attach_info <- "data/samples_names_conversion.txt"

name_conversion <- read.delim(attach_info, header = F)
colnames(name_conversion) <- c("sample", "sample_short", "sample_paper", "sex", "assay")

ann_contributions <- plyr::join(contributions, name_conversion, "sample", type = "left")

ann_contributions <- ann_contributions %>%
  dplyr::mutate(sex = as.factor(ifelse(assay == "whole-gut", "whole-gut", as.character(sex))))

males <- levels(as.factor(ann_contributions$sample[ann_contributions$sex == "male"]))
females <- levels(as.factor(ann_contributions$sample[ann_contributions$sex == "female"]))
whole_gut <- levels(as.factor(ann_contributions$sample[ann_contributions$sex == "whole-gut"]))

group <- data.frame(group = c(rep("Male", length(males)), rep("Female", length(females)), rep("whole-gut", length(whole_gut))))

rownames(group) <- c(
  males,
  females,
  whole_gut
)

cont <- ann_contributions %>%
  dplyr::select(-total, -sample_short,-sample_paper, -sex, -assay) %>%
  as.data.frame()

wide <- cont %>% tidyr::spread(sample, score)

new <- wide %>%
  dplyr::mutate(sum = rowSums(.[2:ncol(.)])) %>%
  # dplyr::filter(sum >= 0.1) %>%
  droplevels()


new[is.na(new)] <- 0

sigs <- new$signature
new$signature <- NULL
new$sum <- NULL

contribution_matrix <- as.matrix(new)

rownames(contribution_matrix) <- sigs

write.table(contribution_matrix, file = "data/mut_sigs_contributions_nature_2013.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# colours <- c("#EBF4F5", "#7BB3D4", "#65ABC2", "#4594B3", "#366C9E")

# colours <- c("#FFFFFF", "#ADD8E6", "#BFEFFF", "#B2DFEE", "#9AC0CD");

colours <- c("#f1f9fe", "#d4ecfb", "#9bd3f6", "#7ec6f3", "#53b3ef")

colours <- rev(heat.colors(5, 0.5))

pheatmap::pheatmap(t(contribution_matrix),
                   color = colours,
                   cluster_row = T, cluster_cols = T,
                   treeheight_row = 0, treeheight_col = 0,
                   main = "Mutational signature contribution per sample",
                   angle_col = 90, annotation_row = group
)

library(d3heatmap)
d3heatmap(t(contribution_matrix), scale = "row")
