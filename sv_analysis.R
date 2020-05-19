# library(svBreaks)
library(tidyverse)
library(stringr)

# all_samples <- paste0(desktop, '/final_analysis/filtered/summary/merged/all_samples_merged.txt'
# infile <- paste0(desktop, '/final_analysis/filtered/summary/merged/all_bps_merged.txt'

rootDir <- ifelse(dir.exists('/Users/Nick_curie/'), '/Users/Nick_curie/', '/Users/Nick/iCloud/')
# devtools::install(paste0(rootDir, 'Desktop/script_test/svBreaks'))
devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/svBreaks'))
devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/mutationProfiles'))

bedDir = paste0(rootDir, 'Desktop/misc_bed')

all_samples <- paste0(rootDir, 'Desktop/parserTest/all_combined_23719/all_samples_merged.txt')
infile <- paste0(rootDir, 'Desktop/parserTest/all_combined_23719/all_bps_merged.txt')

attach_info <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt')
all_sample_names <- read.delim(attach_info, header = F)
colnames(all_sample_names) <- c("sample", "sample_short", "sample_paper", "sex", "assay")

excluded_samples <- c("B241R41-2",  "A373R7", "A512R17")
excluded_samples <- c(excluded_samples, "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R09", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
excluded_samples <- c(excluded_samples, "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
excluded_samples <- c(excluded_samples, "D106R1", "D106R3", "D106R5", "D106R7", "D106R9", "D106R11", "D106R13", "D106R15", "D106R17", "D106R19", "D106R21", "D106R23", "D106R25", "D106R27", "D106R29", "D106R31", "D106R33"  )
excluded_samples <- c(excluded_samples, "D197R09", "D197R11", "D197R13", "D197R15")
excluded_samples <- c(excluded_samples, "D265R01", "D265R03", "D265R05", "D265R07", "D265R09", "D265R11", "D265R13")

all_sample_names <- all_sample_names %>% dplyr::filter(!sample %in% excluded_samples)

red <- "#FC4E07"
# blue <- "#00AFBB"
blue <- "#259FBF"
yellow <- "#E7B800"
grey <- "#333333"

chromosomes <- c("2L", "2R", "3L", "3R", "X", "Y", "4")

all_hits <- svBreaks::getData(!sample %in% excluded_samples, infile=infile)

n_hits <- svBreaks::geneHit(!sample %in% excluded_samples, plot = F, all_samples = all_samples)

n_hits <- n_hits %>% 
  dplyr::mutate(key = paste0(sample, "_", event)) %>% 
  dplyr::select(sample, event, key)

non_notch <- all_hits %>% 
  dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::filter(!key2 %in% n_hits$key) %>% 
  droplevels() %>% 
  dplyr::select(-key2) %>% 
  ungroup()

non_notch %>% 
  dplyr::filter(bp_no=='bp1') %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::summarise(median(n()))

notch <- all_hits %>% 
  dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::filter(key2 %in% n_hits$key) %>% 
  droplevels() %>% 
  dplyr::select(-key2) %>% 
  ungroup()



snvs <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_snvs.txt')
indels <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_indels.txt')

male_snvs <- mutationProfiles::getData(infile = snvs, type='snv', sex=='male', !sample %in% excluded_samples,
                                       assay %in% c("neoplasia", "delta-neoplasia"), chrom %in% chromosomes,
                                       attach_info = attach_info)

male_indels <- mutationProfiles::getData(infile = indels, type='indel', sex=='male', !sample %in% excluded_samples,
                                         assay %in% c("neoplasia", "delta-neoplasia"), chrom %in% chromosomes,
                                         attach_info = attach_info )

combined_muts <- plyr::join(male_snvs, male_indels, type='full')



########################
######   FIG. 2   ######
########################

 

all_samples_df <- read.delim(all_samples, header = T) %>% 
  dplyr::filter(!sample %in% excluded_samples)

all_samples_df <- plyr::join(all_samples_df, all_sample_names, "sample", type = 'left') %>% 
  dplyr::rename(sample_old = sample,
                sample = sample_paper) %>% 
  dplyr::select(sample, everything())


# Fig 2a: SVs affecting Notch
svBreaks::geneHit(!sample %in% excluded_samples, all_samples_df = all_samples_df, show_sample = T, drivers = c("N", "kuz"), plot=T)

# Fig 2b: Svs in Notch region
svBreaks::notchHits(!sample %in% excluded_samples, all_samples_df = all_samples_df, show_samples = T, bp_density = T)

# Fig 2b2: Breakpoint density over Notch region
# svBreaks::notchHits(!sample %in% excluded_samples, all_samples = all_samples, show_samples = F, bp_density = T)
 

# S.fig 3
svBreaks::notchHits(!sample %in% excluded_samples, all_samples_df = all_samples_df, from=3.12, to=3.173, show_samples = T, bp_density = T, adjust = 0.1, ticks = 2.5)
# svBreaks::notchHits(!sample %in% excluded_samples, all_samples = all_samples, from=3.12, to=3.173, show_samples = F, bp_density = T, ticks = 5)


## Write breakpoints +/- 5kb for Notch SVs
n_hits <- svBreaks::geneHit(!sample %in% excluded_samples, plot = F, all_samples = all_samples)

# expand_by <- 5000
expand_by = 500
today <- format(Sys.time(), "%d_%m_%y")
N_breakpoints <- n_hits %>% 
  dplyr::distinct(sample, event, type, .keep_all=TRUE) %>% 
  dplyr::mutate(chrom = chromosome1) %>% 
  dplyr::select(chrom, start, end) %>% 
  tidyr::gather("bp", "pos", c("start", "end")) %>% 
  dplyr::mutate(rstart = pos - expand_by,
                rend = pos + expand_by) %>% 
  dplyr::select(-bp, -pos) %>% 
  droplevels() %>% 
  svBreaks::writeBed(df=., outDir = paste0(rootDir, 'Desktop/misc_bed/breakpoints'), name = paste0('notch_breakpoints_', expand_by, '_', today, '.bed'))

# Run `bash extract_breakpoints.sh -i ../notch_breakpoints_5k_10.3.20.bed -n 100`

# Calculate median NBR length:

read.delim(paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_5k.bed'), header = F) %>% 
  dplyr::mutate(length = V3 - V2) %>% 
  dplyr::summarise(median_length = median(length, na.rm = TRUE))

mappable_regions = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_mappable.bed'
chrom_lengths = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/chrom.sizes.txt'
overlaps <- bpRegioneR(regionA = paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_500_mappable.bed'),
           regionB = paste0(rootDir, 'Desktop/misc_bed/motifs/motif_1.mappable.bed'),
           mappable_regions = mappable_regions, n=10000, from=2700000, to=3400000, chrom = 'X', plot=F)

plot_bpRegioner(overlaps, bins = 10)

# G4s
bpRegioneR(regionA = paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_500_mappable.bed'),
              regionB = paste0(rootDir, 'Desktop/misc_bed/features/nonBform/g4s.mappable.bed'),
              mappable_regions = mappable_regions, n=100, from=2700000, to=3400000, chrom = 'X', plot=T)
# SIRs
bpRegioneR(regionA = paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_500_mappable.bed'),
                         regionB = paste0(rootDir, 'Desktop/misc_bed/features/nonBform/SIRs.merged.mappable.bed'),
                         mappable_regions = mappable_regions, n=100, from=2700000, to=3400000, chrom = 'X', plot=T)
# Cruciform
bpRegioneR(regionA = paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_500_mappable.bed'),
           regionB = paste0(rootDir, 'Desktop/misc_bed/features/nonBform/cruciform.merged.mappable.bed'),
           mappable_regions = mappable_regions, n=100, from=2700000, to=3400000, chrom = 'X', plot=T)

# Min G4s
bpRegioneR(regionA = paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_500_mappable.bed'),
           regionB = paste0(rootDir, 'Desktop/misc_bed/features/nonBform/min_g4s.bed'),
           mappable_regions = mappable_regions, n=100, from=2700000, to=3400000, chrom = 'X', plot=T)

# Fig 2d,e: Two characterised events in Notch 

# Fig 2f: TEs at Notch breakpoints
te_in <- paste0(rootDir, 'Documents/Curie/Documents/SV_paper/sv_characterisation/TE_involvement.txt')
notch_tes <- read.delim(te_in)
notch_tes <- notch_tes %>% 
  dplyr::filter(!sample %in% excluded_samples) %>% 
  dplyr::select(sample, bp1, bp2, type) %>% 
  dplyr::mutate(type2 = type, 
                bp1 = as.logical(ifelse(stringr::str_detect(bp1, ''), 1, 0)),
                bp2 = as.logical(ifelse(stringr::str_detect(bp2, ''), 1, 0))) %>% 
  dplyr::mutate(has_te = ifelse(bp1 | bp2, TRUE, FALSE)) %>% 
  dplyr::group_by(type2) %>%
  dplyr::count(type2, has_te) %>% 
  dplyr::mutate(total = sum(n), perc = (n/total)*100) %>% 
  dplyr::filter(has_te) %>% 
  as.data.frame()

# plot_tes_at_breakpoints <- function(x){
#   colours = svBreaks::sv_colours()
#   
#   p <- ggplot(x)
#   p <- p + geom_bar(aes(fct_reorder(type2, -perc), perc, colour = type2, fill = type2), alpha=0.7, stat='identity')
#   p <- p + scale_fill_manual("SV type\n", values = colours)
#   p <- p + scale_colour_manual(values = colours)
#   p <- p + scale_y_continuous(limits=c(0,100))
#   p <- p + guides(colour = FALSE)
#   p <- p + cleanTheme() +
#     theme(
#       panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
#       axis.title.x = element_blank()
#     )
#   p
# }
# 
# plot_tes_at_breakpoints(notch_tes)

# notch <- convert_names(notch)
# 
# te_tally <- notch %>% 
#   dplyr::group_by(summary) %>%
#   dplyr::summarise(count = n())
# 
# te_tally$class <- 'TE'
# te_tally$summary <- as.character(te_tally$summary)
# te_tally[is.na(te_tally)]<-"No-TE"
# 
# te_tally$summary <- as.factor(as.character(te_tally$summary))
# 
# order <- levels(fct_reorder(te_tally$summary, -te_tally$count))
# 
# ggbarplot(te_tally, "summary", "count", order = order,
#           fill = "summary", color = "summary", palette = "jco", alpha = 0.6,
#           title="TE counts at Notch breakpoints", ylab="Count", xlab=F, ggtheme = theme_minimal(), x.text.angle = 90)


########################
######   FIG. 4   ######
########################

# Fig 4a - Rainfall plot of breakpoints
all_hits <- svBreaks::getData(!sample %in% excluded_samples, infile=infile)

svBreaks::bpRainfall(bp_data = all_hits, chroms = chromosomes)

n_hits <- svBreaks::geneHit(!sample %in% excluded_samples, plot = F, all_samples = all_samples)

n_hits <- n_hits %>% 
  dplyr::mutate(key = paste0(sample, "_", event)) %>% 
  dplyr::select(sample, event, key)

non_notch <- all_hits %>% 
  dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::filter(!key2 %in% n_hits$key) %>% 
  droplevels() %>% 
  dplyr::select(-key2) %>% 
  ungroup()

non_notch %>% 
  dplyr::filter(bp_no=='bp1') %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::summarise(median(n()))
  
notch <- all_hits %>% 
  dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::filter(key2 %in% n_hits$key) %>% 
  droplevels() %>% 
  dplyr::select(-key2) %>% 
  ungroup()


# Fig. 4b Svs type percentages 
svs_condensed_types <- non_notch %>%
  dplyr::mutate(type2 = factor(type2)) %>% 
  dplyr::group_by(type2) %>%
  dplyr::distinct(chrom, event, .keep_all=TRUE) %>%
  dplyr::tally() %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1)) %>% 
  as.data.frame()

order <- levels(fct_reorder(svs_condensed_types$type2, -svs_condensed_types$perc))

svs_condensed_types %>% 
  ggplot2::ggplot(.) +
  ggplot2::geom_bar(aes(fct_reorder(type2, -perc), perc), colour = grey, fill = grey, alpha=0.6, stat='identity') + 
  ggplot2::scale_y_continuous("Percentage", limits=c(0,100)) + 
  ggplot2::scale_colour_manual(values = colours) +
  cleanTheme() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y =element_text(size=15)
  )

noR1 <- non_notch %>% dplyr::filter(sample != "A373R1")
# noR1 <- transform_types(noR1)
noR1_samples <- c("A373R1", excluded_samples)
# Get SV tpes per sample

sv_data <- svTypes(!sample %in% noR1_samples, bp_data = non_notch, plot = F)

sv_df <- sv_data[[1]]

sv_df2 <- swapSampleNames(df=sv_df)
levels(sv_df2$sample)

sv_df2$sample <- fct_reorder(sv_df2$sample, -sv_df2$sv_count)

# Sup fig 6 - af of mutations by type
gghistogram(non_notch, x = "af", color = "type2", fill= "type2", alpha = 0.5,
            add = "median", bins=50, title="VAF distribution of SV types",
            ylab="Number of SVs", xlab = "VAF",
            ggtheme = theme_minimal())


# Fig 4c - Protein coding mutations (non-Notch)
snv_indel_df <- read.csv(paste0(rootDir, 'Desktop/final_analysis/annotated_mutations_filt.csv'))

# tally_impacts <- ann_data %>%
#   dplyr::filter(!is.na(impact)) %>%
#   dplyr::group_by(sample, impact) %>%
#   tally() %>%
#   dplyr::group_by(sample) %>%
#   dplyr::mutate(total = sum(n))

# order <- levels(fct_reorder(tally_impacts$sample, tally_impacts$total))

tally_impacts <- snv_indel_df %>%
  dplyr::filter(!is.na(impact),
                sex == 'male',
                !sample_old %in% noR1_samples) %>%
  dplyr::group_by(sample, impact) %>%
  dplyr::mutate(type_count = n()) %>% 
  dplyr::group_by(sample) %>%
  dplyr::add_tally()  %>% 
  dplyr::distinct(impact, .keep_all=T) %>%
  as.data.frame()

tally_impacts$sample <- factor(tally_impacts$sample, levels = levels(sv_df2$sample))

# tally_impacts$sample <- fct_reorder(tally_impacts$sample, -sv_df$sv_count)

# Median protein-coding muts
median(tally_impacts$n)
table(tally_impacts$impact)

# Fig 4e - coding consequences
tally_impacts %>% 
  dplyr::group_by(impact) %>% 
  dplyr::summarise(total = sum(type_count)) %>% 
  ggplot2::ggplot(.) + 
  ggplot2::geom_bar(aes(fct_reorder(impact, -total), total), colour = grey, fill = grey, alpha=0.6, stat = 'identity') +
  cleanTheme() +
  theme(
    axis.title.y = element_text(size=18),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
    axis.title.x = element_blank(),
    legend.position = 'none'
  )


colours <- sv_colours()

p1 <- ggplot(sv_df2)
p1 <- p1 + geom_bar(aes(fct_reorder(sample, -sv_count), type_count, colour = type2, fill = type2), alpha=0.7, stat = 'identity')
p1 <- p1 + scale_fill_manual("SV type\n", values = colours)
p1 <- p1 + scale_colour_manual(values = colours)
p1 <- p1 + guides(color = FALSE)
p1 <- p1 + cleanTheme() +
  theme(
    axis.title.y = element_blank(),
    legend.position = 'none',
    # panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
    axis.title.x=element_blank()
  )

tally_impacts$sample <- factor(tally_impacts$sample, levels = levels(sv_df2$sample))
colours <- snv_colours()

p2 <- ggplot(tally_impacts)
p2 <- p2 + geom_bar(aes(sample, type_count, colour = impact, fill = impact), alpha=0.7, stat = 'identity')
p2 <- p2 + scale_fill_manual("SV type\n", values = colours)
p2 <- p2 + scale_colour_manual(values = colours)
# p2 <- p2 + guides(color = FALSE, impact=FALSE, type_count=FALSE)
p2 <- p2 + cleanTheme() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
    legend.position = 'none'
    )


# cowplot::plot_grid(p1 + coord_flip(),
#                    p2 + coord_flip(),
#                    align="h", rel_widths = c(1, 1))


cowplot::plot_grid(p2, p1, align="v", axis = 'b', nrow = 2)

# Fig 3e - non-Notch TEs at breakpoints
# tes_at_breakpoints <- read.delim(all_samples) %>% 
#   dplyr::mutate(type = factor(type)) %>% 
#   dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
#   dplyr::mutate(type2 = as.character(ifelse(stringr::str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>% 
#   dplyr::mutate(type2 = factor(type2)) %>% 
#   dplyr::filter(!sample %in% excluded_samples) %>% 
#   dplyr::mutate(te_involvement_1 = grepl("bp1_", notes),
#                 te_involvement_2 = grepl("bp2_", notes)) %>%
#   dplyr::group_by(sample, event) %>% 
#   dplyr::filter(!key2 %in% n_hits$key) %>% 
#   dplyr::mutate(count_1 = as.numeric(ifelse(te_involvement_1, as.numeric(str_replace(notes, ".*=(\\d+).*", "\\1")), 0))) %>% 
#   dplyr::mutate(count_2 = as.numeric(ifelse(te_involvement_2, as.numeric(str_replace(notes, ".*=(\\d+).*", "\\1")), 0))) %>% 
#   dplyr::mutate(count_1 = ifelse(is.na(count_1), 0, count_1),
#                 count_2 = ifelse(is.na(count_2), 0, count_2)) %>% 
#   dplyr::mutate(count = count_1 + count_2) %>% 
#   dplyr::select(sample, event, type2, notes, te_involvement_1, te_involvement_2, count_1, count_2, count) %>% 
#   dplyr::distinct(event, sample, type2, .keep_all=T) %>% 
#   dplyr::mutate(has_te = ifelse(te_involvement_1 | te_involvement_2 && count_1 > 1 | count_2 > 1, TRUE, FALSE)) %>%  
#   dplyr::group_by(type2) %>%
#   dplyr::count(type2, has_te) %>% 
#   dplyr::mutate(total = sum(n), perc = (n/total)*100) %>% 
#   as.data.frame()
# 
# colours = svBreaks::sv_colours()
# 
# tes_at_breakpoints <- tes_at_breakpoints %>% dplyr::filter(has_te)
# 
# plot_tes_at_breakpoints(tes_at_breakpoints)

# Tally types of sv in Notch
# n_hits %>% View()


########################
######   FIG. 5   ######
########################
bedDir = "/Users/Nick_curie/Desktop/misc_bed"

source(paste0(rootDir, 'Desktop/Analysis_pipelines/figs_for_paper_1/R/enrichment_funs.R'))

# mutationProfiles::rainfall(snv_data = combined_muts, chroms = chromosomes)


####################
## SNVS + INDELS  ##
####################

# devtools::install(paste0(rootDir, 'Desktop/script_test/mutationProfiles'))
devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/mutationProfiles'))

snvs <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_snvs.txt')
indels <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_indels.txt')

male_snvs <- mutationProfiles::getData(infile = snvs, type='snv', sex=='male', !sample %in% excluded_samples,
                                       assay %in% c("neoplasia", "delta-neoplasia"), chrom %in% chromosomes,
                                       attach_info = attach_info)

male_indels <- mutationProfiles::getData(infile = indels, type='indel', sex=='male', !sample %in% excluded_samples,
                                         assay %in% c("neoplasia", "delta-neoplasia"), chrom %in% chromosomes,
                                         attach_info = attach_info )

# male_snvs <- male_snvs %>%
#   dplyr::rename(sample_old = sample,
#                 sample = sample_paper) %>%
#   dplyr::select(sample, everything())
# 
# male_indels <- male_indels %>%
#   dplyr::rename(sample_old = sample,
#                 sample = sample_paper) %>%
#   dplyr::select(sample, everything())

combined_muts <- plyr::join(male_snvs, male_indels, type='full')


## Old - muts in genes by expression
# r <- hits_in_genes(sample != "A373R1", df = non_notch, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), out_file='bps_in_genes.txt')
# bps_in_genes <- r[[2]]
# bps_in_genes$mut <- 'sv'
# 
# 
# r <- hits_in_genes(sample != "A373R1", df = male_snvs, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), out_file='snvs_in_genes.txt')
# snvs_in_genes <- r[[2]]
# snvs_in_genes$mut <- 'snv'
# 
# r <- hits_in_genes(sample != "A373R1", df = male_indels, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), out_file='indels_in_genes.txt')
# indels_in_genes <- r[[2]]
# indels_in_genes$mut <- 'indel'
# 
# comb <- rbind(bps_in_genes, snvs_in_genes, indels_in_genes)
# 
# comb <- comb %>% 
#   tidyr::complete(group, feature, mut) %>% 
#   dplyr::mutate(Log2FC = ifelse(is.na(Log2FC), 0, Log2FC),
#                 label = ifelse(sig %in% c(NA,'-'), '', sig)) %>% 
#   dplyr::mutate(observed = replace_na(observed, 0)) %>% 
#   as.data.frame()
# 
# 
# comb$mut = factor(comb$mut, levels=c("sv","snv", "indel"), labels=c("SV","SNV", "INDEL")) 
# # comb$group = factor(comb$group, levels=c("All expressed genes", "Top 20% FPKM", "Non-expressed genes"), labels=c("All expressed genes", "Top 20% FPKM", "Non-expressed genes")) 
# comb$group = factor(comb$group, levels=c("All expressed genes", "Non-expressed genes"), labels=c("ISC-expressed genes", "Non-expressed genes")) 
# 
# # Figure 5 a - mutations in gene features
# gf_exp_grouped <- ggbarplot(comb, x = "feature", y = "Log2FC",
#                             fill = "mut", color = "mut", group = "group", palette =c(blue, yellow, red, grey), alpha = 0.6, 
#                             label = 'label', lab.vjust = -.1,
#                             orientation = "horizontal", x.text.angle = 90, position = position_dodge(0.8), 
#                             xlab = FALSE, ylab = "Log2(FC)", ggtheme = theme_pubr())
# 
# gf_exp_grouped <- ggpar(gf_exp_grouped, ylim = c(-2.2,2.2))
# facet(gf_exp_grouped, facet.by = 'group', nrow = 1)


# Gene features
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/umr_gene_features/mappable/paper'))

comb <- bind_rows(d)

comb <- comb %>% 
  tidyr::complete(group, feature) %>% 
  dplyr::mutate(Log2FC = ifelse(is.na(Log2FC), 0, Log2FC),
                label = ifelse(sig %in% c(NA,'-'), '', sig)) %>% 
  dplyr::mutate(observed = replace_na(observed, 0)) %>% 
  as.data.frame()

comb$group = factor(comb$group, levels=c("sv","snv", "indel"), labels=c("SV","SNV", "INDEL")) 
comb$feature = factor(comb$feature, levels=c("gene", "5UTR", "3UTR", "CDS", "Intronic", "Intergenic"), labels=c("Geneic", "5'UTR", "3'UTR", "CDS", "Intron", "Intergenic")) 

# Figure 5 a - mutations in gene features
gf_exp_grouped <- ggbarplot(comb, x = "feature", y = "Log2FC",
                            fill = "group", color = "group", group = "group", palette =c(blue, yellow, red, grey), alpha = 0.6, 
                            label = 'label', lab.vjust = 2,
                            orientation = "horizontal", x.text.angle = 90, position = position_dodge(0.8), 
                            xlab = FALSE, ylab = "Log2(FC)", ggtheme = theme_pubr())

ggpar(gf_exp_grouped, ylim = c(-2.2,2.2))


facet(gf_exp_grouped, facet.by = 'group', nrow = 1)




# Modencode enriched
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/enriched/merged/'))
plot_all(d)

# Manon DamID
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/Manon_2020/deseq/'))
plot_all(d)

# repeat seqs
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/repeats/mappable/'))
plot_all(d)

# TF binding
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/jaspar/'))
plot_all(d)

# Generic features
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/features/'))
plot_all(d)


# reldist

# Bedtools reldist
# library(bedr)

bps_as_snvs <- non_notch %>% 
  dplyr::rename(pos = bp) %>% 
  as.data.frame()
  
all_mutations <- bind_rows(list(bps_as_snvs, combined_muts)) %>%
  dplyr::filter(chrom %in% chromosomes,sample != "A373R1") %>% 
  dplyr::select(chrom, pos) %>% 
  dplyr::distinct(chrom, pos)


dist <- run_reldist(a = all_mutations,
                    b = paste0(bedDir, '/repeats/mappable/combined/repeats_merged.bed'),
                    chroms = chromosomes, verbose = FALSE, type='snv', print = FALSE, sim=T)

plot_reldist(dist)



########################
######   FIG. 6   ######
########################

# Fig 6c & d
source(paste0(rootDir, 'Desktop/script_test/svSupport/scripts/plotFreq.R'))

samples_with_no_N_sv <- c(excluded_samples, 'A373R9', 'B241R43', 'B241R49', 'D197R05')

plot_tumour_evolution(all_samples = all_samples, 
                      !sample %in% samples_with_no_N_sv, 
                      tes = FALSE, 
                      annotate_with = paste0(rootDir, 'Desktop/gene2bed/dna_damage_merged.bed'))


tumour_evolution <- plot_tumour_evolution(all_samples = all_samples, 
                                          !sample %in% samples_with_no_N_sv, 
                                          tes = FALSE, annotate_with = paste0(rootDir, 'Desktop/gene2bed/merged_groups.bed'),
                                          plot=F)

notch_hits <- tumour_evolution[[1]]
all_muts <- tumour_evolution[[2]]

n_col <- '#56D66F'
# sv_col <- '#259FBF'
# indel_col <- '#F5D658'
snv_col <- '#F5D658'

all_muts <- all_muts %>%
  dplyr::mutate(class = ifelse(class=='SV' & gene_hit != 'Other', 'Notch', as.character(class)),
                time = 1 - cell_fraction)

all_muts$group = factor(all_muts$class, levels=c("Notch","SV","SNV", "INDEL"), labels=c("Notch","SV","SNV", "INDEL")) 

all_muts$medN <- 1-(median(notch_hits$highest_n))

ggplot(all_muts, aes(time, colour = group)) + 
  geom_line(stat='ecdf', size=1.5, alpha=.7) + 
  scale_y_continuous('ECDF') +
  scale_x_continuous('Pseudotime (1 - cell fraction)') +
  cleanTheme() +
  theme(
    axis.text = element_text(size=rel(1.1)),
    axis.title = element_text(size=rel(1.4)),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")
    ) +
  geom_vline(aes(xintercept = 1-(median(notch_hits$highest_n))),
                      linetype = "dashed", size = .7, alpha=.9) + 
  scale_fill_manual(values = c(n_col, blue, snv_col, red )) +
  scale_colour_manual(values = c(n_col, blue, snv_col, red )) + 
  guides(size = FALSE, colour = FALSE) 

# all_muts <- all_muts %>% 
#   dplyr::mutate(class = ifelse(class=='SV' & gene_hit != 'Other', 'Notch', as.character(class)),
#                 time = 1 - cell_fraction)  
# 
# all_muts$group = factor(all_muts$class, levels=c("Notch","SV","SNV", "INDEL"), labels=c("Notch","SV","SNV", "INDEL")) 
# 
# p <- gghistogram(all_muts, x = "time", color = 'group', fill= 'group', palette =c('darkgrey',  blue, yellow, red), alpha = 0.5, rug= TRUE,
#                  add = "median", bins=10, title="Tumour evolution",
#                  ylab="Number of mutations", xlab = "Pseudotime",
#                  ggtheme = theme_minimal())
# 
# facet(p, facet.by = "group", ncol = 1, scales = 'free_y')

# Pre and post Notch hits
snv_indel_df <- read.csv(paste0(rootDir, 'Desktop/final_analysis/annotated_mutations_filt.csv')) %>% 
  dplyr::rename(sample_new = sample,
                sample = sample_old) 


notch_svs <- notch_hits %>% dplyr::select(sample, highest_n)
combined <- plyr::join(all_muts, notch_svs, c("sample"), type = "left")


coding2notch <- plyr::join(snv_indel_df, notch_svs, c("sample"), type = "left") %>% 
  dplyr::filter(!sample %in% samples_with_no_N_sv) %>% 
  dplyr::mutate(clonality = ifelse(cell_fraction > highest_n, 'pre', 'post'),
                time = 1 - cell_fraction)


ggplot(coding2notch) +
  geom_density(aes(time, colour=impact, fill=impact), alpha=0.2) +
  geom_rug(aes(time, colour=impact) )


combined %>% 
  dplyr::filter(group != 'Notch', !sample %in% samples_with_no_N_sv) %>% 
  dplyr::group_by(sample, group) %>% 
  # dplyr::add_count(name = 'total') %>% 
  dplyr::mutate(rel_to_notch = ifelse(cell_fraction > highest_n, 'pre', 'post')) %>% 
  dplyr::group_by(rel_to_notch) %>% 
  dplyr::tally() 


# Annotate notch inactivating mut on all hits
rel2notch <- plyr::join(combined_muts, notch_svs, c("sample"), type = "left") %>% 
  dplyr::filter(!sample %in% samples_with_no_N_sv) %>% 
  dplyr::mutate(clonality = ifelse(cell_fraction > highest_n, 'pre', 'post') )
  # dplyr::select(-c(pos:decomposed_tri, af, caller, status, snpEff_anno, sample_short, id, sex, assay, type) ) %>% 
  # dplyr::group_split(clonality)



rel2notch[[2]] %>% 
  group_by(gene) %>% 
  tally() %>% 
  arrange(-n)


d1 <- bpRegionEnrichment(sample != "A373R1", clonality =='pre', bp_data = rel2notch, dataType='mutationProfiles', bedDir = paste0(bedDir, '/enriched/merged/'), plot=F) %>% 
  dplyr::mutate(group = 'pre-notch')

d2 <- bpRegionEnrichment(sample != "A373R1", clonality == 'post', bp_data = rel2notch, dataType='mutationProfiles', bedDir = paste0(bedDir, '/enriched/merged/'), plot=F) %>% 
  dplyr::mutate(group = 'post-notch')

plot_all(l = list(d1, d2)) 


mappable_regions = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_mappable.bed'
chrom_lengths = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/chrom.sizes.txt'
bpRegioneR(regionA = paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_5k_mappable.bed'),
           regionB = paste0(rootDir, 'Desktop/misc_bed/motifs/motif_1.mappable.bed'),
           mappable_regions = mappable_regions, n=5000, from=2700000, to=3400000)


# Get svs in Notch
notch <- all_hits %>% 
  dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::filter(key2 %in% n_hits$key) %>% 
  dplyr::select(-key2) %>% 
  droplevels() %>% 
  ungroup()


write_breakpoints <- function(..., x, file_name='breakpoints.bed', outDir=paste0(rootDir, 'Desktop/misc_bed/breakpoints')){
  x <- x %>% 
    dplyr::filter(...) %>% 
    dplyr::rename(start = bp) %>% 
    dplyr::mutate(end = start+1) %>%
    dplyr::distinct(event, sample, chrom, start, .keep_all=T) %>% 
    dplyr::mutate(info = paste(sample, event, type, bp_no, sep="_")) %>% 
    dplyr::select(chrom, start, end, info) %>% 
    svBreaks::writeBed(df= ., name = file_name, outDir=outDir)
}

write_breakpoints(x=notch, file_name = 'notch_breakpoints_30.4.20.bed')





# Fix for non-existing dir
# R39 <- notch %>% dplyr::filter(sample == 'B241R39')
# svBreaks::writeBed(df=R39, svBreaks = T, name = 'R39_notch_svs.bed', outDir = paste0(desktop, '/misc_bed/breakpoints/')
# svBreaks::writeBed(df = all_dels,  svBreaks = T, name = 'all_prec_dels.bed', outDir = paste0(desktop, '/misc_bed/breakpoints')

svBreaks::writeBed(df=notch, svBreaks = T, name = 'Notch_svs.bed', outDir = paste0(rootDir, 'Desktop/misc_bed/breakpoints/'))
svBreaks::writeBed(df=non_notch, svBreaks = T, name = 'non-Notch_svs.bed', outDir = paste0(rootDir, 'Desktop/misc_bed/breakpoints/'))

# Fig 2a: sv count per sample
# svBreaks::bpStats(bp_data = non_notch, sample != "A37 3R1")
svBreaks::svTypes(bp_data = non_notch)

noR1 <- all_hits %>% dplyr::filter(!sample %in% c("A373R1", "D106R25"))

svBreaks::bpRainfall(bp_data = all_hits, chroms = c("2L", "2R", "3L", "3R", "X", "Y", "4"))

non_notch_slimmed <- non_notch %>% 
  dplyr::mutate(type2 = as.character(ifelse(stringr::str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type))),
                allele_frequency = af) %>% 
  dplyr::mutate(type2 = factor(type2))

non_notch_noR1 <- non_notch_slimmed %>% dplyr::filter(sample != 'A373R1') %>% droplevels()

nn <- svBreaks::tally_hits(bp_data = non_notch_noR1, plot=F, freq = T, all_samples = all_samples)

n$group = "Notch"
nn$group = 'Other'

n$class <- factor(n$class)
missingclasses = list()
df <- merge(n,nn, all = T)

sv_classes <- levels(non_notch_slimmed$type2)

for(i in 1:length(sv_classes)){
  if(!(sv_classes[i] %in% levels(n$class))){
    missingclasses[i] <- sv_classes[i]
  }
}

missingclasses <- plyr::compact(missingclasses)
dat <- data.frame(class = c(levels(n$class), paste(missingclasses)),
                  group = "Notch",
                  frequency = c(n$frequency, rep(0, length(missingclasses))),
                  count = c(n$count, rep(0, length(missingclasses)))) 

df <- merge(dat,nn, all = T)

df <- df %>% 
  dplyr::mutate(class = as.character(ifelse(class == "BND", "INV", as.character(class))))

p <- ggplot(df) + geom_bar(aes(fct_reorder(class, -frequency), frequency, fill = group, group=group), stat="identity", width=.9, alpha=0.7, position = "dodge")
# geom_bar(stat="identity", width=.5, position = "dodge")
p <- p + scale_y_continuous('Percent')
p <- p + cleanTheme() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20)
  )
p <- p + scale_fill_manual(values = c("#0073C299", "#E7B800"))
p <- p + scale_colour_manual(values = c("#0073C299", "#E7B800"))
# p <- p + facet_wrap(~group, ncol = 1)
p



svBreaks::bpRainfall(bp_data = all_hits)


bpFeatureEnrichmentPlot(confidence == 'precise')




### mechanisms notch vs non-notch
nn_mechs <- svBreaks::mechansimSize(bp_data = non_notch, plot = F)
n_mechs <- mechansimSize(bp_data = notch, plot = F)

n_mechs$group <- 'Notch'
nn_mechs$group <- 'Other'
n_nn_mechs <- merge(n_mechs, nn_mechs, all = T)

mechs <- n_nn_mechs %>% 
  dplyr::group_by(group) %>%
  dplyr::mutate(mech_count = n()) %>%
  dplyr::mutate(freq = mech_count / sum(mech_count)) %>% 
  dplyr::ungroup()

# proportion of mechanism in Notch vs non-notch
p <- ggplot(mechs)
p <- p + geom_bar(aes(x=fct_reorder(mechanism, -freq), y = ..prop..*100, fill = group, colour=group, group = group),  width=.9, alpha=0.7, position = "dodge")
p <- p + scale_y_continuous('Percent')
# p <- p + facet_wrap(~group, ncol = 1)
p <- p + slideTheme() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20)
  )
p <- p + scale_fill_manual(values = c("#0073C299", "#E7B800"))
p <- p + scale_colour_manual(values = c("#0073C299", "#E7B800"))


p

ggbarplot(mechs, x = "mechanism", y = "freq", sort.val = "desc",
                    fill = "group", color = "group", palette = "jco", alpha = 0.8, 
                    x.text.angle = 90)


# ggbarplot(mechs, "mechanism", "freq",
#           fill = "group", color = "group", palette = c("#00AFBB", "#E7B800"),
#           label = TRUE,
#           position = position_dodge(0.9))




p2 <- ggplot(mechs) + geom_bar(aes(y = freq, x = group, fill = mechanism), stat="identity")
p2 <- p2 + scale_y_continuous('Frequency')
p2 <- p2 + slideTheme() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20)
  )

p2

# microhomology length
svBreaks::micromologyPlot(bp_data = all_hits, stringr::str_detect(type, 'COMPLEX'))

notch_bps <- paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_svs.bed')
non_notch_bps <- paste0(rootDir, 'Desktop/misc_bed/breakpoints/non-Notch_svs.bed')

### distance to tss
featureDir = paste0(desktop, '/misc_bed/tss_locs/'
featureDir = '/Users/Nick_curie/Desktop/misc_bed/tss_locs/Meers_2108'
featureDir = paste0(desktop, '/misc_bed/repStress/'

featureDir = paste0(desktop, '/misc_bed/seqs_surrounding_snvs/'

featureDir = "~/Desktop/misc_bed/features/nonBform/"
featureDir = "~/Desktop/misc_bed/damID/"
featureDir = '~/GitHub/BardinLab/meme/out/notch_CFS/fimo_out/'

nhej <- all_hits %>%  dplyr::filter(stringr::str_detect(mechanism, 'NHEJ')) %>% droplevels()
alt <- all_hits %>%  dplyr::filter(stringr::str_detect(mechanism, 'Alt-EJ')) %>% droplevels()
fostes <- all_hits %>%  dplyr::filter(stringr::str_detect(mechanism, 'FoSTeS')) %>% droplevels()

all_dels <- all_hits %>% dplyr::filter(stringr::str_detect(type, 'DEL'),
                                       confidence=='precise') %>% droplevels()

notch_dels <- notch %>% dplyr::filter(stringr::str_detect(type, 'DEL'),
                                         confidence=='precise') %>% droplevels()
d <- svBreaks::distOverlay2(df = non_notch,
                            featureDir = featureDir,
                            position = 'edge',
                            lim = 2,
                            n = 100,
                            plot=T,  threshold = 0.01, write=F, out_dir = paste0(desktop, '/misc_bed/bpSimulations')


### notch to motfis ####

featureDir = '~/GitHub/BardinLab/meme/out/notch_CFS/fimo_out/'
d <- svBreaks::distOverlay2(df = notch,
                            featureDir = featureDir,
                            position = 'edge',
                            lim = 1,
                            n = 500,
                            plot=F,  threshold = 0.01, write=F, out_dir = paste0(desktop, '/misc_bed/bpSimulations')


### RPA70 ###

## Do we see an enrichment of motif1 in RPA-bound regions?

bpRegioneR(regionA = paste0(desktop, '/misc_bed/damID/mappable/RPA70_damid.mappable.bed',
           regionB = paste0(desktop, '/misc_bed/motifs/motif_1.mappable.bed',
           n=500)

d <- svBreaks::bpRegionEnrichment(
                                  bedDir = featureDir)

all_no_TRA <- all_hits %>% dplyr::filter(!stringr::str_detect(type, 'TRA'), sample!='A373R1') %>% droplevels()

svBreaks::writeBed(df=all_no_TRA, svBreaks = T, name = 'all_breakpoints.bed', outDir = paste0(desktop, '/misc_bed/')

all_dels <- all_hits %>% dplyr::filter(stringr::str_detect(type, 'DEL'),
                                         confidence=='precise') %>% droplevels()

svBreaks::writeBed(df=all_dels, svBreaks = T, name = 'all_precise_dels.bed', outDir = paste0(desktop, '/misc_bed/')

featureDir = paste0(desktop, '/misc_bed/damID/'
featureDir = 'inst/extdata/features/RepOris/'

d <- svBreaks::bpRegionEnrichment(!sample %in% excluded_samples, slop = 5,
                                  bedDir = featureDir)


d <- svBreaks::distOverlay2(
                            df=non_notch_del,
                            featureDir = featureDir,
                            position = 'edge',
                            lim = 5,
                            n = 2)


prec_svs <- all_hits %>% dplyr::filter(confidence=='precise')

featureDir = paste0(desktop, '/misc_bed/enriched/'
featureDir = 'inst/extdata/features/nonBform/'

complex_svs <- all_hits %>% 
  dplyr::filter(stringr::str_detect(type, 'COMPLEX')) %>% 
  dplyr::filter(!stringr::str_detect(type, 'TRA'),
                sample != 'A373R1')

dels <- non_notch %>% 
  dplyr::filter(stringr::str_detect(type, 'DEL'),
                sample != 'A373R1')

non_notch_noR1_prec_del <- non_notch %>% 
  dplyr::filter(sample != 'A373R1',
                confidence=='precise',
                stringr::str_detect(type, 'DEL')) %>% 
  droplevels()

non_notch_del <- non_notch %>% 
  dplyr::filter(stringr::str_detect(type, 'DEL')) %>% 
  droplevels()

svBreaks::writeBed(df=non_notch_noR1_prec_del, svBreaks = T, name = 'all_dels_non_notch.bed', outDir = paste0(desktop, '/misc_bed/breakpoints/')
svBreaks::writeBed(df=non_notch, svBreaks = T, name = 'all_svs_non_notch.bed', outDir = paste0(desktop, '/misc_bed/breakpoints/')


featureDir = paste0(desktop, '/misc_bed/geneFeatures/'
featureDir = paste0(desktop, '/misc_bed/damID/'
featureDir = paste0(desktop, '/misc_bed/genes/TSS/'
featureDir = paste0(desktop, '/misc_bed/repeats/'

d <- svBreaks::distOverlay2(df=all_hits,
                            featureDir = featureDir,
                            position = 'edge',
                            lim = 25,
                            n = 50,
                            threshold = 0.001,
                            plot=T, histo=F, write=F, out_dir = paste0(desktop, '/misc_bed/bpSimulations')


svBreaks::bpRegionEnrichment(bp_data = complex_svs,
                             bedDir = featureDir, slop=1000)


svBreaks::plotdistanceOverlay2(distances = d[[1]], histo=F)

# Notch exclded 
non_notch_dels <- non_notch %>%  dplyr::ungroup() %>% dplyr::filter(stringr::str_detect(type, 'DEL'),
                                               confidence=='precise') %>% dplyr::select(chrom, bp, bp2) %>% droplevels()

svBreaks::writeBed(df=non_notch_dels, svBreaks = T, name = 'all_precise_dels_non_notch.bed', outDir = paste0(desktop, '/misc_bed/')

d <- svBreaks::distOverlay(!sample %in% excluded_samples,
                           stringr::str_detect(type, 'DEL'),
                           confidence == 'precise',
                           feature_file = paste0(desktop, '/misc_bed/damID/mappable/RPA70_damid.mappable.bed',
                           position = 'edge',
                           lim = 10,
                           n = 20)

# most highly expressed genes & longest genes are enriched for complex event breakpoints
svBreaks::bpRegionEnrichment(bp_data = complex_svs, 
                            bedDir = paste0(desktop, '/misc_bed/genes/', slop=100)


cn_events <- non_notch %>% 
  dplyr::filter(confidence!='precise') %>% 
  droplevels() %>% 
  svBreaks::writeBed(df=., svBreaks = T, name = 'non_notch_cn_events.bed', outDir = paste0(desktop, '/misc_bed/breakpoints/')


## Wholegut DELs

whole_gut <- c("D050R10", "D050R12", "D050R14", 
               "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")

wg_hits <- svBreaks::getData(sample %in% whole_gut, stringr::str_detect(type, 'DEL'))
svBreaks::writeBed(df=wg_hits, svBreaks = T, name = 'whole_git_dels.bed', outDir = paste0(desktop, '/misc_bed/')




## Write breakpoints +/- 5kb for Notch SVs

n_hits <- svBreaks::geneHit(!sample %in% excluded_samples, plot = F, all_samples = all_samples)

keycol <- "condition"
valuecol <- "measurement"
gathercols <- c("control", "cond1", "cond2")

N_breakpoints <- n_hits %>% 
  dplyr::distinct(sample, event, type, .keep_all=TRUE) %>% 
  dplyr::mutate(chrom = chromosome1) %>% 
  dplyr::select(chrom, start, end) %>% 
  tidyr::gather("bp", "pos", c("start", "end")) %>% 
  dplyr::mutate(rstart = pos - 5000,
                rend = pos + 5000) %>% 
  dplyr::select(-bp, -pos) %>% 
  droplevels() %>% 
  svBreaks::writeBed(df=., outDir = paste0(desktop, '/misc_bed/breakpoints', name = 'notch_breakpoints_10k.bed')

# See additional steps to create 'notch_10kb_merged_breakpoints.bed'
bpRegioneR(regionA = paste0(desktop, '/misc_bed/breakpoints/Notch_CFS/notch_10kb_merged_breakpoints.bed',
           regionB = paste0(desktop, '/misc_bed/motifs/motif_1.mappable.bed', n=1000, from=2700000, to=3400000)


svBreaks::bpRegionEnrichment(bed_file = paste0(desktop, '/misc_bed/motifs/motif_1.mappable.bed', bedDir = paste0(desktop, '/misc_bed/damID/')

loh <- 'D106R11'
noloh <- 'D106R23'

gsub("\\R", "", loh)

bpRegionEnrichment(bp_data = notch, bedDir = paste0(desktop, '/misc_bed/repeats/', slop = 1000)

d <- svBreaks::distOverlay2(df = notch,
                           featureDir = paste0(desktop, '/misc_bed/damID/',
                           # feature_file = paste0(desktop, '/trf/out.bed',
                           position = 'edge',
                           threshold = 0.001,
                           lim = 1,
                           n = 20)
