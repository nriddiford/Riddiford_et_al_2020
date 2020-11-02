### INIT START #####

library(devtools)

rootDir <- ifelse(dir.exists('/Users/Nick_curie/'), '/Users/Nick_curie/', '/Users/Nick/iCloud/')

devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/svBreaks'))
devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/mutationProfiles'))

# library(svBreaks)
# library(mutationProfiles)

library(tidyverse)
library(stringr)
library(bedr)

# setwd('~/Desktop/script_test/SV_paper_20/')

dfsObj <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/muts.RData")
load(dfsObj)

######
## Filter snvs for those in DGPR
######


dgrp_snps_file <- '/Users/Nick_curie/Documents/Curie/Documents/SV_paper/analysis/dgrp_snps.bed'
dgrp_snps <- read.delim(dgrp_snps_file, header = T) %>% 
  dplyr::mutate(key = paste(chrom, start, sep='_')) %>% 
  dplyr::select(key)

test_overlaps <- c('2L_1938632', '2L_10499828')
  
overlaps <- male_snvs %>% 
  dplyr::filter(paste(chrom, pos, sep = '_') %in% dgrp_snps$key)

cat("There are", nrow(overlaps), "overlaps")

sim_data <- as.data.frame(bpSim(n=nrow(dgrp_snps))) %>% 
  dplyr::mutate(key = paste(chrom, pos, sep='_')) %>% 
  dplyr::select(key)

male_snvs %>% 
  dplyr::filter(paste(chrom, pos, sep = '_') %in% sim_data$key) %>% 
  nrow()

data <- data.frame(type = c('snvs', 'dgrp', 'overlaps'),
                   count = c(nrow(male_snvs), nrow(dgrp_snps), nrow(overlaps)))

ggplot(data, aes(type, count)) +
  geom_bar(aes(fct_reorder(type, -count), count), stat='identity')



##########################
## Fig 6 for WG samples ##
##########################

source('R/plotFreq.R')

wholegut_samples <- c('D050R10', 'D050R12', 'D050R14', 'D050R16', 'D050R18', 'D050R20', 'D050R22', 'D050R24')
wholegut_samples_n <- c('Fd10', 'Fd12', 'Fd14', 'Fd16', 'Fd18', 'Fd20', 'Fd22', 'Fd24')

# read.delim(all_samples) %>% 
#   dplyr::filter(sample %in% wholegut_samples) %>% 
#   write.table(., file = 'data/all_WG_samples_merged_filt.txt', quote=FALSE, sep='\t', row.names = FALSE)

# Need to annotate this ^ file with whole gut dels: 'all_WG_samples_merged_filt_annotated.txt'
wholegut_svs_af <- read.delim('data/all_WG_samples_merged_filt_annotated.txt') %>%
  dplyr::mutate(wg_del = grepl("wg_del", notes)) %>%
  dplyr::filter(!wg_del) %>% 
  dplyr::mutate(type = 'SV') %>% 
  dplyr::select(sample, type, allele_frequency)

wholegut_snvs_af <- read.delim('data/annotated_snvs_wholegut.txt') %>% 
  dplyr::mutate(type = 'SNV') %>% 
  dplyr::select(sample, type, allele_frequency)

wholegut_indels_af <- read.delim('data/annotated_indels_wholegut.txt') %>% 
  dplyr::mutate(type = 'INDEL') %>% 
  dplyr::select(sample, type, allele_frequency)

wholegut_muts_combined <- dplyr::bind_rows(wholegut_svs_af, wholegut_snvs_af, wholegut_indels_af) %>% 
  dplyr::filter(allele_frequency > 0)


##### Compare real data allele frequency dist to whole gut

# read.delim(all_samples) %>% 
#   dplyr::filter(!paste(sample, event, sep = '_') %in% excluded_events) %>% 
#   write.table(., file = 'data/all_samples_merged_filt.txt', quote=FALSE, sep='\t', row.names = FALSE)


tumour_evolution <- plot_tumour_evolution(all_samples = 'data/all_samples_merged_filt.txt', 
                                          tes = FALSE, annotate_with = 'data/dna_damage_merged.bed', plot=F)

all_muts <- tumour_evolution[[2]]

all_muts <- all_muts %>%
  dplyr::mutate(class = ifelse(class=='SV' & gene_hit != 'Other', 'Notch', as.character(class)),
                time = 1 - cell_fraction) %>% 
  dplyr::filter(class != 'Notch')

all_muts$group = factor(all_muts$class, levels=c("SV","SNV", "INDEL"), labels=c("SV","SNV", "INDEL")) 
  

wholegut_muts_combined$group = factor(wholegut_muts_combined$type, levels=c("SV","SNV", "INDEL"), labels=c("SV","SNV", "INDEL")) 

ggplot(wholegut_muts_combined, aes(1-allele_frequency)) +
  geom_density() +
  facet_wrap(~group)

all_muts$assay <- 'tumour'
all_muts_af <- all_muts %>% 
  dplyr::select(sample, group, allele_frequency, assay)

wholegut_muts_combined$assay <- 'wholegut'

wholegut_muts_combined <- wholegut_muts_combined %>%
  dplyr::select(sample, group, allele_frequency, assay)

combined_groups_af <- bind_rows(all_muts_af, wholegut_muts_combined)

med_af <- combined_groups_af %>%
  group_by(assay, group) %>%
  summarise(median=median(allele_frequency))

#### test ####
ggplot(combined_groups_af, aes(1-allele_frequency)) +
  geom_density(aes(fill=assay), alpha=0.5) +
  facet_wrap(~group, nrow=3)

#############







# Allele freqs dist per sample
ggplot(wholegut_muts_combined, aes(1-allele_frequency)) + 
  geom_histogram(aes(fill=type), alpha = 0.7, binwidth = .01) + 
  facet_wrap(~sample, ncol=1)

wholegut_svs_df <- read.delim('data/all_WG_samples_merged_filt_annotated.txt') %>% 
  dplyr::mutate(wg_del = grepl("wg_del", notes)) %>%
  dplyr::filter(!wg_del)

wholegut_snvs_df <- mutationProfiles::getData(infile = 'data/annotated_snvs_wholegut.txt', expression_data = 'data/Buchon_summary_ISCs.txt')
wholegut_indels_df <- mutationProfiles::getData(infile = 'data/annotated_indels_wholegut.txt', expression_data = 'data/Buchon_summary_ISCs.txt', type='indel')

# See overlap with 
wholegut_snvs_df %>% 
  dplyr::filter(paste(chrom, pos, sep = '_') %in% dgrp_snps$key) %>% 
  nrow()

plot_tumour_evolution(all_samples = 'data/all_WG_samples_merged_filt.txt', sample %in% wholegut_samples_n, 
                      tes = FALSE, annotate_with = 'data/dna_damage_merged.bed',
                      all_samples_snvs = 'data/annotated_snvs_wholegut.txt',
                      all_samples_indels = 'data/annotated_indels_wholegut.txt')

snvs <- 'data/annotated_snvs_wholegut.txt'
indels <- 'data/annotated_indels_wholegut.txt'

wholegut_snvs <- mutationProfiles::getData(infile = snvs, expression_data = 'data/Buchon_summary_ISCs.txt', type='snv', sex=='male', chrom %in% chromosomes, attach_info = attach_info)
wholegut_indels_df <- mutationProfiles::getData(infile = indels, expression_data = 'data/Buchon_summary_ISCs.txt', type='indel')

wholegut_svs_df <- read.delim('data/all_WG_samples_merged_filt_annotated.txt') %>% 
  dplyr::mutate(wg_del = grepl("wg_del", notes)) %>%
  dplyr::filter(!wg_del)

old_samples <- c('D050R18', 'D050R20', 'D050R22', 'D050R24')

male_snv_count <- male_snvs %>% 
  filter(af > 0) %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = 'tumour')

wgsnv_count <- wholegut_snvs_df %>% 
  filter(af > 0) %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = ifelse(sample %in% old_samples, 'wholegut_old', 'wholegut_young'))

snv_counts <- bind_rows(male_snv_count, wgsnv_count)
snv_counts$group <- 'snv'

male_indel_count <- male_indels %>% 
  filter(af > 0) %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = 'tumour')

wg_indel_count <- wholegut_indels_df %>% 
  filter(af > 0) %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = ifelse(sample %in% old_samples, 'wholegut_old', 'wholegut_young'))

indel_counts <- bind_rows(male_indel_count, wg_indel_count)
indel_counts$group <- 'indel'

male_sv_count <- all_hits %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = 'tumour')

wg_sv_count <- wholegut_svs_df %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = ifelse(sample %in% old_samples, 'wholegut_old', 'wholegut_young'))

sv_counts <- bind_rows(male_sv_count, wg_sv_count)
sv_counts$group <- 'sv'

combined_muts <-  bind_rows(snv_counts, indel_counts, sv_counts)
combined_muts$group = factor(combined_muts$group, levels=c("sv", "snv", "indel"), labels=c("SV", "SNV", "INDEL")) 

combined_muts %>% 
  dplyr::filter(n < 300) %>% 
  ggplot2::ggplot(., aes(type, n)) + 
  ggplot2::geom_boxplot(aes(fill = type), alpha=0.7) +
  scale_y_continuous("Mutation count") + 
  ggplot2::facet_wrap(~group, scales = 'free') +
  ggplot2::guides(type = FALSE) +
  cleanTheme() +
  ggplot2::theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.title.y =element_text(size=15)
  )


#################
##   WG DELS  ###
#################

wholegut_svs_df <- read.delim('data/all_WG_samples_merged_filt_annotated.txt') %>% 
  dplyr::mutate(wg_del = grepl("wg_del", notes)) %>%
  dplyr::filter(!wg_del)



################
##  SIM data ###
################

attach_info <- 'data/samples_names_conversion.txt'

all_sample_names <- read.delim(attach_info, header = F)
colnames(all_sample_names) <- c("sample", "sample_short", "sample_paper", "sex", "assay")

all_sim_names <- all_sample_names %>% dplyr::filter(grepl("visor", sample))

infile <- '~/Desktop/SV_paper_20/svParser/summary/merged/all_bps_merged.txt'
all_sim_samples <- '~/Desktop/SV_paper_20/svParser/summary/merged/all_samples_merged.txt'

all_hits <- svBreaks::getData(infile=infile, attach_info = attach_info)

error_rates <- 'data/error_rates.txt'
error_rates_df <- read.delim(error_rates, header = T)

error_rates_df$sample = factor(error_rates_df$sample, levels=c("truth", "visorR1", "visorR3", "visorR5", "visorR7", "visorR9", "visorR11", "visorR13", "visorR15", "visorR17", "visorR19"), labels=c("Tr", "s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5"))

total_true_positives <- 8

error_rates_df %>%
  dplyr::filter(sample != 'Tr') %>% 
  dplyr::group_by(sample, type, error_group) %>%
  dplyr::tally() %>% 
  # tidyr::complete(sample, type, error_group) %>% 
  dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>% 
  dplyr::mutate(seq_depth = ifelse(grepl('s', sample), 'shallow', 'deep')) %>%
  ggplot(., aes(sample, n, fill=error_group)) +
  # geom_bar(position = position_dodge(width = .8), width=.8, stat='identity') +
  geom_bar(stat='identity') +
  # geom_hline(yintercept = total_true_positives, linetype="dashed", color = "red") +
  facet_wrap(seq_depth~type, scales = 'free')
           


####
## Plot CNVs along genome
###

devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/cnvPlotteR'))

allPlot(path = '/Volumes/perso/Analysis/Analysis/CNV-Seq/visor/results/w_50000/')

