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


red <- "#FC4E07"
blue <- "#259FBF"
yellow <- "#E7B800"
grey <- "#333333"
# green <- "#7DBF75"
green <- '#61C456F3'

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
wholegut_svs_df <- read.delim('data/all_WG_samples_merged_filt_annotated.txt') %>%
  dplyr::mutate(wg_del = grepl("wg_del", notes),
                cell_fraction = ifelse(chromosome1 %in% c('X', 'Y'), allele_frequency,
                                       ifelse(allele_frequency*2>1, 1, allele_frequency*2))) %>%
  dplyr::filter(!wg_del) 


wholegut_young_samples <- c('D050R10', 'D050R12', 'D050R14', 'D050R16')
## Look at wholgut DUPS

wg_sv_types <- wholegut_svs_df %>% 
  dplyr::filter(!sample %in% c(wholegut_young_samples, 'D050R24')) %>% 
  dplyr::mutate(allele_frequency = as.double(as.character(allele_frequency))) %>% 
  dplyr::mutate(cell_fraction = ifelse(chromosome1 %in% c('X', 'Y'), allele_frequency,
                                       ifelse(allele_frequency*2>1, 1, allele_frequency*2)),
                assay='wholegut',
                length = length.Kb.) %>% 
  dplyr::group_by(type) %>% 
  dplyr::tally() %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1),
                group = 'wholegut') %>% 
  dplyr::select(type, group, n, perc) 
  

gw_sv_types <- non_notch %>%
  dplyr::filter(!sample %in% c('A373R1', excluded_samples)) %>% 
  dplyr::mutate(type2 = factor(type2)) %>% 
  dplyr::group_by(type2) %>%
  dplyr::distinct(chrom, event, .keep_all=TRUE) %>%
  dplyr::tally() %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1),
                group = 'genome-wide',
                type = type2) %>% 
  dplyr::select(-type2, type, group, n, perc) %>% 
  as.data.frame() %>% 
  droplevels()


combined_sv_types <- plyr::join(wg_sv_types, gw_sv_types, type='full') %>% 
  dplyr::group_by(group) %>% 
  tidyr::complete(type, fill = list(n=0, perc=0))
  
combined_sv_types$type <- factor(combined_sv_types$type, levels = c("TRA", "BND", "DUP", "COMPLEX", "DEL", "TANDUP"), labels=c("Translocation","Inversion", "Duplication", "Complex", "Deletion", "Tandem Duplication"))

cols <- c("genome-wide" = grey, "wholegut" = red)

# Fig 4b: SV types
combined_sv_types %>% 
  ggplot2::ggplot(.) +
  ggplot2::geom_bar(aes(type, perc, fill = group), colour = 'transparent', alpha=0.7, stat='identity', position = 'dodge') + 
  ggplot2::scale_y_continuous("Percentage", limits=c(0,100)) + 
  ggplot2::scale_fill_manual("Group\n", values = cols) +
  cleanTheme() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=15),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y =element_text(size=15)
  )


## See dup lengths
wg_dups <- wholegut_svs_df %>% 
  dplyr::filter(!sample %in% c(wholegut_young_samples, 'D050R24')) %>% 
  dplyr::mutate(allele_frequency = as.double(as.character(allele_frequency))) %>% 
  dplyr::mutate(assay='wholegut',
                length = length.Kb.) %>% 
  dplyr::filter(type == 'DUP')

  ggplot2::ggplot(.) +
  ggplot2::geom_density(aes(length))

tumour_dups <- all_hits %>% 
  dplyr::filter(!sample %in% c('A373R1', excluded_samples)) %>% 
  dplyr::mutate(allele_frequency = as.double(as.character(af))) %>% 
  dplyr::mutate(assay='tumour') %>% 
  dplyr::filter(type == 'DUP')
  

bind_rows(wg_dups, tumour_dups) %>% 
  ggplot2::ggplot(.) +
  ggplot2::geom_histogram(aes(length, fill=assay, colour=assay), binwidth = 1, alpha = 0.5, position='identity') +
  ggplot2::geom_vline(xintercept = 5, linetype="dashed", color = "black")


## Allele freqs ##

wholegut_svs_af <- wholegut_svs_df %>% 
  dplyr::mutate(type = 'SV') %>% 
  dplyr::select(sample, type, cell_fraction)

wholegut_snvs_af <- read.delim('data/annotated_snvs_wholegut.txt') %>% 
  dplyr::mutate(type = 'SNV',
                cell_fraction = ifelse(chromosome %in% c('X', 'Y'), allele_frequency,
                                       ifelse(allele_frequency*2>1, 1, allele_frequency*2))) %>% 
  dplyr::select(sample, type, cell_fraction)

wholegut_indels_af <- read.delim('data/annotated_indels_wholegut.txt') %>% 
  dplyr::mutate(type = 'INDEL',
                cell_fraction = ifelse(chromosome %in% c('X', 'Y'), allele_frequency,
                                       ifelse(allele_frequency*2>1, 1, allele_frequency*2))) %>% 
  dplyr::select(sample, type, cell_fraction)

wholegut_muts_combined <- dplyr::bind_rows(wholegut_svs_af, wholegut_snvs_af, wholegut_indels_af) %>% 
  dplyr::filter(cell_fraction > 0)


##### Compare real data allele frequency dist to whole gut

tumour_evolution <- plot_tumour_evolution(all_samples = 'data/all_samples_merged_filt.txt', 
                                          tes = FALSE, annotate_with = 'data/dna_damage_merged.bed', plot=F)

all_muts <- tumour_evolution[[2]]

all_muts <- all_muts %>%
  dplyr::filter(!sample_old %in% c('A373R1', excluded_samples)) %>% 
  dplyr::mutate(time = 1 - cell_fraction) 

all_muts$group = factor(all_muts$class, levels=c("SV","SNV", "INDEL"), labels=c("SV","SNV", "INDEL")) 
wholegut_muts_combined$group = factor(wholegut_muts_combined$type, levels=c("SV","SNV", "INDEL"), labels=c("SV","SNV", "INDEL")) 

all_muts$assay <- 'tumour'
all_muts_af <- all_muts %>% 
  dplyr::select(sample, group, cell_fraction, assay)

wholegut_muts_combined$assay <- 'wholegut'

wholegut_muts_combined <- wholegut_muts_combined %>%
  dplyr::select(sample, group, cell_fraction, assay)

combined_groups_af <- bind_rows(all_muts_af, wholegut_muts_combined) %>% 
  dplyr::mutate(cell_fraction = ifelse(cell_fraction == 1, 0.99, cell_fraction)) %>% 
  dplyr::mutate(time = 1 - cell_fraction) %>% 
  dplyr::filter(cell_fraction > 0)


med_af <- combined_groups_af %>%
  dplyr::filter(group == 'INDEL') %>% 
  dplyr::group_by(assay, group) %>%
  dplyr::summarise(median = median(time))


#### test ####
combined_groups_af %>% 
  dplyr::filter(group == 'INDEL') %>% 
  ggplot2::ggplot(., aes(time)) +
  ggplot2::geom_histogram(aes(fill=assay, colour=assay), binwidth = 0.01, alpha = 0.3, position='identity') +
  ggplot2::geom_vline(data=med_af, aes(xintercept=median, colour=assay),
             linetype="dashed", size=0.5)


snvs <- 'data/annotated_snvs_wholegut.txt'
indels <- 'data/annotated_indels_wholegut.txt'

wholegut_snvs_df <- mutationProfiles::getData(infile = snvs, expression_data = 'data/Buchon_summary_ISCs.txt', type='snv', sex=='male', chrom %in% chromosomes, attach_info = attach_info)
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
  dplyr::filter(!sample %in% c('A373R1', 'D050R24')) %>% 
  dplyr::filter(type != 'wholegut_old') %>% 
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

# wholegut_svs_df <- read.delim('data/all_WG_samples_merged_filt_annotated.txt') %>% 
#   dplyr::mutate(wg_del = grepl("wg_del", notes)) %>%
#   dplyr::filter(!wg_del)
# 
# 

################
##  SIM data ###
################

####
## Plot CNVs along genome
###

devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/cnvPlotteR'))

files <- dir("~/Desktop/visorCNV/", pattern = "cnv")

for (f in files){
  n <- strsplit(f, "/")[[1]]
  s <- strsplit(n, "\\.")[[1]][1]
  print(regionPlot(cnv_file=paste0("~/Desktop/visorCNV/",f), from=4950000, to=5099998, bp1=5000000,bp2=5049998,chrom="X", tick=50000, ylim = c(-3,3), title=paste0(s, ": 50.0Kb DEL on X")))
}

#####
## Plot error rates
####

attach_info <- 'data/samples_names_conversion.txt'

all_sample_names <- read.delim(attach_info, header = F)
colnames(all_sample_names) <- c("sample", "sample_short", "sample_paper", "sex", "assay")

all_sim_names <- all_sample_names %>% dplyr::filter(grepl("visor", sample))

infile <- '~/Desktop/SV_paper_20/svParser/visor/summary/merged/all_bps_merged.txt'
all_sim_samples <- '~/Desktop/SV_paper_20/svParser/visor/summary/merged/all_samples_merged.txt'

# all_hits <- svBreaks::getData(infile=infile, attach_info = attach_info)

error_rates <- 'data/error_rates.txt'
error_rates_df <- read.delim(error_rates, header = T)

error_rates_df$sample = factor(error_rates_df$sample, levels=c("visorR1", "visorR3", "visorR5", "visorR7", "visorR9", "visorR11", "visorR13", "visorR15", "visorR17", "visorR19"), labels=c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5"))

total_true_positives <- 8

error_rates_df %>%
  dplyr::group_by(sample, type, error_group)

error_rates_df %>%
  dplyr::group_by(sample, type, depth, error_group) %>%
  dplyr::tally() %>% 
  # tidyr::complete(sample, type, error_group) %>% 
  dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>% 
  # dplyr::mutate(seq_depth = ifelse(grepl('s', sample), 'shallow', 'deep')) %>%
  ggplot(., aes(sample, n, fill=error_group)) +
  # geom_bar(position = position_dodge(width = .8), width=.8, stat='identity') +
  geom_bar(stat='identity') +
  # geom_hline(yintercept = total_true_positives, linetype="dashed", color = "red") +
  facet_wrap(depth~type, scales = 'free')
           
  
install.packages('lemon')
library(lemon)
error_rates_df %>%
  dplyr::group_by(type, error_group, depth, purity, rep) %>% 
  dplyr::tally() %>%
  tidyr::complete(type, error_group, depth, purity, rep) %>% 
  dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
  dplyr::summarise(mean_muts = mean(n), 
                   sd_muts = sd(n)) %>%
  dplyr::mutate(vaf = paste0(purity, '%')) %>%
  ggplot(., aes(type, mean_muts, fill=error_group)) +
  geom_col(alpha=0.8, position = position_dodge(width = .7), width=.7) +
  geom_errorbar(aes(ymin = mean_muts - sd_muts, ymax = mean_muts + sd_muts), position = position_dodge(width = .7), width=0.2) +
  facet_rep_wrap(depth~vaf, nrow = 2, repeat.tick.labels = T) +
  guides(colour = FALSE, fill=FALSE) +
  cleanTheme() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15)
  )


#   p <- p + cleanTheme() +
#     theme(
#       panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
#       axis.title.x = element_blank()
#     )
  

# Quantification of error rates
error_rates_df %>% 
  dplyr::group_by(sample, type, error_group) %>% 
  dplyr::tally() %>% 
  dplyr::ungroup() %>% 
  dplyr::filter( error_group != 'fp') %>%
  dplyr::group_by(sample, type) %>%
  dplyr::mutate(total = sum(n)) %>% 
  dplyr::mutate(perc = (n/total)*100) %>% 
  dplyr::filter(error_group != 'fn') %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(av = mean(perc))






