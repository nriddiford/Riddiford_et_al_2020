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
  

samples_with_no_N_sv <- c(excluded_samples, 'A373R9', 'B241R43', 'B241R49', 'D197R05')

read.delim(all_samples) %>% 
  dplyr::filter(sample %in% wholegut_samples) %>% 
  write.table(., file = 'data/all_WG_samples_merged_filt.txt', quote=FALSE, sep='\t', row.names = FALSE)

plot_tumour_evolution(all_samples = 'data/all_WG_samples_merged_filt.txt', sample %in% wholegut_samples, 
                      tes = FALSE, annotate_with = 'data/dna_damage_merged.bed')

