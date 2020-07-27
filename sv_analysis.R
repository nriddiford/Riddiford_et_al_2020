### INIT START #####

library(devtools)
# devtools::install_github("nriddiford/svBreaks", force = T)
# devtools::install_github("nriddiford/mutationProfiles", force = T)

devtools::load_all(path = '~/Desktop/script_test/svBreaks')
devtools::load_all(path = '~/Desktop/script_test/mutationProfiles')

library(svBreaks)
library(mutationProfiles)

library(tidyverse)
library(stringr)
library(bedr)

dfsObj <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/muts.RData")
load(dfsObj)
# rootDir <- ifelse(dir.exists('/Users/Nick_curie/'), '/Users/Nick_curie/', '/Users/Nick/iCloud/')
rootDir <- getwd()

source('R/misc_funs.R')

### INIT END #####


# devtools::install(paste0(rootDir, 'Desktop/script_test/svBreaks'))
# devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/svBreaks'))
# devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/mutationProfiles'))

# bedDir = paste0(rootDir, 'Desktop/misc_bed')
# all_samples <- paste0(rootDir, 'Desktop/parserTest/all_combined_23719/all_samples_merged.txt')

all_samples <- 'data/all_samples_merged.txt'

excluded_events <- read.delim(all_samples) %>% 
  dplyr::filter(status == 'FALSE') %>% 
  dplyr::mutate(exclude = paste(sample, event, sep = '_')) %>%  
  dplyr::pull(exclude)

red <- "#FC4E07"
blue <- "#259FBF"
yellow <- "#E7B800"
grey <- "#333333"
# green <- "#7DBF75"
green <- '#61C456F3'
chromosomes <- c("2L", "2R", "3L", "3R", "X", "Y", "4")

######
######

## Need to update this to include INVersion in notch...
# infile <- paste0(rootDir, 'Desktop/parserTest/all_combined_23719/all_bps_merged.txt')
infile <- 'data/all_bps_merged.txt'

# attach_info <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt')
attach_info <- 'data/samples_names_conversion.txt'
all_sample_names <- read.delim(attach_info, header = F)
colnames(all_sample_names) <- c("sample", "sample_short", "sample_paper", "sex", "assay")

excluded_samples <- c("B241R41-2",  "A373R7", "A512R17")
excluded_samples <- c(excluded_samples, "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R09", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
excluded_samples <- c(excluded_samples, "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
excluded_samples <- c(excluded_samples, "D106R1", "D106R3", "D106R5", "D106R7", "D106R9", "D106R11", "D106R13", "D106R15", "D106R17", "D106R19", "D106R21", "D106R23", "D106R25", "D106R27", "D106R29", "D106R31", "D106R33"  )
excluded_samples <- c(excluded_samples, "D197R09", "D197R11", "D197R13", "D197R15")
excluded_samples <- c(excluded_samples, "D265R01", "D265R03", "D265R05", "D265R07", "D265R09", "D265R11", "D265R13")

all_sample_names <- all_sample_names %>% dplyr::filter(!sample %in% excluded_samples)

all_hits <- svBreaks::getData(!sample %in% excluded_samples,
                              !paste(sample, event, sep = '_') %in% excluded_events,
                              infile=infile, attach_info = 'data/samples_names_conversion.txt')

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

notch <- all_hits %>% 
  dplyr::mutate(key2 = paste0(sample, "_", event)) %>% 
  dplyr::group_by(sample, event) %>% 
  dplyr::filter(key2 %in% n_hits$key) %>% 
  droplevels() %>% 
  dplyr::select(-key2) %>% 
  ungroup()

# snvs <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_snvs.txt')
# indels <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_indels.txt')
snvs <- 'data/annotated_snvs.txt'
indels <- 'data/annotated_indels.txt'

male_snvs <- mutationProfiles::getData(infile = snvs, type='snv', sex=='male', !sample %in% excluded_samples,
                                       assay %in% c("neoplasia", "delta-neoplasia"), chrom %in% chromosomes,
                                       attach_info = attach_info)

male_indels <- mutationProfiles::getData(infile = indels, type='indel', sex=='male', !sample %in% excluded_samples,
                                         assay %in% c("neoplasia", "delta-neoplasia"), chrom %in% chromosomes,
                                         attach_info = attach_info )

combined_muts <- plyr::join(male_snvs, male_indels, type='full') %>% 
  dplyr::rename(allele_frequency = af, chr = chrom)

combined_muts <- addPurity(df=combined_muts) %>% 
  dplyr::rename(af = allele_frequency, chrom = chr)

dfsObj <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/muts.RData")
# save(male_snvs, male_indels, combined_muts, excluded_samples, excluded_events, all_sample_names, notch, non_notch, file = dfsObj)
save(file = dfsObj, all_samples, excluded_events, chromosomes, infile, attach_info, all_sample_names, excluded_samples, all_hits, non_notch, notch, male_snvs, male_indels, combined_muts)


### ANALYSIS START #####

########################
######   FIG. 2   ######
########################

all_samples_df <- read.delim(all_samples, header = T) %>% 
  dplyr::filter(!sample %in% excluded_samples)

all_samples_df <- plyr::join(all_samples_df, all_sample_names, "sample", type = 'left') %>% 
  dplyr::rename(sample_old = sample,
                sample = sample_paper) %>% 
  dplyr::select(sample, everything())

non_inactivating_events <- c('D1_63', 'P7_2148')

all_samples_df <- all_samples_df %>% 
  dplyr::filter(!paste(sample, event, sep = '_') %in% non_inactivating_events)


# Fig 2a: SVs affecting Notch
svBreaks::geneHit(!sample %in% excluded_samples, all_samples_df = all_samples_df, show_sample = T, drivers = c("N"), plot=T)

# Supp Table 1 - Notch events
nnn <- svBreaks::geneHit(!sample %in% excluded_samples, all_samples_df = all_samples_df, show_sample = T, drivers = c("N", "kuz"), plot=F) %>% 
  dplyr::mutate(notch_event = paste(sample, event, sep = '_')) %>%
  dplyr::pull(notch_event)
  
all_samples_df %>% 
  dplyr::mutate(key = paste(sample, event, sep = '_')) %>%
  dplyr::filter(key %in% nnn) %>% 
  dplyr::select(-sample_old, -key, -sample_short) %>% 
  write.table(., file = 'tables/S_table_1.txt', quote=FALSE, sep='\t', row.names = FALSE)


sv_count <- all_hits %>% 
  dplyr::filter(bp_no=='bp1') %>% 
  dplyr::distinct(sample, event, .keep_all=TRUE) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::tally()

sum(sv_count$n)
median(sv_count$n)

# Supp Table 2 - gw events
all_samples_df %>% 
  dplyr::select(-sample_old) %>% 
  dplyr::mutate(status = tidyr::replace_na(status, '')) %>% 
  write.table(., file = 'tables/S_table_2.txt', quote=FALSE, sep='\t', row.names = FALSE)

# Supp Table 5 - Excluded CN events
all_samples_df %>% 
  dplyr::mutate(status = tidyr::replace_na(status, '')) %>% 
  dplyr::mutate(key = paste(sample_old, event, sep = '_')) %>%
  dplyr::filter(key %in% excluded_events) %>% 
  dplyr::select(-sample_old, -key, -sample_short, -sex, -assay) %>% 
  write.table(., file = 'tables/S_table_5.txt', quote=FALSE, sep='\t', row.names = FALSE)


# Fig 2b: Svs in Notch region
svBreaks::notchHits(!sample %in% excluded_samples, all_samples_df = all_samples_df, show_samples = T, bp_density = T)

# S.fig 3
svBreaks::notchHits(!sample %in% excluded_samples, all_samples_df = all_samples_df, from=3.12, to=3.173, show_samples = T, bp_density = T, adjust = 0.1, ticks = 2.5)
# svBreaks::notchHits(!sample %in% excluded_samples, all_samples = all_samples, from=3.12, to=3.173, show_samples = F, bp_density = T, ticks = 5)

## Write breakpoints +/- 5kb for Notch SVs
n_hits <- svBreaks::geneHit(!sample %in% excluded_samples, plot = F, all_samples_df = all_samples_df, drivers = c("N"))

notch_tss <- 3.135669 * 1e6
dist <- 2000

n_hits %>% 
  dplyr::mutate(distbp1 = abs(start - notch_tss) < dist,
                distbp2 = abs(end - notch_tss) < dist) %>% 
  dplyr::select(sample, event, distbp1, distbp2) %>% 
  dplyr::filter(distbp1 | distbp2) 

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
  svBreaks::writeBed(df=., outDir = 'data', name = paste0('notch_breakpoints_', expand_by, '_', today, '.bed'))

# Run `bash extract_breakpoints.sh -i ../notch_breakpoints_5k_10.3.20.bed -n 100`

# Calculate median NBR length:

read.delim(paste0(rootDir, 'Desktop/misc_bed/breakpoints/Notch_CFS/notch_breakpoint_regions_5k.bed'), header = F) %>% 
  dplyr::mutate(length = V3 - V2) %>% 
  dplyr::summarise(median_length = median(length, na.rm = TRUE))

mappable_regions = 'data/dmel6_mappable.bed'
chrom_lengths = 'data/chrom.sizes.txt' 
overlaps <- bpRegioneR(regionA = 'data/notch_breakpoint_regions_500_mappable.bed',
           regionB = 'data/motif_1.mappable.bed',
           mappable_regions = mappable_regions,
           n=100, from=2700000, to=3400000, chrom = 'X', plot=F)

# Fig 2d: Enrichment of motif1 in notch breakpoint regions
plot_bpRegioner(overlaps, bins = 10)

# G4s
bpRegioneR(regionA = 'data/notch_breakpoint_regions_500_mappable.bed',
           regionB = 'data/g4s.mappable.bed', mappable_regions = mappable_regions,
           n=100, from=2700000, to=3400000, chrom = 'X', plot=F)


# Fig 2e,f: Two characterised events in Notch 

# # Fig 2f: TEs at Notch breakpoints
# te_in <- paste0(rootDir, 'Documents/Curie/Documents/SV_paper/sv_characterisation/TE_involvement.txt')
# notch_tes <- read.delim(te_in)
# notch_tes <- notch_tes %>% 
#   dplyr::filter(!sample %in% excluded_samples) %>% 
#   dplyr::select(sample, bp1, bp2, type) %>% 
#   dplyr::mutate(type2 = type, 
#                 bp1 = as.logical(ifelse(stringr::str_detect(bp1, ''), 1, 0)),
#                 bp2 = as.logical(ifelse(stringr::str_detect(bp2, ''), 1, 0))) %>% 
#   dplyr::mutate(has_te = ifelse(bp1 | bp2, TRUE, FALSE)) %>% 
#   dplyr::group_by(type2) %>%
#   dplyr::count(type2, has_te) %>% 
#   dplyr::mutate(total = sum(n), perc = (n/total)*100) %>% 
#   dplyr::filter(has_te) %>% 
#   as.data.frame()

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
gw_sv_types <- non_notch %>%
  dplyr::mutate(type2 = factor(type2)) %>% 
  dplyr::group_by(type2) %>%
  dplyr::distinct(chrom, event, .keep_all=TRUE) %>%
  dplyr::tally() %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1),
                group = 'genome-wide') %>% 
  as.data.frame()

notch_sv_types <- notch %>% 
  dplyr::filter(bp_no=='bp1') %>% 
  dplyr::mutate(type2 = factor(type2)) %>% 
  dplyr::group_by(type2) %>%
  dplyr::distinct(event, .keep_all=TRUE) %>%
  dplyr::tally() %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1),
                group = 'Notch') %>% 
  as.data.frame()

combined_sv_types <- plyr::join(gw_sv_types, notch_sv_types, type='full')

order <- levels(fct_reorder(gw_sv_types$type2, -gw_sv_types$perc))

combined_sv_types$type2 <- factor(combined_sv_types$type2, levels = order, labels=c("Translocation","Inversion", "Duplication", "Complex", "Deletion", "Tandem Duplication"))

combined_sv_types <- complete(combined_sv_types, group, type2, fill = list(n=0, perc = 0))

cols <- c("genome-wide" = grey, "Notch" = green)

# Fig 4b: SV types
combined_sv_types %>% 
  ggplot2::ggplot(.) +
  ggplot2::geom_bar(aes(type2, perc, fill = group), colour = 'transparent', alpha=0.7, stat='identity', position = 'dodge') + 
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


# noR1 <- non_notch %>% dplyr::filter(sample != "A373R1")
# noR1 <- transform_types(noR1)
# noR1_samples <- c("A373R1", excluded_samples)
# Get SV tpes per sample

sv_data <- svTypes(!sample %in% excluded_samples, bp_data = non_notch, plot = F)

sv_df <- sv_data[[1]]

sv_df2 <- swapSampleNames(df=sv_df) %>% 
  dplyr::filter(sample_old != 'A373R1')

levels(sv_df2$sample)

sv_df2 %>% 
  dplyr::group_by(type2) %>% 
  dplyr::summarise(Freq = sum(type_count))


sv_df2$sample <- fct_reorder(sv_df2$sample, -sv_df2$sv_count)

# Sup fig 4a - af of mutations by type
gghistogram(non_notch, x = "af", color = "type2", fill= "type2", alpha = 0.5,
            add = "median", bins=50, title="VAF distribution of SV types",
            ylab="Number of SVs", xlab = "VAF", ggtheme = theme_minimal())

# Fig 4c - Protein coding mutations (non-Notch)
snv_indel_df <- read.csv('data/annotated_mutations_filt.csv')

# Supp Table 3 - Protein coding mutations
snv_indel_df %>% 
  dplyr::select(-sample_old, -sample_short, -assay, -sex) %>% 
  write.table(., file = 'tables/S_table_3.txt', quote=FALSE, sep='\t', row.names = FALSE)

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
                sample_old != 'A373R1',
                !sample_old %in% excluded_samples) %>%
  dplyr::mutate(impact = ifelse(str_detect(impact, 'inframe'), 'Synonymous', as.character(impact))) %>%
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

# tally_impacts %>% 
#   dplyr::group_by(impact) %>% 
#   dplyr::summarise(total = sum(type_count)) %>% 
#   ggplot2::ggplot(.) + 
#   ggplot2::geom_bar(aes(fct_reorder(impact, -total), total), colour = grey, fill = grey, alpha=0.6, stat = 'identity') +
#   cleanTheme() +
#   theme(
#     axis.title.y = element_text(size=18),
#     panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=15),
#     axis.title.x = element_blank(),
#     legend.position = 'none'
#   )


colours <- sv_colours()

p1 <- ggplot(sv_df2)
p1 <- p1 + geom_bar(aes(fct_reorder(sample, -sv_count), type_count, colour = 'transparent', fill = type2), alpha=0.7, stat = 'identity')
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
p2 <- p2 + geom_bar(aes(sample, type_count, colour = 'transparent', fill = impact), alpha=0.7, stat = 'identity')
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

# Fig4c/d mutations per sanmple
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

source('R/enrichment_funs.R')

####################
## SNVS + INDELS  ##
####################

## Fig 5 S. 1 - muts in genes by expression
r <- hits_in_genes(sample != "A373R1", df = non_notch, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), out_file='bps_in_genes.txt')
bps_in_genes <- r[[2]]
bps_in_genes$mut <- 'sv'

r <- hits_in_genes(sample != "A373R1", df = male_snvs, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), out_file='snvs_in_genes.txt')
snvs_in_genes <- r[[2]]
snvs_in_genes$mut <- 'snv'

r <- hits_in_genes(sample != "A373R1", df = male_indels, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), out_file='indels_in_genes.txt')
indels_in_genes <- r[[2]]
indels_in_genes$mut <- 'indel'

comb <- rbind(bps_in_genes, snvs_in_genes, indels_in_genes)

comb <- comb %>%
  tidyr::complete(group, feature, mut) %>%
  dplyr::mutate(Log2FC = ifelse(is.na(Log2FC), 0, Log2FC),
                label = ifelse(sig %in% c(NA,'-'), '', sig)) %>%
  dplyr::mutate(observed = replace_na(observed, 0)) %>%
  as.data.frame()


comb$mut = factor(comb$mut, levels=c("sv","snv", "indel"), labels=c("SV","SNV", "INDEL"))
# comb$group = factor(comb$group, levels=c("All expressed genes", "Top 20% FPKM", "Non-expressed genes"), labels=c("All expressed genes", "Top 20% FPKM", "Non-expressed genes"))
comb$group = factor(comb$group, levels=c("All expressed genes", "Non-expressed genes"), labels=c("ISC-expressed genes", "Non-expressed genes"))

# Figure 5 a - mutations in gene features
gf_exp_grouped <- ggbarplot(comb, x = "feature", y = "Log2FC",
                            fill = "mut", color = "transparent", group = "group", palette =c(blue, yellow, red, grey), alpha = 0.6,
                            label = 'label', lab.vjust = -.1,
                            orientation = "horizontal", x.text.angle = 90, position = position_dodge(0.8),
                            xlab = FALSE, ylab = "Log2(FC)", ggtheme = theme_pubr())

gf_exp_grouped <- ggpar(gf_exp_grouped, ylim = c(-2.2,2.2))
facet(gf_exp_grouped, facet.by = 'group', nrow = 1)


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

# Fig5 a - mutations in gene features
gf_exp_grouped <- ggbarplot(comb, x = "feature", y = "Log2FC",
                            fill = "group", color = "transparent", group = "group", palette =c(blue, yellow, red, grey), alpha = 0.6, 
                            label = 'label', lab.vjust = 2,
                            orientation = "horizontal", x.text.angle = 90, position = position_dodge(0.8), 
                            xlab = FALSE, ylab = "Log2(FC)", ggtheme = theme_pubr())

ggpar(gf_exp_grouped, ylim = c(-2.2,2.2))

facet(gf_exp_grouped, facet.by = 'group', nrow = 1)


# Fig5 d Modencode enriched
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/enriched/merged/'))
plot_all(d)

# Manon DamID
d <- run_all(sample != "A373R1", bps = non_notch, snvs = male_snvs, indels = male_indels, dir = paste0(bedDir, '/Manon_2020/deseq/'))
plot_all(d)

# Fig5 b repeats
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


dist <- run_reldist(a = all_mutations, b = 'data/repeats_merged.bed',
                    chroms = chromosomes, verbose = FALSE, type='snv', print = FALSE, sim=T)

plot_reldist(dist)



########################
######   FIG. 6   ######
########################

# Fig 6c & d
source('R/plotFreq.R')

samples_with_no_N_sv <- c(excluded_samples, 'A373R9', 'B241R43', 'B241R49', 'D197R05')

read.delim(all_samples) %>% 
  dplyr::filter(!paste(sample, event, sep = '_') %in% excluded_events) %>% 
  write.table(., file = 'data/all_samples_merged_filt.txt', quote=FALSE, sep='\t', row.names = FALSE)

plot_tumour_evolution(all_samples = 'data/all_samples_merged_filt.txt', !sample %in% samples_with_no_N_sv, 
                      tes = FALSE, annotate_with = 'data/dna_damage_merged.bed')

tumour_evolution <- plot_tumour_evolution(all_samples = 'data/all_samples_merged_filt.txt', 
                                          !sample %in% samples_with_no_N_sv, 
                                          tes = FALSE, annotate_with = 'data/dna_damage_merged.bed', plot=F)

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
  scale_y_continuous('Cumulative mutations') +
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


addPurity <- function(df, purity_file='data/tumour_purity.txt'){
  purity <- read.delim(purity_file, header = F)
  colnames(purity) <- c('sample', 'purity')
  
  df <- plyr::join(df, purity, by='sample')
  df$corrected_af <- df$allele_frequency + (1-df$purity) * df$allele_frequency
  
  
  df <- df %>% 
    dplyr::mutate(corrected_af = ifelse(corrected_af>1, 1, corrected_af)) %>% 
    dplyr::rename(af_old = allele_frequency,
                  allele_frequency = corrected_af) %>% 
    dplyr::mutate(cell_fraction = ifelse(chr %in% c('X', 'Y'), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1, allele_frequency*2)))
  
  return(df)
}

# Pre and post Notch hits
snv_indel_df <- read.csv('data/annotated_mutations_filt.csv') %>% 
  dplyr::rename(sample_new = sample,
                sample = sample_old,
                allele_frequency = af) 

# Correct af and cell fraction for tumour purity
snv_indel_df <- addPurity(df=snv_indel_df) %>% 
  dplyr::rename(sample_old = sample,
                sample = sample_new) 

snv_indel_df %>% 
  dplyr::select(-c(sample_short:af_old, purity)) %>% 
  write.csv(., file = 'data/annotated_mutations_filt_pur_adj.csv', row.names = FALSE)

notch_svs <- notch_hits %>% dplyr::select(sample, highest_n)

coding2notch <- plyr::join(snv_indel_df, notch_svs, c("sample"), type = "left") %>% 
  dplyr::filter(!sample_old %in% c('A373R1',samples_with_no_N_sv)) %>% 
  dplyr::mutate(clonality = as.factor(ifelse(cell_fraction > highest_n, 'pre', 'post')),
                time = 1 - cell_fraction) 

clonality_counts <- coding2notch %>%  
  dplyr::group_by(clonality, impact) %>%
  tally() %>% 
  dplyr::filter(!is.na(impact)) %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1))


clonality_counts <- complete(clonality_counts, clonality, impact, fill = list(n=0, perc = 0))

clonality_counts$clonality = factor(clonality_counts$clonality, levels=c("pre", "post"), labels=c("Pre-Notch", "Post-Notch")) 

colours <- snv_colours()

clonality_counts %>% 
  ggplot2::ggplot(., aes(fct_reorder(impact, -perc), perc, group=clonality, fill=impact, label=n)) +
  ggplot2::geom_bar(aes(alpha=clonality), stat='identity', position = 'dodge') +
  ggplot2::scale_alpha_discrete(range = c(0.7, 0.3)) +
  ggplot2::scale_y_continuous("Percentage", limits=c(0,100)) + 
  ggplot2::geom_text(size = 3, position=position_dodge(width=0.9), vjust=-1.25) +
  ggplot2::scale_fill_manual("SV type\n", values = colours) +
  cleanTheme() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=15),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y =element_text(size=15)
  )

gw_point_muts <- combined_muts %>% 
  dplyr::rename(sample_old = sample,
              sample = sample_paper) 

# All muts = coding and non-coding
combined <- plyr::join(gw_point_muts, notch_svs, c("sample"), type = "left")

gw_rel_to_Notch <- combined %>% 
  dplyr::filter(!sample_old %in% c('A373R1', samples_with_no_N_sv)) %>%
  # dplyr::filter(! is.na(grouped_trans)) %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(clonality = ifelse(cell_fraction > highest_n, 'pre', 'post'))

gw_clonality_counts <- gw_rel_to_Notch %>%  
  dplyr::mutate(length = nchar(as.character(ref)) - nchar(as.character(alt)),
                type = ifelse(!is.na(grouped_trans), 'snv', 
                              ifelse(length < 0, 'ins', 'del'))) %>% 
  # dplyr::select(sample, type, ref, alt, clonality) %>% View() 
  dplyr::group_by(clonality, type) %>%
  tally() %>% 
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1))

gw_clonality_counts <- complete(gw_clonality_counts, clonality, type, fill = list(n=0, perc = 0))

gw_clonality_counts$clonality = factor(gw_clonality_counts$clonality, levels=c("pre", "post"), labels=c("Pre-Notch", "Post-Notch")) 
gw_clonality_counts$type = factor(gw_clonality_counts$type, levels=c("del", "snv", "ins"), labels=c("Deletion", "SNV", "Insertion")) 

gw_clonality_counts %>% 
  ggplot2::ggplot(., aes(fct_reorder(type, -perc), perc, group=clonality, fill=type, label=n)) +
  ggplot2::geom_bar(aes(alpha=clonality), stat='identity', position = 'dodge') +
  ggplot2::scale_alpha_discrete(range = c(0.7, 0.3)) +
  ggplot2::scale_y_continuous("Percentage", limits=c(0,100)) + 
  ggplot2::geom_text(size = 3, position=position_dodge(width=0.9), vjust=-1.25) +
  cleanTheme() +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=15),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y =element_text(size=15)
  )
