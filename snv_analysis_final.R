
# library(mutationProfiles)
library(ggpubr)
library(stringr)

excluded_samples <- c("B241R41-2",  "A373R7", "A512R17")
excluded_samples <- c(excluded_samples, "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R09", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
excluded_samples <- c(excluded_samples, "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
excluded_samples <- c(excluded_samples, "D106R1", "D106R3", "D106R5", "D106R7", "D106R9", "D106R11", "D106R13", "D106R15", "D106R17", "D106R19", "D106R21", "D106R23", "D106R25", "D106R27", "D106R29", "D106R31", "D106R33"  )
excluded_samples <- c(excluded_samples, "D197R09", "D197R11", "D197R13", "D197R15")
excluded_samples <- c(excluded_samples, "D265R01", "D265R03", "D265R05", "D265R07", "D265R09", "D265R11", "D265R13")


## GLOBALS ###
rootDir <- ifelse(dir.exists('/Users/Nick_curie/'), '/Users/Nick_curie/', '/Users/Nick/iCloud/')

attach_info <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt')
snvs <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_snvs.txt')
indels <- paste0(rootDir, 'Desktop/script_test/mutationProfiles/data/annotated_indels.txt')

male_snvs <- mutationProfiles::getData(infile = snvs, type='snv', sex=='male', !sample %in% excluded_samples,
                                       assay %in% c("neoplasia", "delta-neoplasia"),
                                       attach_info = attach_info)

male_indels <- mutationProfiles::getData(infile = indels, type='indel', sex=='male', !sample %in% excluded_samples,
                                         assay %in% c("neoplasia", "delta-neoplasia"),
                                         attach_info = attach_info )

male_snvs <- male_snvs %>%
  dplyr::rename(sample_old = sample,
                sample = sample_paper) %>%
dplyr::select(sample, everything())

male_indels <- male_indels %>%
  dplyr::rename(sample_old = sample,
                sample = sample_paper) %>%
  dplyr::select(sample, everything())


# Median SNVs + Ts/Tv ratio
snvStats(snv_data = male_snvs)
# Median INDELS + Ts/Tv ratio
snvStats(snv_data = male_indels)

########################
#####  S. FIG. 6   #####
########################

# S. fig 6 - all SNVs + INDEL per sample

plot_snvs(snv_data = male_snvs)
# plot_snvs(snv_data = male_indels)
indel_lengths(indel_data = male_indels)

# Attach function impact to data

library(dndscv)
source(paste0(rootDir, 'Desktop/script_test/mutationProfiles/script/dnds_analysis.R'))
# path_cds_table = 'data/dm6_cds_table.tsv'
# path_genome_fasta = '/Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa'
# buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "data/dm6_refcds.rda", onlychrs = c("2L", "2R", "3L", "3R", "4", "X", "Y"))

# load("data/dm6_refcds.rda")
# print(RefCDS[[1]])

muts <- combineData(male_snvs, male_indels)

# SNVS + INDELS
all_mutations <- selctionDetect(snv_data = muts)
# save(all_mutations, file = "data/snv_indel_final_dnds")
# all_mutations <- load("data/snv_indel_final_dnds")

# Run per-sample INDELS + SNVS
# per_sample_dNdS <- selctionDetect(per_sample = T, snv_data = muts)
# per_sample_dNdS <- load("data/snv_indel_per_sample_dnds")


########################
## For combined data  ##
########################


# View combined data
quickView(all_mutations)

# Recode INDEL impacts
ann_data <- all_mutations$annotmuts %>%
  dplyr::mutate(impact = ifelse(impact == 'no-SNV',
                                gsub("^.*-", "", ntchange),
                                impact)) %>%
  dplyr::mutate(impact = ifelse(str_detect(impact, 'inframe'), 'Synonymous', impact)) %>%
  dplyr::rename(sample = sampleID)

ann_data <- add_info(ann_data, attach_info = attach_info, paper=TRUE)

snv_data <- male_snvs %>%
  dplyr::select(sample, chrom, pos, af, fpkm, cell_fraction)

indel_data <- male_indels %>%
  dplyr::select(sample, chrom, pos, af, fpkm, cell_fraction)

combined_data <- rbind(snv_data, indel_data)

allele_freqs <- combined_data %>%
  dplyr::select(sample, chrom, pos, af, fpkm, cell_fraction) %>%
  dplyr::rename(chr = chrom)

# Add allele freqs
ann_data <- plyr::join(ann_data, allele_freqs, c("sample", "chr", "pos"), type = "left")

# Write to csv

# write.csv(ann_data, file = '~/Desktop/final_analysis/annotated_mutations_filt.csv', row.names = FALSE)

tally_impacts <- ann_data %>%
  dplyr::filter(!is.na(impact)) %>%
  dplyr::group_by(sample, impact) %>%
  tally() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(total = sum(n))

# Percetnage of mu types
mut_type_summary <- tally_impacts %>%
  dplyr::group_by(impact) %>%
  dplyr::tally() %>%
  dplyr::mutate(perc = plyr::round_any( (100*n/sum(n)), 1)) %>%
  as.data.frame()

p <- ggplot(mut_type_summary)
p <- p + geom_bar(aes(fct_reorder(impact, -perc), perc), stat = 'identity')
p

order <- levels(fct_reorder(tally_impacts$sample, tally_impacts$total))

# Median number of protein coding mutations
median(tally_impacts$total)
table(tally_impacts$impact)


########################
#####  S. FIG. 8   #####
########################

mutSpectrum(snv_data = male_snvs, max_y = 10)





## Fig
# Plot all muts per sample, coloured by impact


## SHare order with SVs by type:
sv_order <- c("B241R45", "B241R63", "B241R53", "HUM-7", "B241R35", "B241R51",
              "B241R49", "B241R39", "A512R21", "B241R55", "HUM-1", "A573R25",
              "B241R43", "A573R29", "B241R61", "A373R5", "A573R33", "B241R37",
              "A573R27", "A573R31", "D106R23", "D106R31", "A512R23", "A373R9",
              "B241R57", "A373R11", "D106R29", "B241R59", "D106R27", "D106R33",
              "B241R41-1", "A373R13", "A373R3", "D106R25", "A373R1")
## Fig
# Plot protein-coding mutations per sample, coloured by impact
ggbarplot(tally_impacts, "sample", "n",
          fill = "impact", color = "impact",  alpha = 0.6, order = order, orientation = "horiz",
          lab.hjust = -1, lab.vjust = 0.5, title="Protein-coding mutations", ylab="Number of Mutations", xlab=T) + theme_minimal()









# Allele Freqs per impact. Combine all non-syn muts
snv_af_data <- ann_data %>%
  dplyr::filter(!is.na(impact)) %>%
  dplyr::mutate(impact = ifelse(impact == 'Synonymous', impact, 'Non-synonymous'))


red <- "#FC4E07"
blue <- "#00AFBB"
yellow <- "#E7B800"
grey <- "grey"


gghistogram(snv_af_data, x = "cell_fraction", color = "impact", fill= "impact", palette = c(blue, 'grey'), alpha = 0.5,
            add = "median", bins=50,
            title="VAF distribution", ylab="Number of Mutations", xlab = "Variant Allele Frequency")


ann_data



#####
# filt_snvs <- read.delim('~/Desktop/snvs_per_sample.txt')
# colnames(filt_snvs) <- 'sample'
#
# filt_snvs <- separate(data = filt_snvs, col = sample, into = c("sample", "count"), sep = "\\:")
#
#
# filt_snvs <- filt_snvs %>%
#   dplyr::mutate(sample = gsub("\\_.*", "", sample))
#
#
# filt_snvs$count <- as.numeric(filt_snvs$count)
# median(filt_snvs$count)
#
# ann_data_filt <- add_info(filt_snvs)
#
#
# ggbarplot(ann_data_filt, "sample", "count", sort.val = 'desc',
#           fill = "sex", color = "sex", palette = c(yellow, blue, red, 'grey'), alpha = 0.6,
#           x.text.angle = 90, title="Consensus SNVs post high-stringency filtering per sample", ylab="Number of Mutations", xlab=F)






########################
#### For per-sample data
########################

# data <- compressTable(per_sample_dNdS)

ann_data <- add_info(per_sample_dNdS)

red <- "#FC4E07"
blue <- "#00AFBB"
yellow <- "#E7B800"

ggdotchart(ann_data[ann_data$name=='wall',], x = "sample", y = "mle",
           color = "sex",
           # group = "sex",
           palette = c(red, blue, yellow),
           sorting = "descending",
           add = "segments",
           dot.size = 4,
           title="dN/dS",
           ggtheme = theme_pubr() ) + theme_cleveland()








# plot snvs per smaple

plot_snvs(snv_data = male_snvs)


# rainfall per sample

for(s in levels(male_snvs$sample)){
  print(s)
  d <- male_snvs %>%
    dplyr::filter(sample == s) %>%
    droplevels()
  tryCatch({
    cat("Plotting rainfall plot for sample", s, "\n")
    print(rainfall(snv_data=d, title=s, write=T))
  }, error=function(e){})
}
