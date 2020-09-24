# readBed
#'
#' Simple function to read in bedfiles and return as a 3 col df
#' @keywords bed
#' @import dplyr
#' @export
readBed <- function(bed_in = NULL, annotated_file=FALSE){
  if(annotated_file){
    f <- read.delim(bed_in, header=T)
    df <- f %>%
      dplyr::mutate(chrom = chromosome,
                    start = pos,
                    end = start + 1) %>%
      dplyr::select(chrom, start, end)
  } else{
    df <- read.delim(bed_in, header=F)
    if(is.null(df$V3)){
      df$V3 <- df$V2 + 2
    }
    if(!(is.null(df$V4))){
      colnames(df[,c(1,2,3,4,5)]) <- c("chrom", "start", "end", "gene", "id")
      names(df)[1:5] <- c("chrom", "start", "end","gene", "id")
    } else {
      colnames(df[,c(1,2,3)]) <- c("chrom", "start", "end")
      names(df)[1:3] <- c("chrom", "start", "end")
    }
   
  }
  return(df)
}


geneIn <- function(gene, gene_list) {
  # hits <- sapply(as.character(gene_list), function(x) tolower(gene) %in% tolower(strsplit(x, ", ")[[1]]), USE.NAMES=FALSE)
  affected_genes <- gene_list
  genes <- gene
  # if(!strsplit(affected_genes, ", ")[[1]]) return(FALSE)
  i <- which(tolower(genes) %in% tolower(strsplit(affected_genes, ", ")[[1]]))
  if(length(i)){
    l = list()
    for (h in 1:length(i)){
      l[h] = as.character(gene[[h]])
    }
    return(paste(l, collapse = ','))
  }
  return("Other")
}


get_SVs <- function(..., all_samples, bed_file, drivers){
  cat("Reading SV data...\n")
  
  all_data <- read.delim(all_samples, header = T)
  all_data$affected_genes <- as.character(all_data$affected_genes)

  gene_hits <- all_data %>%
    dplyr::filter(
      !status %in% c('F', 'aF')) %>%
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>%
    dplyr::mutate(allele_frequency = as.double(as.character(allele_frequency)),
                  type_decomposed = as.character(ifelse(str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>%
    dplyr::group_by(sample, event) %>%
    dplyr::mutate(affected_genes = paste0(affected_genes, collapse = ", ")) %>% 
    # dplyr::mutate(gene_hit = ifelse(geneIn(drivers, affected_genes),'Notch', 'Other')) %>%
    # dplyr::mutate(special_hit = as.character(ifelse(type %in% c("COMPLEX", "DEL") && any(geneIn(bed_file$gene, affected_genes)), 'Hit', 'Other'))) %>%
    dplyr::mutate(gene_hit = as.character(geneIn(drivers, affected_genes)), "Other") %>% 
    dplyr::mutate(gene_hit = as.character(ifelse(gene_hit != "Other", "Notch", "Other"))) %>% 
    dplyr::mutate(special_hit = ifelse(type %in% c("COMPLEX", "DEL"), as.character(geneIn(bed_file$gene, affected_genes)), "Other")) %>% 
    
    dplyr::mutate(cell_fraction = ifelse(chromosome1 %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1, allele_frequency*2))) %>%
    dplyr::mutate(class = 'SV') %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit) %>%
    droplevels()

  return(as.data.frame(gene_hits))
}

addPurity <- function(df, purity_file='data/tumour_purity.txt'){
  purity <- read.delim(purity_file, header = F)
  colnames(purity) <- c('sample', 'purity')
  
  df <- plyr::join(df, purity, by='sample')
  df$corrected_af <- df$allele_frequency + (1-df$purity) * df$allele_frequency

  df <- df %>% 
    dplyr::mutate(corrected_af = ifelse(corrected_af>1, 1, corrected_af)) %>% 
    dplyr::rename(af_old = allele_frequency,
                  allele_frequency = corrected_af)
  
  return(df)
}


get_SNVs <- function(..., all_samples_snvs, bed_file, drivers, purity_adj){
  cat("Reading SNV data...\n")
  
  snv_data <- read.delim(all_samples_snvs, header = T)
  colnames(snv_data) <- c('sample', 'chromosome',	'pos',	'ref',	'alt',	'trinuc',	'trans',	'decomp_trinuc',	'grouped_trans',	'allele_frequency',	'caller',	'variant_type',	'status',	'snpEff_anno',	'feature',	'gene',	'id')
  
  if(purity_adj) snv_data <- addPurity(df = snv_data)
  
  snv_hits <- snv_data %>% 
    dplyr::filter(...) %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1, allele_frequency*2)),
                  class="SNV") %>%
    dplyr::mutate(gene_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & gene %in% drivers, 'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & tolower(gene) %in% bed_file$gene, as.character(gene), 'Other'))) %>% 
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit) %>% 
    droplevels()
  
  return(as.data.frame(snv_hits))
}


get_INDELs <- function(..., all_samples_indels, bed_file, drivers, purity_adj){
  cat("Reading INDEL data...\n")
  
  indel_data <- read.delim(all_samples_indels, header = T)
  # colnames(snv_data) <- c('sample', 'chromosome',	'pos',	'ref',	'alt',	'trinuc',	'trans',	'decomp_trinuc',	'grouped_trans',	'allele_frequency',	'caller',	'variant_type',	'status',	'snpEff_anno',	'feature',	'gene',	'id')
  
  if(purity_adj) indel_data <- addPurity(df = indel_data)
  
  indel_hits <- indel_data %>% 
    dplyr::filter(...) %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1,allele_frequency*2)),
                  class="INDEL") %>%
    # dplyr::group_by(sample, chromosome, pos) %>% 
    dplyr::mutate(gene_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & gene %in% drivers, 'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & tolower(gene) %in% bed_file$gene, as.character(gene), 'Other'))) %>% 
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit) %>% 
    droplevels()
  
  return(as.data.frame(indel_hits))
}


get_TEs <- function(..., all_samples_tes, bed_file, drivers){
  cat("Reading TE data...\n")
  te_data <- read.delim(all_samples_tes, header = F)
  
  colnames(te_data) <- c("sample", "bp_no", "caller", "genotype", "chromosome1", "bp1", "gene1", "feature1", "chromosome2", "bp2", "gene2", "feature2", "type", "somefield", "allele_frequency", "confidence")
  
  te_data <- te_data %>% 
    dplyr::filter(...,
                  bp_no == 'bp1') %>% 
    tidyr::separate(sample, into = c("sample", "event"), sep = "_") %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome1 %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
    dplyr::mutate(shortName = sample,
                  class = 'TE') %>% 
    dplyr::mutate(shortName = factor(shortName)) %>%
    dplyr::mutate(gene_hit = as.character(ifelse(!feature1 %in% c('intergenic', 'intron') & gene1 %in% drivers, 'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(!feature1 %in% c('intergenic', 'intron') & tolower(gene1) %in% bed_file$gene, as.character(gene1), 'Other'))) %>% 
    droplevels()
  
  te_data <- convert_names(x=te_data, short_to_long=T)
  te_data <- te_data %>% 
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit)
  
  return(te_data)
}


highest_N <- function(gene_hits){
  new_df <- gene_hits %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(highest_n = ifelse(any(gene_hit == 'Notch'), max(cell_fraction[gene_hit=='Notch']), 0)) %>% 
    dplyr::distinct(sample, highest_n) %>% 
    ungroup() %>% 
    as.data.frame()
  
  return(new_df)
}


convert_names <- function(..., df, attach_info){
  
  all_sample_names <- read.delim(attach_info, header = F)
  colnames(all_sample_names) <- c("sample", "sample_short", "sample_paper", "sex", "assay")
  
  df <- plyr::join(df, all_sample_names, "sample", type = 'left') %>% 
    dplyr::filter(...) %>% 
    dplyr::rename(sample_old = sample,
                  sample = sample_paper) %>% 
    dplyr::select(sample, everything())

  return(df)
}


include <- c("A373R1", "A373R11", "A373R13", "A373R3", "A373R5", "A373R9", 
             "A512R19", "A512R21", "A512R23", "A573R25", "A573R27", "A573R29", 
             "A573R31", "A573R33", "B241R35", "B241R37", "B241R39", "B241R41-1", 
             "B241R43", "B241R45", "B241R49", "B241R51", "B241R53", "B241R55", 
             "B241R57", "B241R59", "B241R61", "B241R63", "D106R23", "D106R25", 
             "D106R27", "D106R29", "D106R31", "D106R33", "HUM-1", "HUM-4", 
             "HUM-7")

plot_tumour_evolution <- function(..., all_samples = 'data/all_samples_merged.txt',
                              all_samples_snvs = 'data/annotated_snvs.txt',
                              all_samples_indels = 'data/annotated_indels.txt',
                              attach_info = 'data/samples_names_conversion.txt',
                              annotate_with = 'data/dna_damage_repair.bed', drivers=c('N', 'kuz'),
                              purity_adj=TRUE, show_sample = TRUE, snvs=TRUE, indels=TRUE, tes=TRUE, plot=TRUE){
  
  bed_file <- readBed(bed_in = annotate_with)
  colnames(bed_file) <- c("chrom", "start", "end", "gene", "id")
  gene_hits <- get_SVs(all_samples=all_samples, bed_file=bed_file, drivers=drivers)
  Notch_svs <- highest_N(gene_hits)
  
  if(snvs) snv_hits <- get_SNVs( all_samples_snvs = all_samples_snvs, bed_file=bed_file, drivers=drivers, purity_adj=purity_adj)
  Notch_snvs <- highest_N(snv_hits)
  
  if(indels) indel_hits <- get_INDELs(all_samples_indels = all_samples_indels, bed_file=bed_file, drivers=drivers, purity_adj=purity_adj)
  Notch_indels <- highest_N(indel_hits)
  
  df <- do.call("rbind", list(indel_hits, snv_hits, gene_hits))
  
  # Get max Notch accross all mut types
  Notch_hits <- merge(Notch_svs,  Notch_snvs,   by = "sample", all = TRUE)
  Notch_hits <- merge(Notch_hits, Notch_indels, by = "sample", all = TRUE)
  colnames(Notch_hits) <- c('sample', 'sv', 'snv', 'indel')
  
  Notch_hits <- Notch_hits %>% 
    dplyr::mutate(highest_n=pmax(sv, snv, indel))
  
  
  Notch_hits[is.na(Notch_hits)] <- 0
  
  df <- convert_names(..., df=df, attach_info=attach_info)
  
  Notch_hits <- convert_names(..., df=Notch_hits, attach_info=attach_info)
  
  df$sample <- factor(df$sample, levels(fct_reorder(Notch_hits$sample, -Notch_hits$highest_n)))
  
  df <- df %>% 
    dplyr::filter(..., !is.na(sample)) %>% 
    dplyr::mutate(time = 1-cell_fraction) %>% 
    dplyr::mutate(class = factor(as.character(class))) %>% 
    droplevels()
  
  Notch_hits <- Notch_hits %>% 
    dplyr::filter(...) %>% 
    dplyr::mutate(time = 1-highest_n) %>% 
    dplyr::mutate(sample=fct_reorder(sample, time)) %>%
    droplevels()
  
  if(!plot) return(list(Notch_hits, df))
  
  special_hits <- df %>% 
    dplyr::filter(special_hit != "Other")
  
  special_hits$sample <- factor(special_hits$sample, levels(fct_reorder(Notch_hits$sample, -Notch_hits$highest_n)))
  
  special_hits %>% 
    dplyr::mutate(time = 1-cell_fraction) %>% 
    dplyr::select(sample, class, special_hit, time) %>% 
    print()
  
  # Custom reorder
  df$class = factor(df$class, levels=c("SV","SNV", "INDEL"), labels=c("SV","SNV", "INDEL")) 
  special_hits$class = factor(special_hits$class, levels=c("SV","SNV", "INDEL"), labels=c("SV","SNV", "INDEL")) 
   
  red <- "#FC4E07"
  blue <- "#00AFBB"
  yellow <- "#F5D658"
  grey <- "grey"
  
  cols <- c(blue, yellow, red)
  
  p <- ggplot(Notch_hits)
  p <- p + geom_violin(data=df, aes(sample, time, fill = class, colour=class), alpha = 0.3, show.legend = FALSE,  adjust=0.3, scale='width')
  p <- p + geom_jitter(data=df, aes(sample, time, fill = class, colour=class), width=0.2, height = 0.01, size=.5, alpha = 0.3, show.legend = FALSE)
  
  p <- p + geom_point(data=special_hits, aes(sample, time), size=2, shape=4, show.legend = FALSE)
  p <- p + geom_errorbar(data=Notch_hits, aes(fct_reorder(sample, time), ymin=time, ymax=time), alpha=0.6, size=0.7, color='black', show.legend = FALSE)
  p <- p + coord_flip()
  p <- p + labs(x = NULL)
  p <- p + scale_fill_manual(values=cols)
  p <- p + scale_colour_manual(values=cols)
  p <- p +  scale_y_continuous('Pseudotime (1 - cell fraction)', labels=seq(0,1,by=.25), expand = c(0, 0.01))
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.x = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.title.x = element_text(size = 15),
      axis.text = element_text(size = 12),
      panel.spacing = unit(2, "lines"))
 
  if(!show_sample) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
  p <- p + facet_wrap(~class, nrow=1)
  print(p)
  
  return(Notch_hits)
}


