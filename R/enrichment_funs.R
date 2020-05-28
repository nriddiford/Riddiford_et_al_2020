library(ggrepel)
hits_in_genes <- function(..., df, bedDir, out_file){
  dataType='svBreaks'
  if("pos" %in% colnames(df)) dataType='mutationProfiles'
  
  non_notch_genes <- df %>% 
    dplyr::filter(..., gene != 'intergenic')
  # all_genes <- bpRegionEnrichment(bp_data = non_notch_genes, dataType=dataType, bedDir = paste(bedDir, '/umr_gene_features/paper/', sep=''), plot=F)
  
  # all_genes$group <- factor("All genes")
  # comb2 <- all_genes
  
  just_genes <- non_notch_genes %>% 
    dplyr::distinct(gene, .keep_all=TRUE)
  
  top_20_percent <- as.integer(0.2*nrow(just_genes))
  bottom_20_percent <- as.integer(0.8*nrow(just_genes))
  
  non_notch_ranked <- non_notch_genes %>% 
    dplyr::arrange(-fpkm) %>% 
    dplyr::mutate(exp_rank = rank(desc(fpkm)))
  # dplyr::filter(rank(desc(fpkm)) <= top_10_percent) %>% View()
  
  exp_genes <- dplyr::filter(non_notch_genes, fpkm > 1)
  
  bps_in_exp_genes <- bpRegionEnrichment(bp_data = exp_genes, dataType=dataType, bedDir = bedDir, plot=F)
  bps_in_exp_genes$group <- factor('All expressed genes')
  comb2 <- bps_in_exp_genes
  
  # high_exp_genes <- dplyr::filter(non_notch_ranked, exp_rank < top_20_percent)
  # if(nrow(high_exp_genes)>0) {
  #   bps_in_high_exp_genes <- bpRegionEnrichment(bp_data = high_exp_genes, dataType=dataType, bedDir = bedDir, plot=F)
  #   bps_in_high_exp_genes$group <- factor('Top 20% FPKM')
  #   comb2 <- rbind(comb2, bps_in_high_exp_genes)
  # }
  
  # low_exp_genes <- dplyr::filter(non_notch_ranked, exp_rank >= bottom_20_percent & fpkm > 1)
  # if(nrow(low_exp_genes)>0) {
  #   bps_in_low_exp_genes <- bpRegionEnrichment(bp_data = low_exp_genes, dataType=dataType, bedDir = bedDir, plot=F)
  #   bps_in_low_exp_genes$group <- factor('Bottom 20% FPKM')
  #   comb2 <- rbind(comb2, bps_in_low_exp_genes)
  # }
  
  non_exp_genes <- dplyr::filter(non_notch_genes, fpkm <= 1)
  bps_in_non_exp_genes <- bpRegionEnrichment(bp_data = non_exp_genes, dataType=dataType, bedDir = bedDir, plot=F)
  
  bps_in_non_exp_genes$group <- factor('Non-expressed genes')
  
  comb2 <- rbind(comb2, bps_in_non_exp_genes)
  
  # comb2 <- rbind(all_genes, bps_in_high_exp_genes, bps_in_low_exp_genes, bps_in_non_exp_genes)
  
  comb2 <- comb2 %>% 
    tidyr::complete(group, feature) %>% 
    dplyr::mutate(Log2FC = ifelse(is.na(Log2FC), 0, Log2FC),
                  label = ifelse(sig %in% c(NA,'-'), '', sig)) %>% 
    dplyr::select(group:eScore, label) %>% 
    as.data.frame()
  
  gf_exp_grouped <- ggbarplot(comb2, x = "feature", y = "Log2FC",
                              fill = "group", color = "group", group = "group", palette =c(blue, yellow, red, grey), alpha = 0.6, 
                              label = 'label', lab.vjust = -.1,
                              orientation = "horiz", x.text.angle = 90,
                              xlab = FALSE, ylab = "Log2(FC)")
  
  gf_exp_grouped <- ggpar(gf_exp_grouped, ylim = c(-2.2,2.2))
  p1 <- facet(gf_exp_grouped, facet.by = 'group', nrow = 1)
  
  df1 <- comb2 %>% 
    dplyr::select(feature:padj, group) %>% 
    dplyr::filter(!is.na(test)) %>% 
    write_tsv(., path = paste0(rootDir, 'Desktop/final_analysis/', out_file))
  
  return(list(p1, df1))
}


run_all <- function(..., bps, snvs, indels, dir, plot=F, slop=0){
  bps_table <- bpRegionEnrichment(..., bp_data = bps, dataType='svBreaks', bedDir = dir, plot=F, slop=slop) %>% 
    dplyr::mutate(group = 'sv')
  snvs_table <- bpRegionEnrichment(..., bp_data = snvs, dataType='mutationProfiles', bedDir = dir, plot=F, slop=slop) %>% 
    dplyr::mutate(group = 'snv')
  indels_table <- bpRegionEnrichment(..., bp_data = indels, dataType='mutationProfiles', bedDir = dir, plot=F, slop=slop) %>% 
    dplyr::mutate(group = 'indel')
  
  
  if(!plot) return(list(bps_table, snvs_table, indels_table))
  
  bp_df <- bps_table %>% 
    dplyr::select(feature, Log2FC) %>% dplyr::rename(sv = Log2FC)
  
  snv_df <- snvs_table %>% 
    dplyr::select(feature, Log2FC) %>% dplyr::rename(snv = Log2FC)
  
  indel_df <- indels_table %>% 
    dplyr::select(feature, Log2FC) %>% dplyr::rename(indel = Log2FC)
  
  comb <- plyr::join(bp_df, snv_df, indel_df, by='feature', type='left')
  
  library(pheatmap)
  
  g <- comb$feature
  comb$feature <- NULL
  
  m <- comb %>% as.matrix()
  
  rownames(m) <- g
  
  # colours <- c("#f1f9fe", "#d4ecfb", "#9bd3f6", "#7ec6f3", "#53b3ef")
  
  colours <- rev(heat.colors(5, 0.5))
  
  pheatmap::pheatmap(m,
                     color = colours,
                     cluster_row = T, cluster_cols = T,
                     treeheight_row = 0, treeheight_col = 0,
                     angle_col = 90
  )
  
}


plot_all <- function(l, escore_threshold=5, maxp=50){
  # devtools::install_github("stefanedwards/lemon", ref='v0.3.3')
  library(lemon)
  
  df <- plyr::join_all(d, type='full')
  
  cols <- c("enriched" = blue,
    "depleted" = red,
    "ns" = grey)
  
  
  maxLog2 <- max(abs(df$Log2FC))
  maxLog2 <- round_any(maxLog2, 1, ceiling)
  
  df$colour <- ifelse(df$eScore >= escore_threshold, 
                      ifelse(df$test=='enrichment', 'enriched', 'depleted'),
                      'ns')
  
  df$label <- ifelse(df$eScore >= escore_threshold, df$feature, '')
  # df$colour = factor(df$colour, levels=c("yes","no"), labels=c("***","ns")) 
  
 
  # cols = c(blue, red, grey)
  if(!nrow(df[df$eScore >= escore_threshold,])) cols = c(grey)
  
  df$group = factor(df$group, levels=c("sv","snv", "indel"), labels=c("SV","SNV", "INDEL")) 
  # df$colour = factor(df$colour, levels=c("enriched","depleted", "ns"), labels=c("enriched","depleted", "ns")) 
  
  
  df$maxp <- ifelse(-log10(df$padj)>maxp, maxp, -log10(df$padj))
  
  p <- ggplot(df, aes(Log2FC, maxp))
  p <- p + geom_point(aes(colour = colour), size=3, alpha=0.7)
  p <- p + scale_fill_manual(values = cols)
  p <- p + scale_colour_manual(values = cols)
  p <- p + geom_text_repel(aes(label = label), size=3, inherit.aes = TRUE)
  p <- p + guides(size = FALSE, colour = FALSE) 
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.title =element_text(size=rel(1.2))
    )
  p <- p + scale_y_continuous("-Log10(padj)")
  p <- p + scale_x_continuous("Log2(FC)", limits=c(-maxLog2, maxLog2), breaks=seq(-maxLog2, maxLog2, by=1))
  p <- p + facet_rep_wrap(~group, ncol=3, repeat.tick.labels = T)
  p
}


run_reldist <- function(a, b, verbose = FALSE, print=TRUE, chroms = chromosomes, type='bed', sim=FALSE, merge_regions=FALSE){
  if(!is.null(nrow(a))) {
    a_bed <- getBed(f = a, verbose = verbose, type=type, chr %in% chroms) %>% 
      droplevels()
    rownames(a_bed) <- paste0(a_bed$chr, ":", a_bed$start, "-", a_bed$end)
  } else {
    a_bed <- getBed(f = a, verbose = verbose)
  }
  
  # a_bed <- bedr.merge.region(a_bed, check.chr = FALSE, verbose = verbose)
  if(!is.null(nrow(b))) {
    b_bed <- getBed(f = b, verbose = verbose, type=type, chr %in% chroms) %>% 
      droplevels()
    rownames(b_bed) <- paste0(b_bed$chr, ":", b_bed$start, "-", b_bed$end)
  } else {
    b_bed <- getBed(f = b, verbose = verbose, chr %in% chroms)
    if(merge_regions) b_bed <- bedr.merge.region(b_bed, check.chr = FALSE, verbose = verbose)
  }

  cat("Running reldist...\n")
  
  dist <- reldist(rownames(a_bed), rownames(b_bed), check.chr = FALSE, verbose = verbose);
  dist$group <- 'real'
  if(print) print(plot_reldist(dist))
  if(sim){
    cat("Simulating data for comparison\n")
    s <- svBreaks::bpSim(nSites = nrow(a_bed)*10)
    s_bed <- getBed(f = s, verbose = F, type='snv') %>% 
      dplyr::filter(chr %in% chroms) %>% 
      droplevels()
    rownames(s_bed) <- paste0(s_bed$chr, ":", s_bed$start, "-", s_bed$end)
    simd <- reldist(rownames(s_bed), rownames(b_bed), check.chr = FALSE, verbose = FALSE);
    simd$group <- 'simulated'
    c <- plyr::join(dist, simd, type='full')
    return(c)
  }
  return(dist)
}



getBed <- function(..., f, verbose=FALSE, type='bed'){
  if(type=='sv'){
    cat("Reading SVs from dataframe\n")
    b <- f %>% 
      dplyr::mutate(chr = as.character(chrom),
                    start = as.integer(bp),
                    end = as.integer(bp) + 1) %>% 
      dplyr::distinct(chr, start, end, .keep_all=TRUE) %>% 
      dplyr::select(chr, start, end) %>% 
      droplevels() %>% 
      as.data.frame()
    
  } else if(type=='snv'){
    cat("Reading SNVs from dataframe\n")
    b <- f %>% 
      dplyr::mutate(chr = as.character(chrom),
                    start = as.integer(pos),
                    end = start + 1) %>% 
      dplyr::distinct(chr, start, end, .keep_all=TRUE) %>% 
      dplyr::select(chr, start, end) %>% 
      droplevels() %>% 
      as.data.frame()
    
  } else{
    b <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
    colnames(b) <- c('chr', 'start', 'end')
  }
  b <- b %>% 
    dplyr::filter(...) %>% 
    dplyr::select(chr, start, end) %>% 
    droplevels() %>% 
    as.data.frame()
  
  b_sort <- bedr.sort.region(b, check.chr = FALSE, verbose = verbose)
  b_bedr <- bedr(engine = "bedtools", input = list(i = b_sort), method = "sort", check.chr = FALSE, verbose = verbose)
  return(b_bedr)
}


plot_reldist <- function(d){
  
  blue <- '#259FBF'
  cols <- c(blue, grey)
  
  p <- ggplot(d, aes(reldist, fraction, colour=group), alpha = .7)
  p <- p + geom_point()
  p <- p + geom_smooth(aes(fill=group), alpha=0.3, span=0.4)
  p <- p + scale_fill_manual(values = cols)
  p <- p + scale_colour_manual(values = cols)
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.title =element_text(size=rel(1.2))
    )
  p <- p + guides(fill=FALSE, colour=FALSE)
  print(p)
}
