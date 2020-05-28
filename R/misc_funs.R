addPurity <- function(df, purity_file='/Users/Nick_curie/Desktop/script_test/svSupport/data/tumour_purity.txt'){
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
