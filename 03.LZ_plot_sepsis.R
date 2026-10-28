library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(data.table)
library(vroom)
#library(ggplot2)
#library(biomaRt)

##############sepsis GWAS############################
#####################################################

prepare_locus_with_ld <- function(chr, pos, snp_id, 
                                  gwas_dir, ld_dir,
                                  flank = 5e5, 
                                  ens_db = EnsDb.Hsapiens.v86) {
  
  # Construct file paths
  gwas_file <- file.path("gwasfile")
  ld_file   <- file.path("ldfile")
  bim_file  <- file.path("bimfile")
  
  # Read GWAS
  sepsis_gwas <- read.table(gwas_file, sep = "\t", header = TRUE)
  
  # Parse SNP column
  sepsis_gwas <- sepsis_gwas %>%
    separate(SNP, into = c("chrom", "pos", "a2", "a1"), sep = ":", remove = FALSE) %>%
    mutate(
      chrom = sub("^chr", "", chrom),
      pos = as.integer(pos),
      keep = TRUE
    )
  
  
  # Read LD matrix and BIM file
  ld_df <- read.table(ld_file)
  bim_df <- read.table(bim_file)
  
  sepsis_gwas<- sepsis_gwas %>% dplyr::filter(SNP %in% bim_df$V2)
  
  # Create locus object
  loc <- locus(
    data = sepsis_gwas,
    index_snp = snp_id,
    flank = flank,
    ens_db = ens_db,
    LD = "ld"  # We'll define this column
  )
  
  
  # Find the index SNP
  snp_index <- which(bim_df$V2 == snp_id)
  if (length(snp_index) == 0) {
    stop("Index SNP not found in BIM file")
  }
  
  # Extract LD for that SNP
  ld_for_plot <- ld_df[, snp_index, drop = FALSE]
  ld_for_plot$SNP <- bim_df$V2
  colnames(ld_for_plot)[1] <- "ld"
  
  # Merge LD info
  loc$data <- loc$data %>%
    left_join(ld_for_plot[, c("SNP", "ld")], by = "SNP")
  
  return(loc)
}

snps<-c("chr14:94378610:C:T", "chr19:33096952:C:T", "chr2:632609:A:G", "chr6:32224045:G:A")
chromosomes<-c(14,19,2,6)
positions<-c(94378610, 33096952, 632609, 32224045)


for (i in 1:length(snps)){
  chr <- chromosomes[i]
  pos <- positions[i]
  snp_id <- snps[i]
  
  loc <- prepare_locus_with_ld(
    chr = chr,
    pos = pos,
    snp_id = snp_id,
    gwas_dir = "../gwas_summ_eur/",
    ld_dir = "../locusZoom_plots/sepsis_LD/"
  )
  out_dir<-"../locusZoom_plots/"
  filename=paste0(out_dir,"LZ_chr",chr,"_pos",pos,"_",snp_id,".pdf")
  pdf(filename,
      width = 8, height = 6)
  locus_plot(loc, labels = "index", filter_gene_biotype = "protein_coding", maxrows = 4) 
  dev.off() 
}

specific_r2 <- ld_matrix[var1, var2]

