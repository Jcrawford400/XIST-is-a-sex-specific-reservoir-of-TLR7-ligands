###############################################################################
### First, download the human transcriptome from RefSeq
###############################################################################
# Load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genomes")
library(Biostrings)
BiocManager::install("AnnotationHub")
install.packages("biomartr")
library("biomartr")
install.packages("stringr")
library("stringr")
library(magrittr)
library(dplyr)

library("biomartr")
options(timeout = 30000)
# Get human transcriptome from RefSeq
HS.rna <- getRNA(
  db = "refseq",
  organism = "Homo sapiens",
  reference = FALSE,
  release = NULL,
  path = file.path("_ncbi_downloads", "RNA")
)
Human_rna <- read_cds(file     = HS.rna, 
                      obj.type = "Biostrings" ) 
###############################################################################
### To save time later on, we will remove all but 1 version of genes with 
### multiple "variants" annotated
###############################################################################
all_genes = data.frame("Gene Name" = character(0), "Sequence" = character(0), "UU Count" = numeric(0), "UU Ratio" = numeric(0), "Section" = numeric(0), "Length" = numeric(0))
back_half = unlist(strsplit(Human_rna@ranges@NAMES[[1]], split= "="))[2]
gene_name = unlist(strsplit(back_half, split = "]"))[1]
all_genes[[1,1]] = gene_name
all_genes[[1,2]] = toString(Human_rna[[1]])
index = 1

for (i in 2:174978){
  print(i)
  back_half = unlist(strsplit(Human_rna@ranges@NAMES[[i]], split= "="))[2]
  gene_name = unlist(strsplit(back_half, split = "]"))[1]
  back_half_preceding = unlist(strsplit(Human_rna@ranges@NAMES[[i-1]], split= "="))[2]
  gene_name_preceding = unlist(strsplit(back_half_preceding, split = "]"))[1]
  if (gene_name != gene_name_preceding){
    index = index + 1
    all_genes[[index,1]] = gene_name # Add the gene name to our dataframe
    all_genes[[index,2]] = toString(Human_rna[[i]]) # Add the sequence
  }
}
# Now we have a list of 45,768 mostly unique genes, with 1 transcript each. 
# Now we will sort the names alphabetically and then purify further
all_genes = all_genes[order(all_genes$Gene.Name),]
all_genes_2 =  data.frame("Gene Name" = character(0), "Sequence" = character(0), "UU Count" = numeric(0), "Expression" = numeric(0), "UU Richness" = numeric(0), "Length" = numeric(0))
index = 0
for (i in 2:45768){
  print(i)
  gene = all_genes[[i,1]]
  preceding_gene_name = all_genes[[i-1,1]]
  if (gene != preceding_gene_name){
    index = index + 1
    all_genes_2[[index,1]] = all_genes[[i,1]] # Add the gene name to our dataframe
    all_genes_2[[index,2]] = all_genes[[i,2]] # Add the sequence
    all_genes_2[[index,3]] = all_genes[[i,3]]
    all_genes_2[[index,4]] = all_genes[[i,4]]
    all_genes_2[[index,5]] = all_genes[[i,5]]
    all_genes_2[[index,6]] = all_genes[[i,6]]
  } else if (gene == preceding_gene_name){
    print("Removed")
  }
}
###############################################################################
### I want to count the UU frequency of all genes in the genome and average them
###############################################################################
uu_freq_each_gene = c()
for (i in 1:39000){
  print(i)
  num_uus = 0
  sequence = all_genes_2[[i,2]]
  for (j in 1:nchar(sequence)){
    if (substr(sequence, j, j+1) == "TT"){
      num_uus = num_uus + 1
      all_genes_2[[i,3]] = num_uus
    }
  }
  uu_freq = num_uus/nchar(sequence)
  uu_freq_each_gene = c(uu_freq_each_gene,uu_freq)
  
}
average_uu_all_transcripts = mean(uu_freq_each_gene)

############## Count local uu frequencies in XIST ##################

xist_sequence = all_genes_2[[37992,2]]
local_uu_freqs = c()
last_nuc = nchar(xist_sequence) - 500
for (i in 1:last_nuc){    # nchar(xist_sequence)-500
  local_count = 0
  end = i+499
  for (j in (i:end)){
    if (substr(xist_sequence, j, j+1) == "TT"){
      local_count = local_count + 1
    }
  }
  local_per = local_count
  local_uu_freqs = c(local_uu_freqs, local_per)
}

path <- "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data"
filename <- paste0(path, '/gapdh_local_uu_freqs_motif.xlsx')
write_xlsx(as.data.frame(local_uu_freqs), path = filename)

################### Now count the max local uu richness of every gene ##########################

for (i in 1:39000){ # For every gene...
  print(i)
  uu_freq = 0
  high_score = 0
  for (j in 1:499){
    if (substr(all_genes_2[[i,2]], j, j+1) == "TT"){ # Then first, count the UUs in the first 500 nucleotides
      uu_freq = uu_freq + 1
    }
  }
  high_score = uu_freq
  for (k in 1:nchar(all_genes_2[[1,2]]) - 500){ # Then, for every nucleotide up to the length of the gene minus 500
    if (substr(all_genes_2[[1,2]], k, k+1) == "TT"){ # If that nucleotide was the start of a UU, subtract 1 from the UU freq
      uu_freq = uu_freq - 1
    }
    if (substr(all_genes_2[[1,2]], k+499, k+500) == "TT"){ # If the 500th nucleotide in this shifting frame is the start of a UU, then add 1 to the UU freq
      uu_freq = uu_freq + 1
    }
    if (uu_freq > high_score){ # Every time you move over 1 nucleotide, if the count went up above its previous high, save that as the new high score
      high_score = uu_freq
    }
  }
  all_genes_2[[i,5]] = high_score # Once you've checked every window of 500 nucleotides, save the high score in the dataframe
}

# Find genes withh the guccuucaa motif #
genes_containing_motif = data.frame(name=character(0), expression=numeric(0))
gtccttcaa_counter = 0
# Full-length is 120689
for (x in 1:39000){
  print(x)
  if (all_genes_2[[x,1]] %in% expression_levels_3$name){
    if (grepl("GTCCTTCAA",all_genes_2[[x,2]], fixed = TRUE)){
      gtccttcaa_counter = gtccttcaa_counter + 1
      genes_containing_motif[[gtccttcaa_counter,1]] = all_genes_2[[x,1]] # Gene Name
      for (z in 1:15860){
        if (expression_levels_2[[z,1]] == all_genes_2[[x,1]]){ # When you find the gene
          genes_containing_motif[[gtccttcaa_counter,2]] = expression_levels_2[[z,3]]
          break
        }
      }
    }
  }
}