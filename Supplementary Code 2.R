###################### Jonathan Crawford #######################################
### XIST Paper Supplementary Code

### Figure 1A: Look for sex-biased genes that contain known TLR7 Ligand GUCCUUCAA
### First, download the human transcriptome

# Load necessary packages
install.packages("stringr")
library("stringr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genomes")
library(Biostrings)
BiocManager::install("AnnotationHub")
install.packages("biomartr")
library("biomartr")

library(magrittr)


# Download transcriptome with getRNA
HS.rna <- getRNA(
  db = "refseq",
  organism = "Homo sapiens",
  reference = FALSE,
  release = NULL,
  path = file.path("_ncbi_downloads", "RNA")
)
Human_rna <- read_cds(file     = HS.rna, 
                      obj.type = "Biostrings" ) 
# Now clean the transcriptome

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
# Now we will sort the names alphabetically and then purify further. 
all_genes = all_genes[order(all_genes$Gene.Name),]
all_genes_2 =  data.frame("Gene Name" = character(0), "Sequence" = character(0), "UU Count" = numeric(0))
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
  } else if (gene == preceding_gene_name){
    print("duplicate")
  }
}
# all genes 2 contains every unique transcript and its sequence

# Now load the list of sex biased genes from avg_sex_bias_mele (Derived from Mele et al)
# new_data_set = read.table(file = "clipboard", 
                          # sep = "\t", header=TRUE)
library(readxl)
path = "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data/"
new_data_set = read_xlsx(paste0(path, "avg_sex_bias_mele.xlsx"))

names(new_data_set)[1:(ncol(new_data_set)-1)] = names(new_data_set)[2:(ncol(new_data_set))]
sex_spec_genes_containing_motif = data.frame(name=character(0), sex_bias=numeric(0), neg_log_q = numeric(0))
gtccttcaa_counter = 0
# Full-length is 120689
for (x in 1:39000){
  if (all_genes_2[[x,1]] %in% new_data_set$gene_name){
    if (grepl("GTCCTTCAA",all_genes_2[[x,2]], fixed = TRUE)){
      gtccttcaa_counter = gtccttcaa_counter + 1
      sex_spec_genes_containing_motif[[gtccttcaa_counter,1]] = all_genes_2[[x,1]] # Gene Name
      for (z in 1:826){
        if (new_data_set[[z,1]] == all_genes_2[[x,1]]){ # When you find the gene
          print(gtccttcaa_counter)
          sex_spec_genes_containing_motif[[gtccttcaa_counter,2]] = new_data_set[[z,7]]
          sex_spec_genes_containing_motif[[gtccttcaa_counter,3]] = -log10(new_data_set[[z,9]]) # -log of adjusted p value
          break
        }
      }
    }
  }
}

# Start Figure 1B -- counting UUs in each sex biased gene
# Now load just the gene names from the list of sex biased genes from avg_sex_bias_mele (Derived from Mele et al)
# list_of_sex_biased_genes = read.table(file = "clipboard", sep = "\t", header=TRUE)
list_of_sex_biased_genes = read_xlsx("C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data/list_sex_biased_genes.xlsx")
list_of_sex_biased_genes = read_xlsx(paste0(path, "list_sex_biased_genes.xlsx"))
list_of_sex_biased_genes$index = 0
list_of_sex_biased_genes$uu_rich = 0
for (i in 1:826){
  print(i)
  gene_name = list_of_sex_biased_genes[[i,1]]
  for (j in 1:39000){
    if (list_of_sex_biased_genes[[i,1]] == all_genes_2[[j,1]]){
      print("found!")
      list_of_sex_biased_genes[[i,2]] = j
      
    }
  }
}

list_indexed_sb_genes = data.frame(gene_name=character(0), index = numeric(0), log2FoldChange=numeric(0), uu_rich = numeric(0), uu_total = numeric(0))
counter=0
for (i in 1:826){
  if (list_of_sex_biased_genes[[i,2]] != 0){
    counter=counter+1
    list_indexed_sb_genes[[counter,1]] = list_of_sex_biased_genes[[i,1]]
    list_indexed_sb_genes[[counter,2]] = list_of_sex_biased_genes[[i,2]]
    list_indexed_sb_genes[[counter,3]] = 0
    list_indexed_sb_genes[[counter,4]] = 0
    list_indexed_sb_genes[[counter,5]] = 0
  }
}
for (i in 1:586){
  print(i)
  for (j in 1:39000){
    if (list_indexed_sb_genes[[i,1]] == all_genes_2[[j,1]]){
      num_uus = 0
      sequence = all_genes_2[[j,2]]
      for (j in 1:nchar(sequence)){
        if (substr(sequence, j, j+1) == "TT"){
          num_uus = num_uus + 1
        }
      }
      list_indexed_sb_genes[[i,5]] = num_uus
    }
  }
}
sex_biased_genes_with_fold_change = read_xlsx(paste0(path, "avg_sex_bias_3_columns.xlsx"))
for (i in 1:586){
  print(i)
  for (j in 1:826){
    if (list_indexed_sb_genes[[i,1]] == sex_biased_genes_with_fold_change[[j,1]]){
      list_indexed_sb_genes[[i,3]] = sex_biased_genes_with_fold_change[[j,3]]
    }
  }
}
# Now we have the data for Figure 2B. To get the data for figure 2C, need to find max UU Richness of every gene

for (i in 1:586){
  print(i)
  for (j in 1:39000){
    if (list_indexed_sb_genes[[i,1]] == all_genes_2[[j,1]]){
      max_uu_richness = 0
      sequence = all_genes_2[[j,2]]
      end_nuc=nchar(sequence) - 500
      for (k in 1:end_nuc){
        uu_count = 0
        fir=as.numeric(k)
        las=fir+499
        for (l in 1:500){
          if (substr(all_genes_2[[j,2]],fir+l,fir+l+1) == "TT"){
            uu_count = uu_count + 1
          }
        }
        if (uu_count > max_uu_richness){
          max_uu_richness=uu_count
        }
      }
      list_indexed_sb_genes[[i,4]] = max_uu_richness
    }
  }
}
# Now we have the data for Figure 1C

library("writexl")
filename <- paste0(path, '/list_indexed_sb_genes_4_24_22.xlsx')
write_xlsx(list_indexed_sb_genes, path = filename)

### Figure 1D: UU Richness of XIST across the sequence
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
filename <- paste0(path, '/local_uu_freqs_4_24_22.xlsx')
write_xlsx(as.data.frame(local_uu_freqs), path = filename)



###############################################################################
### End Figure 1
###############################################################################
### Figure S2: XIST RNA Sequencing 
library(readxl)
#xist_kd_all_data <- read_excel(file.choose()) # XIST KD RNA Sequencing Correct
path="C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data/"
xist_kd_all_data = read_xlsx(paste0(path, "XIST KD RNA Seq Correct.xlsx"))
gene_names_xistkd = xist_kd_all_data$`All Gene Symbol original`
small_df = data.frame(gene = character(0), p_val = numeric(0), adj_p = numeric(0), log2_fold_change = numeric(0), Up_Down = character(0), label = character(0), wt_exp = numeric(0), x2kd_exp = numeric(0))
for (i in 1:11678){
  print(i)
  small_df[i,1] = xist_kd_all_data[i,1]
  small_df[i,2] = xist_kd_all_data[i,18]
  small_df[i,3] = xist_kd_all_data[i,21]
  small_df[i,4] = xist_kd_all_data[i,27]
  small_df[i,5] = 0
  small_df[i,6] = 0
  small_df[i,7] = xist_kd_all_data[i,10]
  small_df[i,8] = xist_kd_all_data[i,11]
}
library(ggplot2)
p <- ggplot(data=small_df, aes(x=log2_fold_change, y=-log10(adj_p))) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2
small_df$diffexpressed <- "NO"
small_df$diffexpressed[small_df$log2_fold_change > 1 & small_df$adj_p < 0.05] <- "UP"
small_df$diffexpressed[small_df$log2_fold_change < -1 & small_df$adj_p < 0.05] <- "DOWN"
p <- ggplot(data=small_df, aes(x=log2_fold_change, y=-log10(small_df$adj_p), col=diffexpressed)) + geom_point() + theme_minimal()
p
p2 <- p + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
p3
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("UP", "DOWN", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

### Add labels ###
small_df$label <- NA
for (i in 1:11678){
  print(i)
  if (small_df[i,9] == "UP" & small_df[i,4] > 3){
    small_df[i,6] = small_df[i,1]
  }
  if (small_df[i,9] == "DOWN" & small_df[i,4] < -3){
    small_df[i,6] = small_df[i,1]
  }
}
library(ggrepel)

p2 <- ggplot(data=small_df, aes(x=log2_fold_change, y=-log10(adj_p), col=diffexpressed, grid.col=NA, label=label )) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()
mycolors <- c("#B3364C", "#5D6AB0", "black")
names(mycolors) <- c("UP", "DOWN", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3
p4 <- p3 + geom_vline(xintercept=c(-1, 1), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed")
p4
mytheme <- theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

print(p4 + mytheme)

###############################################################################
### Make the table for Figure S2D
# First, we need to count the number of UUs in each gene
for (i in 1:39000){
  print(i)
  uu_count = 0
  sequence = all_genes_2[[i,2]]
  length = nchar(all_genes_2[[i,2]]) - 1
  for (j in 1:length){
    if (substr(all_genes_2[[i,2]], j, j+1) == "TT"){
      uu_count = uu_count + 1
    }
  }
  all_genes_2[[i,3]] = uu_count
}



# Now, let's check out the UU richness of each gene that was lost in the knockout
uu_and_fc = data.frame(gene = character(0), expression = numeric(0), uu = numeric(0), guccuucaa = numeric(0))
counter = 0
for (i in 1:11678){
  print(i)
  if (small_df[i,9] == "DOWN" & small_df[i,4] < -2){
    name = small_df[i,1]
    counter = counter + 1
    for (j in 1:39000){
      if (all_genes_2[j,1] == name){
        uu_and_fc[counter,1] = name
        uu_and_fc[counter,2] = small_df[i,7]
        uu_and_fc[counter,3] = all_genes_2[j,3]
        sequence = all_genes_2[[j,2]]
        if (grepl("GTCCTTCAA",sequence, fixed = TRUE)){
          uu_and_fc[[counter,4]] = 1
        }
        #uu_and_fc[counter,4] = all_genes_2[j,5]
        
      }
    }
  }
  
}
filename <- paste0(path, '/uu_and_fc.xlsx')
write_xlsx(uu_and_fc, path = filename)







