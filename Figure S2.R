################################################################################
### Figure S2
################################################################################

list_of_ko_lost_genes = read.table(file = "clipboard", sep = "\t", header=FALSE)
list_of_ko_lost_genes$index = 0
list_of_ko_lost_genes$uu_rich = 0
for (i in 1:587){
  print(i)
  gene_name = list_of_ko_lost_genes[[i,1]]
  for (j in 1:39000){
    if (list_of_ko_lost_genes[[i,1]] == all_genes_2[[j,1]]){
      print("found!")
      list_of_ko_lost_genes[[i,2]] = j
      
    }
  }
}

# Now we want to calculate the max UU richness of every gene on this list

for (i in 1:65){
  print(i)
  index = list_of_ko_lost_genes[[i,2]]
  high_score = 0
  uu_count = 0
  sequence = all_genes_2[[index,2]] # Pull out the sequence
  seq_length = nchar(sequence)
  sto = seq_length - 500
  for (j in 1:sto){
    uu_count = 0
    for (k in 1:500){
      if (substr(sequence, k+j, k+j+1) == "TT"){
        uu_count = uu_count + 1
      }
    }
    if (uu_count > high_score){
      high_score = uu_count
      list_of_ko_lost_genes[[i,3]] = high_score
    }
  }
  
}

install.packages("writexl")
library("writexl")
# write_xlsx(as.data.frame(list_of_ko_lost_genes))

path <- "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data"
filename <- paste0(path, '/local_uu_freqs.xlsx')
write_xlsx(as.data.frame(list_of_ko_lost_genes), path = filename)