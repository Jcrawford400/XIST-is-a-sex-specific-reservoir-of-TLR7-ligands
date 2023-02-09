### 
# XIST Paper Code for Github
###

# Figure 1:  
# First, Load the curated list of sex biased genes from Mele et. al. Then calculate UU-richness
# Load the curated list of sex biased genes
list_of_sex_biased_genes = read.table(file = "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Github/avg_sex_bias_3_columns.csv", sep = "\t", header=TRUE) # This is the 3-column version I made. Curated list of all sex biased genes in Mele dataset
list_of_sex_biased_genes$index = 0
list_of_sex_biased_genes$uu_count = 0
list_of_sex_biased_genes$uu_rich = 0
list_of_sex_biased_genes$guccuucaa = 0
# Calculate the UU dinucleotide content of each gene in the list
# Download the human transcriptome
library(Biostrings)
library("biomartr")
library("stringr")
library(magrittr)

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
# Now we have a list of 45,768 MOSTLY UNIQUE genes, with 1 transcript each. 
# Now we will sort the names alphabetically and then purify.
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
    print("Busted")
  }
}

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

####################################################################
### Count the Max UU Richness of all genes in the transcriptome ###
####################################################################
for (i in 1:39000){
  print(i)
  high_score = 0
  uu_count = 0
  sequence = all_genes_2[[i,2]] # Pull out the sequence
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
      all_genes_2[[i,5]] = high_score
    }
  }
  
}

for (i in 1:826){
  print(i)
  gene_name = list_of_sex_biased_genes[[i,1]]
  for (j in 1:39000){
    if (list_of_sex_biased_genes[[i,1]] == all_genes_2[[j,1]]){
      print("found!")
      list_of_sex_biased_genes[[i,4]] = j
      list_of_sex_biased_genes[[i,5]] = all_genes_2[[j,3]]
      list_of_sex_biased_genes[[i,6]] = all_genes_2[[j,5]]
      list_of_sex_biased_genes[[i,7]] = grepl("GTCCTTCAA",all_genes_2[[j,2]], fixed = TRUE)
      if (is.na(list_of_sex_biased_genes[[i,6]])){
        list_of_sex_biased_genes[[i,6]] = 0
      }
    }
  }
}
# Figure 1A-C come from list_of_sex_biased_genes
# For Fig 1D, we need to use the expression matrix
library(CePa)
dat.gct = read.gct(file = "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Github/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm (1).gct.gz") # Read median transcripts per million off of GTEx site
library(EnsDb.Hsapiens.v79)


dat.gct.2.df = as.data.frame(dat.gct)
dat = dat.gct.2.df

dat[,55] = rowMeans(dat)



rownames(dat) = sub('\\.[0-9]*$', '', rownames(dat))

# First, Convert Ensembl gene IDs in dat.gct to my names
for (i in 1:56200){
  print(i)
  tryCatch({
    rownames(dat)[i] = ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(dat)[i], keytype = "GENEID", columns = c("SYMBOL","GENEID"))[1]
  }, error=function(e){cat("ERROR", conditionMessage(e), "\n")})
}


for (i in 1:56200){
  print(i)
  fexp = mean(dat[i,41], dat[i,31], dat[i,24], dat[i,52], dat[i,53])
  mexp = mean(dat[i,50], dat[i,44])
  if (fexp > 5*dat[i,55]){
    dat[i,55] = fexp
  } else if (mexp > 5*dat[i,55]){
    dat[i,55] = mexp
  }
}


for (i in 1:39000){
  print(i)
  gene_name = all_genes_2[[i,1]]
  for (j in 1:56200){
    if (gene_name == rownames(dat)[j]){
      all_genes_2[i,4] = dat[j,55] # Gotta check this
      break
    }
  }
}


# Expression and UU data comes from all_genes_2

library("writexl")
write_xlsx(as.data.frame(table_1b))

path <- "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Github"
filename <- paste0(path, '/table_1b.xlsx')
write_xlsx(as.data.frame(table_1b), path = filename)

# Now Make Fig 1H, the UU richness of XIST across the transcript
xist_sequence = all_genes_2[[37992,2]]
heatmap = c()
uu_freq = 0
for (j in 1:499){
  if (substr(xist_sequence, j, j+1) == "TT"){ # Then first, count the UUs in the first 500 nucleotides
    uu_freq = uu_freq + 1
  }
}
heatmap = c(heatmap, uu_freq)
for (k in 1:18780){
  if (substr(xist_sequence, k, k+1) == "TT"){ # If that nucleotide was the start of a UU, subtract 1 from the UU freq
    uu_freq = uu_freq - 1
  }
  if (substr(xist_sequence, k+499, k+500) == "TT"){ # If the 500th nucleotide in this shifting frame is the start of a UU, then add 1 to the UU freq
    uu_freq = uu_freq + 1
  }
  heatmap = c(heatmap,uu_freq)
}

################################################################################
### Make Figure 5A 
################################################################################

# Load in AMP expression matrix and format it for seurat

all_amp_data = read.table(file = "clipboard", 
                          sep = "\t", header=FALSE, stringsAsFactors = FALSE) 

cell_barcodes = all_amp_data[1,]
cell_barcodes=cell_barcodes[-1]
gene_names = all_amp_data[,1]
gene_names=gene_names[-1]
amp_mat = all_amp_data[,-1]
amp_mat = amp_mat[-1,]

rows=rownames(amp_mat)
cols=colnames(amp_mat)
dimnames = list(rows,cols)

amp_mat = matrix(as.numeric(unlist(amp_mat)), ncol = ncol(amp_mat), dimnames=dimnames)


colnames(amp_mat) = cell_barcodes
rownames(amp_mat) = gene_names

library(dplyr)
library(Seurat)
library(patchwork)

# Convert to Seurat object
amp_meta = as.data.frame(cell_barcodes)
amp_meta = sub("_", "", amp_meta)
amp_meta = sub("_", "", amp_meta) # Have to run it twice, once to remove each _
amp_meta = as.data.frame(amp_meta)
rownames(amp_meta) = cell_barcodes
# Make Seurat object of AMP data
amp_seu = CreateSeuratObject(counts = amp_mat, meta.data = amp_meta)
amp_seu = NormalizeData(amp_seu, normalization.method = "LogNormalize", scale.factor=10000)
amp_seu <- FindVariableFeatures(amp_seu, selection.method = "vst", nfeatures = 22709)
amp_seu <- ScaleData(amp_seu)
all.genes=row.names(amp_seu)

# Get some basic params
amp_seu[["percent.mt"]] <- PercentageFeatureSet(amp_seu, pattern = "^MT-")
amp_seu <- RunPCA(amp_seu, features = VariableFeatures(object = amp_seu))

# Cluster AMP cells so we can see what we have
amp_seu <- FindNeighbors(amp_seu, dims = 1:10)
amp_seu <- FindClusters(amp_seu, resolution = 0.5)
head(Idents(amp_seu), 5)
amp_seu = RunUMAP(amp_seu, dims=1:21)
DimPlot(amp_seu, reduction = "umap")
cluster0.markers <- FindMarkers(amp_seu, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(amp_seu, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(amp_seu, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(amp_seu, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(amp_seu, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(amp_seu, ident.1 = 5, min.pct = 0.25)
cluster6.markers <- FindMarkers(amp_seu, ident.1 = 6, min.pct = 0.25)
cluster7.markers <- FindMarkers(amp_seu, ident.1 = 7, min.pct = 0.25)
cluster8.markers <- FindMarkers(amp_seu, ident.1 = 8, min.pct = 0.25)
cluster9.markers <- FindMarkers(amp_seu, ident.1 = 9, min.pct = 0.25)

VlnPlot(amp_seu, features = c("MS4A1"))
FeaturePlot(amp_seu, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                                  "CD8A"))
FeaturePlot(amp_seu, features = c("XIST"), slot="scale.data")


View(amp_seu)

# Add the IFN score to the dataset. Use the same IFN score as AMP

ifn = c("ABCA1", "ACTA1", "ACTA2", "ADAR", "AIM2", "APOL6", "ASPRV1", "BATF2", "BST2", "BTN3A1", "CARD16", "CARD17", "CCL8", "CEACAM1", "CHMP5", "CMPK2", "CPT1B", "CXCL10", "DDX58", "DDX60", "DHRS9", "DHX58", "DYNLT1", "EIF2AK2", "FBX06", "GADD45B", "GALM", "GBP1", "GBP2", "GBP3", "GBP4", "GBP6", "HELZ2", "HERC5", "HERC6", "HES4", "HSH2D", "IDO1", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM3", "IRF7", "IRF9", "ISG20", "LAP3", "LGALS9", "LY6E", "MOV10", "MT1A", "MX1", "NBN", "NCOA7", "NMI", "NT5C3A", "OAS1", "OAS2", "OAS3", "OASL", "OTOF", "PARP10", "PARP12", "PARP14", "PARP9", "PHF11", "PSMB9", "RBCK1", "REC8", "RHBDF2", "RNF213", "RSAD2", "SERPING1", "SOCS1", "SP100", "SP110", "SP140", "SPATS2L", "SRBD1", "STAT1", "STAT2", "TAP1", "TAP2", "TCN2", "TIMM10", "TMEM140", "TNFAIP6", "TNFSF10", "TRAFD1", "TRANK1", "TRIM21", "TRIM22", "TRIM38", "TRIM5", "TRIM56", "TRIM6", "TYMP", "UBE2L6", "UNC93B1", "WARS", "XAF1", "ZBP1", "ZC3HAV1", "ZNF684", "ZNFX1")
ifn_score = AddModuleScore(amp_seu, features = list(ifn), name="IFN.Module")

## Get a FeaturePlot of the IFN Score

ifn_score_2 = ifn_score
ifn_mod_list = ifn_score_2@meta.data[["IFN.Module1"]]
ifn_score_2@assays$RNA@scale.data[1,] = ifn_mod_list
ifn_score_2@assays$RNA@data@Dimnames[[1]][1] = 'A1BG'
FeaturePlot(ifn_score_2, features = c('A1BG'), slot="scale.data")

for (i in 1:2838){
  print(i)
  ifn_score_2@assays$RNA@scale.data[1,i] = ifn_score_2@meta.data[["IFN.Module1"]][i]
}

VlnPlot(ifn_score_2, features = c("A1BG"), slot="scale.data", group.by="seurat_clusters")

# Now we have to account for which cells belong to which patient so we can analyze XIST expression vs IFN signature
# Load in the key with disease status
key_with_disease_status = read.table(file.choose()) # celseq.meta file
key_with_disease_status_2 = key_with_disease_status # Make a copy
cell_clusters = read.table(file = "clipboard", 
                           sep = "\t", header=FALSE, stringsAsFactors = FALSE) # This is the cluster_per_cell file
for (i in 1:8297){
  print(i)
  cell_name = key_with_disease_status_2[i,1]
  for (j in 1:2838){
    if (cell_name == cell_clusters[j,1]){
      key_with_disease_status_2[i,6] = cell_clusters[j,2] # Add in what cell type each cell is, to potentially explore this later
      break
    }
  }
}

# Make a metadata variable with patient ID
pt_metadata = data.frame(cell_name=character(0), patient=character(0))
disease_meta = data.frame(cell_name = character(0), status=character(0))
cell_type = data.frame(cell_name = character(0), cluster=character(0))

# We don't want males or patients with fewer than 10 cells in our dataset, so make 2 lists to exclude them. Disease key has healthy control status, so we don't need to make that list
males = c("200-0341", "200-0115", "200-0339", "200-0610", "200-0615", "200-0871") 
exclude = c("200-0107", "200-0115", "200-0339", "200-0341", "200-0872", "200-0040", "200-0061", "200-0062", "200-0105", "200-0611", "200-0610", "200-0615", "200-0871")

for (i in 1:2838){
  print(i)
  cell_name = ifn_score_2@assays$RNA@data@Dimnames[[2]][i]
  for (j in 1:8297){
    if (cell_name == key_with_disease_status_2[[j,1]]){
      sample_id = key_with_disease_status_2[[j,3]]
      disease_status = key_with_disease_status_2[[j,5]]
      if (disease_status == "Control"){
        sample_id = "LD"
      }
      if (sample_id %in% exclude){
        sample_id = "LD"
        print("LD")
        print(i)
      }
      cell_clus = key_with_disease_status_2[[j,6]] 
      pt_metadata[i,1] = cell_name
      pt_metadata[i,2] = sample_id
      disease_meta[i,1] = cell_name
      disease_meta[i,2] = disease_status
      cell_type[i,1] = cell_name
      cell_type[i,2] = cell_clus
      
      break
    }
  }
}


rownames(pt_metadata) = pt_metadata$cell_name
pt_metadata$cell_name = pt_metadata$patient
pt_metadata = subset(pt_metadata, select = -patient)

rownames(disease_meta) = disease_meta$cell_name
disease_meta$cell_name = disease_meta$status
disease_meta = subset(disease_meta, select = -status)

rownames(cell_type) = cell_type$cell_name
cell_type$cell_name = cell_type$cluster
cell_type = subset(cell_type, select = -cluster)


ifn_score_2 = AddMetaData(object = ifn_score_2, metadata = pt_metadata, col.name = 'Patient')
ifn_score_2 = AddMetaData(ifn_score_2, disease_meta, col.name="Disease")
ifn_score_2 = AddMetaData(ifn_score_2, cell_type, col.name="Cluster")


# IFN score goes in the first gene slot because it has to go somewhere. Lose A1BG data, but not using it. 
VlnPlot(ifn_score_2, features = c("A1BG"), slot="scale.data", group.by="Patient")
VlnPlot(ifn_score_2, features = c("A1BG"), slot="scale.data", group.by="Disease")
VlnPlot(ifn_score_2, features = c("A1BG"), slot="scale.data", group.by="Cluster")

VlnPlot(ifn_score_2, features = c("XIST"), slot="scale.data", group.by="Patient")
VlnPlot(ifn_score_2, features = c("XIST"), slot="scale.data", group.by="Disease")
VlnPlot(ifn_score_2, features = c("XIST"), slot="scale.data", group.by="Cluster")

# Make an expression matrix including interferon scores
exp_mat = FetchData(ifn_score, slot="scale.data", vars = c("IFN.Module1", all.genes))
exp_mat = t(exp_mat)

## Make sure we correctly ID's patients contributing > 10 cells
key_with_disease_status = read.table(file.choose()) # Choose celseq_meta.tsv file
patient_list = data.frame(id=character(0), num_cells=numeric(0))
patient_list[1,1] = "200-0061"
patient_list[1,2] = 0

for (i in 1:8297){
  if (key_with_disease_status[i,3] %in% patient_list$id){
    for (j in 1:nrow(patient_list)){
      if (key_with_disease_status[i,3] == patient_list[j,1]){
        patient_list[j,2] = patient_list[j,2] + 1
      }
    }
  } else {
    patient_list[nrow(patient_list)+1,1] = key_with_disease_status[i,3]
    patient_list[nrow(patient_list),2] = 1
  }
}

# Now find out how many cells are in each patient in the actual published data
count_061=0
count_0608=0
count_0609=0
count_109=0
count_871=0
count_613=0
count_611=0
count_114=0
count_614=0
count_872=0
count_616=0
count_617=0
count_362=0
count_961=0
count_612=0
count_610=0
count_107=0
count_391=0
# count_105=0 # Only 5 possible cells
count_341=0
count_040=0
count_062=0
count_339=0
count_841=0
count_874=0
count_108=0
count_124=0
count_1216=0
count_1217=0
count_873=0
count_875=0
count_127=0
# count_104=0 # Only 3 possible cells
count_095=0
# Now let's find out how many cells from each subject are actually in the dataset
for (i in 1:2839){
  print(i)
  cell_code = all_amp_data[1,i]
  for (j in 1:8297){
    if (cell_code == key_with_disease_status[j,1]){
      if (key_with_disease_status[j,3] == "200-0061"){
        count_061=count_061 + 1
      }
      if (key_with_disease_status[j,3] == "200-0608"){
        count_0608=count_0608 + 1
      }
      if (key_with_disease_status[j,3] == "200-0609"){
        count_0609=count_0609 + 1
      }
      if (key_with_disease_status[j,3] == "200-0109"){
        count_109=count_109 + 1
      }
      if (key_with_disease_status[j,3] == "200-0871"){
        count_871=count_871 + 1
      }
      if (key_with_disease_status[j,3] == "200-0613"){
        count_613=count_613 + 1
      }
      if (key_with_disease_status[j,3] == "200-0611"){
        count_611=count_611 + 1
      }
      if (key_with_disease_status[j,3] == "200-0114"){
        count_114=count_114 + 1
      }
      if (key_with_disease_status[j,3] == "200-0614"){
        count_614=count_614 + 1
      }
      if (key_with_disease_status[j,3] == "200-0872"){
        count_872=count_872 + 1
      }
      if (key_with_disease_status[j,3] == "200-0616"){
        count_616=count_616 + 1
      }
      if (key_with_disease_status[j,3] == "200-0617"){
        count_617=count_617 + 1
      }
      if (key_with_disease_status[j,3] == "200-0362"){
        count_362 = count_362+1
      }
      if (key_with_disease_status[j,3] == "200-0961"){
        count_961=count_961+1
      }
      if (key_with_disease_status[j,3] == "200-0612"){
        count_612=count_612+1
      }
      if (key_with_disease_status[j,3] == "200-0610"){
        count_610=count_610+1
      }
      if (key_with_disease_status[j,3] == "200-0107"){
        count_107=count_107+1
      }
      if (key_with_disease_status[j,3] == "200-0391"){
        count_391=count_391+1
      }
      if (key_with_disease_status[j,3] == "200-0341"){
        count_341=count_341+1
      }
      if (key_with_disease_status[j,3] == "200-0040"){
        count_040=count_040+1
      }
      if (key_with_disease_status[j,3] == "200-0062"){
        count_062=count_062+1
      }
      if (key_with_disease_status[j,3] == "200-0339"){
        count_339=count_339+1
      }
      if (key_with_disease_status[j,3] == "200-0841"){
        count_841=count_841+1
      }
      if (key_with_disease_status[j,3] == "200-0874"){
        count_874=count_874+1
      }
      if (key_with_disease_status[j,3] == "200-0108"){
        count_108=count_108+1
      }
      if (key_with_disease_status[j,3] == "200-0124"){
        count_124=count_124+1
      }
      if (key_with_disease_status[j,3] == "200-1216"){
        count_1216=count_1216+1
      }
      if (key_with_disease_status[j,3] == "200-1217"){
        count_1217=count_1217+1
      }
      if (key_with_disease_status[j,3] == "200-0873"){
        count_873=count_873+1
      }
      if (key_with_disease_status[j,3] == "200-0875"){
        count_875=count_875+1
      }
      if (key_with_disease_status[j,3] == "200-0127"){
        count_127=count_127+1
      }
    }
  }
  
}

# At least 10 cells: 609, 061, 108, 109, 114, 1216, 1217, 124, 339, 362, 391, 612, 613, 614, 616, 617, 841, 872, 873, 875, 961
# Remove not female or not SLE: 609, 339, 613, 614, 616, 617, 872
# Patients: 609, 061, 108, 109, 114, 1216, 1217, 124, 362, 391, 612, 841, 873, 875, 961

### Make matrices for individual patients

sle_1217 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_1216 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_961 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_875 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_874 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_873 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_841 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_612 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_362 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_361 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_114 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_109 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_061 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_108 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sle_124 = data.frame(cell_name=character(0), stringsAsFactors = FALSE)

cell_codes = all_amp_data[1,]
for (i in 1:8297){
  if (key_with_disease_status[i,4] != "Epithelial"){
    if (key_with_disease_status[i,3] == "200-0061"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_061 = rbind(sle_061, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0108"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_108 = rbind(sle_108, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0109"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_109 = rbind(sle_109, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0114"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_114 = rbind(sle_114, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0124"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_124 = rbind(sle_124, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0361"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_361 = rbind(sle_361, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0362"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_362 = rbind(sle_362, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0612"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_612 = rbind(sle_612, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0841"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_841 = rbind(sle_841, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0873"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_873 = rbind(sle_873, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0874"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_874 = rbind(sle_874, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0875"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_875 = rbind(sle_875, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-0961"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_961 = rbind(sle_961, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-1216"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_1216 = rbind(sle_1216, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
    if (key_with_disease_status[i,3] == "200-1217"){
      if (key_with_disease_status[i,1] %in% cell_codes){
        sle_1217 = rbind(sle_1217, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
      }
    }
  }
}
sle_061_expression = matrix(0,22448,nrow(sle_061))
sle_109_expression = matrix(0,22448,nrow(sle_109))
sle_114_expression = matrix(0,22448,nrow(sle_114))
sle_361_expression = matrix(0,22448,nrow(sle_361))
sle_362_expression = matrix(0,22448,nrow(sle_362))
sle_612_expression = matrix(0,22448,nrow(sle_612))
sle_841_expression = matrix(0,22448,nrow(sle_841))
sle_873_expression = matrix(0,22448,nrow(sle_873))
sle_874_expression = matrix(0,22448,nrow(sle_874))
sle_875_expression = matrix(0,22448,nrow(sle_875))
sle_961_expression = matrix(0,22448,nrow(sle_961))
sle_1216_expression = matrix(0,22448,nrow(sle_1216))
sle_1217_expression = matrix(0,22448,nrow(sle_1217))
sle_108_expression = matrix(0,22448,nrow(sle_108))
sle_124_expression = matrix(0,22448,nrow(sle_124))

rownames(sle_061_expression) = rownames(exp_mat)
rownames(sle_109_expression) = rownames(exp_mat)
rownames(sle_114_expression) = rownames(exp_mat)
rownames(sle_361_expression) = rownames(exp_mat)
rownames(sle_362_expression) = rownames(exp_mat)
rownames(sle_612_expression) = rownames(exp_mat)
rownames(sle_841_expression) = rownames(exp_mat)
rownames(sle_873_expression) = rownames(exp_mat)
rownames(sle_874_expression) = rownames(exp_mat)
rownames(sle_875_expression) = rownames(exp_mat)
rownames(sle_961_expression) = rownames(exp_mat)
rownames(sle_1216_expression) = rownames(exp_mat)
rownames(sle_1217_expression) = rownames(exp_mat)
rownames(sle_108_expression) = rownames(exp_mat)
rownames(sle_124_expression) = rownames(exp_mat)

sle_061_counter=0
sle_108_counter=0
sle_109_counter=0
sle_114_counter=0
sle_361_counter=0
sle_362_counter=0
sle_612_counter=0
sle_841_counter=0
sle_873_counter=0
sle_874_counter=0
sle_875_counter=0
sle_961_counter=0
sle_1216_counter=0
sle_1217_counter=0
sle_108_counter=0
sle_124_counter=0

for (i in 1:2838){
  print(i)
  if (grepl(colnames(exp_mat)[i], sle_061)){
    sle_061_counter = sle_061_counter + 1
    sle_061_expression[,sle_061_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_108)){
    sle_108_counter = sle_108_counter + 1
    sle_108_expression[,sle_108_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_109)){
    sle_109_counter = sle_109_counter + 1
    sle_109_expression[,sle_109_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_114)){
    sle_114_counter = sle_114_counter + 1
    sle_114_expression[,sle_114_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_361)){
    sle_361_counter = sle_361_counter + 1
    sle_361_expression[,sle_361_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_362)){
    sle_362_counter = sle_362_counter + 1
    sle_362_expression[,sle_362_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_612)){
    sle_612_counter = sle_612_counter + 1
    sle_612_expression[,sle_612_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_841)){
    sle_841_counter = sle_841_counter + 1
    sle_841_expression[,sle_841_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_873)){
    sle_873_counter = sle_873_counter + 1
    sle_873_expression[,sle_873_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_874)){
    sle_874_counter = sle_874_counter + 1
    sle_874_expression[,sle_874_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_875)){
    sle_875_counter = sle_875_counter + 1
    sle_875_expression[,sle_875_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_961)){
    sle_961_counter = sle_961_counter + 1
    sle_961_expression[,sle_961_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_1216)){
    sle_1216_counter = sle_1216_counter + 1
    sle_1216_expression[,sle_1216_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_1217)){
    sle_1217_counter = sle_1217_counter + 1
    sle_1217_expression[,sle_1217_counter] = exp_mat[,i]
  }
  if (grepl(colnames(exp_mat)[i], sle_124)){
    sle_124_counter = sle_124_counter + 1
    sle_124_expression[,sle_124_counter] = exp_mat[,i]
    
  }
}

sle_108_expression = as.data.frame(sle_108_expression)
sle_124_expression = as.data.frame(sle_124_expression)
sle_1217_expression = as.data.frame(sle_1217_expression)
sle_1216_expression = as.data.frame(sle_1216_expression)
sle_961_expression = as.data.frame(sle_961_expression)
sle_875_expression = as.data.frame(sle_875_expression)
sle_874_expression = as.data.frame(sle_874_expression)
sle_873_expression = as.data.frame(sle_873_expression)
sle_841_expression = as.data.frame(sle_841_expression)
sle_612_expression = as.data.frame(sle_612_expression)
sle_362_expression = as.data.frame(sle_362_expression)
sle_361_expression = as.data.frame(sle_361_expression)
sle_114_expression = as.data.frame(sle_114_expression)
sle_109_expression = as.data.frame(sle_109_expression)
sle_061_expression = as.data.frame(sle_061_expression)


sle_108_expression$average=rowMeans(sle_108_expression)
sle_124_expression$average=rowMeans(sle_124_expression)
sle_1217_expression$average=rowMeans(sle_1217_expression)
sle_1216_expression$average=rowMeans(sle_1216_expression)
sle_961_expression$average=rowMeans(sle_961_expression)
sle_875_expression$average=rowMeans(sle_875_expression)
sle_874_expression$average=rowMeans(sle_874_expression)
sle_873_expression$average=rowMeans(sle_873_expression)
sle_841_expression$average=rowMeans(sle_841_expression)
sle_612_expression$average=rowMeans(sle_612_expression)
sle_362_expression$average=rowMeans(sle_362_expression)
sle_361_expression$average=rowMeans(sle_361_expression)
sle_114_expression$average=rowMeans(sle_114_expression)
sle_109_expression$average=rowMeans(sle_109_expression)
sle_061_expression$average=rowMeans(sle_061_expression)

for (i in 1:22448){
  if (rownames(sle_114_expression)[i] == "XIST"){
    print(i)
  }
}

for (i in 1:22448){
  if (rownames(sle_114_expression)[i] == "IL12B"){
    print(i)
  }
}


ifn_scores = c(sle_108_expression$average[1], sle_124_expression$average[1], sle_1217_expression$average[1], sle_1216_expression$average[1], sle_961_expression$average[1], sle_875_expression$average[1], sle_874_expression$average[1], sle_873_expression$average[1], sle_841_expression$average[1], sle_612_expression$average[1], sle_362_expression$average[1], sle_361_expression$average[1], sle_114_expression$average[1], sle_109_expression$average[1], sle_061_expression$average[1])
xist_expression = c(sle_108_expression$average[21551], sle_124_expression$average[21551], sle_1217_expression$average[21551], sle_1216_expression$average[21551], sle_961_expression$average[21551], sle_875_expression$average[21551], sle_874_expression$average[21551], sle_873_expression$average[21551], sle_841_expression$average[21551], sle_612_expression$average[21551], sle_362_expression$average[21551], sle_361_expression$average[21551], sle_114_expression$average[21551], sle_109_expression$average[21551], sle_061_expression$average[21551])







# Generate Fig 5A
plot(xist_expression, ifn_scores)
cor.test(xist_expression, ifn_scores)

### Are other genes correlated with the IFN signature?

pos_counter=0
neg_counter=0
pos_cor_genes = c()
neg_cor_genes = c()
for (i in 1:22448){ 
  print(i)
  expression = c(sle_108_expression$average[i], sle_124_expression$average[i], sle_1217_expression$average[i], sle_1216_expression$average[i], sle_961_expression$average[i], sle_875_expression$average[i], sle_874_expression$average[i], sle_873_expression$average[i], sle_841_expression$average[i], sle_612_expression$average[i], sle_362_expression$average[i], sle_361_expression$average[i], sle_114_expression$average[i], sle_109_expression$average[i], sle_061_expression$average[i])
  result = cor.test(expression, ifn_scores)
  if (! is.na(result$p.value)){
    if (result$p.value < 0.05){
      if (result$estimate < 0){
        neg_counter=neg_counter+1
        neg_cor_genes = c(neg_cor_genes, rownames(sle_108_expression)[i])
      } else if (result$estimate > 0){
        pos_counter=pos_counter+1
        pos_cor_genes = c(pos_cor_genes, rownames(sle_108_expression)[i])
      }
    }
  }
}


plot(exp_mat[21551,], exp_mat[1,])
new = cor.test(exp_mat[21551,], exp_mat[1,])

rownames(exp_mat)[21551]
rownames(exp_mat)[1]


pos_sc_counter=0
neg_sc_counter=0
for (i in 1:22448){
  print(i)
  new = cor.test(exp_mat[i,], exp_mat[1,])
  if (new$p.value < 0.05){
    if (new$estimate > 0){
      pos_sc_counter = pos_sc_counter+1
    }
    if (new$estimate < 0){
      neg_sc_counter = neg_sc_counter+1
    }
  }
}

##############################################################################
### Now Supplementary Figure 2
##############################################################################

library(readxl)
xist_kd_all_data <- read_excel(file.choose()) # XIST KD RNA Sequencing Correct

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

p2 <- ggplot(data=small_df, aes(x=log2_fold_change, y=-log10(adj_p), col=diffexpressed, label=label )) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()
mycolors <- c("red", "blue", "black")
names(mycolors) <- c("UP", "DOWN", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3
p4 <- p3 + geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p4



