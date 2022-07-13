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

# all_amp_data=backup
amp_seu = CreateSeuratObject(counts = amp_mat, meta.data = amp_meta)
amp_seu = NormalizeData(amp_seu, normalization.method = "LogNormalize", scale.factor=10000)
amp_seu <- FindVariableFeatures(amp_seu, selection.method = "vst", nfeatures = 22709)
amp_seu <- ScaleData(amp_seu)
all.genes=row.names(amp_seu)

# Use AMP interferon score
ifn = c("ABCA1", "ACTA1", "ACTA2", "ADAR", "AIM2", "APOL6", "ASPRV1", "BATF2", "BST2", "BTN3A1", "CARD16", "CARD17", "CCL8", "CEACAM1", "CHMP5", "CMPK2", "CPT1B", "CXCL10", "DDX58", "DDX60", "DHRS9", "DHX58", "DYNLT1", "EIF2AK2", "FBX06", "GADD45B", "GALM", "GBP1", "GBP2", "GBP3", "GBP4", "GBP6", "HELZ2", "HERC5", "HERC6", "HES4", "HSH2D", "IDO1", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM3", "IRF7", "IRF9", "ISG20", "LAP3", "LGALS9", "LY6E", "MOV10", "MT1A", "MX1", "NBN", "NCOA7", "NMI", "NT5C3A", "OAS1", "OAS2", "OAS3", "OASL", "OTOF", "PARP10", "PARP12", "PARP14", "PARP9", "PHF11", "PSMB9", "RBCK1", "REC8", "RHBDF2", "RNF213", "RSAD2", "SERPING1", "SOCS1", "SP100", "SP110", "SP140", "SPATS2L", "SRBD1", "STAT1", "STAT2", "TAP1", "TAP2", "TCN2", "TIMM10", "TMEM140", "TNFAIP6", "TNFSF10", "TRAFD1", "TRANK1", "TRIM21", "TRIM22", "TRIM38", "TRIM5", "TRIM56", "TRIM6", "TYMP", "UBE2L6", "UNC93B1", "WARS", "XAF1", "ZBP1", "ZC3HAV1", "ZNF684", "ZNFX1")
ifn_score = AddModuleScore(amp_seu, features = list(ifn), name="IFN.Module")

exp_mat = FetchData(ifn_score, slot="scale.data", vars = c("IFN.Module1", all.genes))
exp_mat = t(exp_mat)

## Let's find all patients contributing > 10 cells

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

ifn_scores = c(sle_108_expression$average[1], sle_124_expression$average[1], sle_1217_expression$average[1], sle_1216_expression$average[1], sle_961_expression$average[1], sle_875_expression$average[1], sle_874_expression$average[1], sle_873_expression$average[1], sle_841_expression$average[1], sle_612_expression$average[1], sle_362_expression$average[1], sle_361_expression$average[1], sle_114_expression$average[1], sle_109_expression$average[1], sle_061_expression$average[1])
xist_expression = c(sle_108_expression$average[21551], sle_124_expression$average[21551], sle_1217_expression$average[21551], sle_1216_expression$average[21551], sle_961_expression$average[21551], sle_875_expression$average[21551], sle_874_expression$average[21551], sle_873_expression$average[21551], sle_841_expression$average[21551], sle_612_expression$average[21551], sle_362_expression$average[21551], sle_361_expression$average[21551], sle_114_expression$average[21551], sle_109_expression$average[21551], sle_061_expression$average[21551])

# Generate Fig 4G
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







###############################################################################

