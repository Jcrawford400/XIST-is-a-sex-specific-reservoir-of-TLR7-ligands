load("C:/Users/13174/Downloads/2023_05_15_edgeRUU.rda")

myvars <- c("gene", "padj", "log2FCtpm", "female_tpm", "male_tpm", "UU.Count", "UU.Richness", "has_GTCCTTCAA")
################################################################################
###
################################################################################

blood <- edgeR_UU$Blood[myvars]
blood$sbrank = 0
blood$uucountrank = 0
blood$uurichrank = 0
blood$tpmrank = 0
blood$ranksum = 0
blood_copy = blood
blood = blood[order(blood$log2FCtpm),]
for (i in 1:nrow(blood)){
  blood$sbrank[i] = i
}
blood = blood[order(blood$UU.Count),]
for (i in 1:nrow(blood)){
  blood$uucountrank[i] = i
}
blood = blood[order(blood$UU.Richness),]
for (i in 1:nrow(blood)){
  blood$uurichrank[i] = i
}
blood = blood[order(blood$female_tpm),]
for (i in 1:nrow(blood)){
  blood$tpmrank[i] = i
}
blood$ranksum = blood$sbrank + blood$uucountrank + blood$uurichrank + blood$tpmrank + blood$ranksum


blood_sex_bias = subset(blood, blood$padj < 0.05)

library("writexl")
path <- "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data"
filename = paste0(path, '/blood_sex_bias.xlsx')
write_xlsx(as.data.frame(blood_sex_bias), path = filename)

filename = paste0(path, '/blood.xlsx')
write_xlsx(as.data.frame(blood), path = filename)



###################################################################################################
### Let's try adipose tissue
###################################################################################################

adipose <- edgeR_UU$`Adipose Tissue`[myvars]
adipose$sbrank = 0
adipose$uucountrank = 0
adipose$uurichrank = 0
adipose$tpmrank = 0
adipose$ranksum = 0
adipose_copy = adipose
adipose = adipose[order(adipose$log2FCtpm),]
for (i in 1:nrow(adipose)){
  adipose$sbrank[i] = i
}
adipose = adipose[order(adipose$UU.Count),]
for (i in 1:nrow(adipose)){
  adipose$uucountrank[i] = i
}
adipose = adipose[order(adipose$UU.Richness),]
for (i in 1:nrow(adipose)){
  adipose$uurichrank[i] = i
}
adipose = adipose[order(adipose$female_tpm),]
for (i in 1:nrow(adipose)){
  adipose$tpmrank[i] = i
}
adipose$ranksum = adipose$sbrank + adipose$uucountrank + adipose$uurichrank + adipose$tpmrank + adipose$ranksum

###################################################################################################
### Let's try Spleen
###################################################################################################

spleen <- edgeR_UU$Spleen[myvars]
spleen$sbrank = 0
spleen$uucountrank = 0
spleen$uurichrank = 0
spleen$tpmrank = 0
spleen$ranksum = 0
spleen_copy = spleen
spleen = spleen[order(spleen$log2FCtpm),]
for (i in 1:nrow(spleen)){
  spleen$sbrank[i] = i
}
spleen = spleen[order(spleen$UU.Count),]
for (i in 1:nrow(spleen)){
  spleen$uucountrank[i] = i
}
spleen = spleen[order(spleen$UU.Richness),]
for (i in 1:nrow(spleen)){
  spleen$uurichrank[i] = i
}
spleen = spleen[order(spleen$female_tpm),]
for (i in 1:nrow(spleen)){
  spleen$tpmrank[i] = i
}
spleen$ranksum = spleen$sbrank + spleen$uucountrank + spleen$uurichrank + spleen$tpmrank + spleen$ranksum


# myvars <- c("gene", "ranksum")
# bloodranksums <- blood[myvars]
# filename = paste0(path, '/bloodranksums.xlsx')
# write_xlsx(as.data.frame(bloodranksums), path = filename)
# 
# bloodranksums$adipose = 0
# for (i in 1:nrow(bloodranksums)){
#   print(i)
#   gene = bloodranksums$gene[i]
#   for (j in 1:nrow(adipose)){
#     if (adipose$gene[j] == gene){
#       bloodranksums$adipose[i] = adipose$ranksum[j]
#     }
#   }
# }


filename = paste0(path, '/spleen.xlsx')
write_xlsx(as.data.frame(spleen), path = filename)
###################################################################################################
### Kidney
###################################################################################################
kidney <- edgeR_UU$Kidney[myvars]
kidney$sbrank = 0
kidney$uucountrank = 0
kidney$uurichrank = 0
kidney$tpmrank = 0
kidney$ranksum = 0
kidney_copy = kidney
kidney = kidney[order(kidney$log2FCtpm),]
for (i in 1:nrow(kidney)){
  kidney$sbrank[i] = i
}
kidney = kidney[order(kidney$UU.Count),]
for (i in 1:nrow(kidney)){
  kidney$uucountrank[i] = i
}
kidney = kidney[order(kidney$UU.Richness),]
for (i in 1:nrow(kidney)){
  kidney$uurichrank[i] = i
}
kidney = kidney[order(kidney$female_tpm),]
for (i in 1:nrow(kidney)){
  kidney$tpmrank[i] = i
}
kidney$ranksum = kidney$sbrank + kidney$uucountrank + kidney$uurichrank + kidney$tpmrank + kidney$ranksum

filename = paste0(path, '/kidney.xlsx')
write_xlsx(as.data.frame(kidney), path = filename)
################################################################################################
### Liver
################################################################################################

liver <- edgeR_UU$Liver[myvars]
liver$sbrank = 0
liver$uucountrank = 0
liver$uurichrank = 0
liver$tpmrank = 0
liver$ranksum = 0
liver_copy = liver
liver = liver[order(liver$log2FCtpm),]
for (i in 1:nrow(liver)){
  liver$sbrank[i] = i
}
liver = liver[order(liver$UU.Count),]
for (i in 1:nrow(liver)){
  liver$uucountrank[i] = i
}
liver = liver[order(liver$UU.Richness),]
for (i in 1:nrow(liver)){
  liver$uurichrank[i] = i
}
liver = liver[order(liver$female_tpm),]
for (i in 1:nrow(liver)){
  liver$tpmrank[i] = i
}
liver$ranksum = liver$sbrank + liver$uucountrank + liver$uurichrank + liver$tpmrank + liver$ranksum

################################################################################################
### Muscle
################################################################################################

muscle <- edgeR_UU$Muscle[myvars]
muscle$sbrank = 0
muscle$uucountrank = 0
muscle$uurichrank = 0
muscle$tpmrank = 0
muscle$ranksum = 0
muscle_copy = muscle
muscle = muscle[order(muscle$log2FCtpm),]
for (i in 1:nrow(muscle)){
  muscle$sbrank[i] = i
}
muscle = muscle[order(muscle$UU.Count),]
for (i in 1:nrow(muscle)){
  muscle$uucountrank[i] = i
}
muscle = muscle[order(muscle$UU.Richness),]
for (i in 1:nrow(muscle)){
  muscle$uurichrank[i] = i
}
muscle = muscle[order(muscle$female_tpm),]
for (i in 1:nrow(muscle)){
  muscle$tpmrank[i] = i
}
muscle$ranksum = muscle$sbrank + muscle$uucountrank + muscle$uurichrank + muscle$tpmrank + muscle$ranksum

################################################################################################
### Blood Vessel
################################################################################################

blood_vessel <- edgeR_UU$`Blood Vessel`[myvars]
blood_vessel$sbrank = 0
blood_vessel$uucountrank = 0
blood_vessel$uurichrank = 0
blood_vessel$tpmrank = 0
blood_vessel$ranksum = 0
blood_vessel_copy = blood_vessel
blood_vessel = blood_vessel[order(blood_vessel$log2FCtpm),]
for (i in 1:nrow(blood_vessel)){
  blood_vessel$sbrank[i] = i
}
blood_vessel = blood_vessel[order(blood_vessel$UU.Count),]
for (i in 1:nrow(blood_vessel)){
  blood_vessel$uucountrank[i] = i
}
blood_vessel = blood_vessel[order(blood_vessel$UU.Richness),]
for (i in 1:nrow(blood_vessel)){
  blood_vessel$uurichrank[i] = i
}
blood_vessel = blood_vessel[order(blood_vessel$female_tpm),]
for (i in 1:nrow(blood_vessel)){
  blood_vessel$tpmrank[i] = i
}
blood_vessel$ranksum = blood_vessel$sbrank + blood_vessel$uucountrank + blood_vessel$uurichrank + blood_vessel$tpmrank + blood_vessel$ranksum

################################################################################################
### Heart
################################################################################################

heart <- edgeR_UU$Heart[myvars]
heart$sbrank = 0
heart$uucountrank = 0
heart$uurichrank = 0
heart$tpmrank = 0
heart$ranksum = 0
heart_copy = heart
heart = heart[order(heart$log2FCtpm),]
for (i in 1:nrow(heart)){
  heart$sbrank[i] = i
}
heart = heart[order(heart$UU.Count),]
for (i in 1:nrow(heart)){
  heart$uucountrank[i] = i
}
heart = heart[order(heart$UU.Richness),]
for (i in 1:nrow(heart)){
  heart$uurichrank[i] = i
}
heart = heart[order(heart$female_tpm),]
for (i in 1:nrow(heart)){
  heart$tpmrank[i] = i
}
heart$ranksum = heart$sbrank + heart$uucountrank + heart$uurichrank + heart$tpmrank + heart$ranksum

################################################################################################
### Skin
################################################################################################

skin <- edgeR_UU$Skin[myvars]
skin$sbrank = 0
skin$uucountrank = 0
skin$uurichrank = 0
skin$tpmrank = 0
skin$ranksum = 0
skin_copy = skin
skin = skin[order(skin$log2FCtpm),]
for (i in 1:nrow(skin)){
  skin$sbrank[i] = i
}
skin = skin[order(skin$UU.Count),]
for (i in 1:nrow(skin)){
  skin$uucountrank[i] = i
}
skin = skin[order(skin$UU.Richness),]
for (i in 1:nrow(skin)){
  skin$uurichrank[i] = i
}
skin = skin[order(skin$female_tpm),]
for (i in 1:nrow(skin)){
  skin$tpmrank[i] = i
}
skin$ranksum = skin$sbrank + skin$uucountrank + skin$uurichrank + skin$tpmrank + skin$ranksum

################################################################################################
### Salivary Gland
################################################################################################

salgland <- edgeR_UU$`Salivary Gland`[myvars]
salgland$sbrank = 0
salgland$uucountrank = 0
salgland$uurichrank = 0
salgland$tpmrank = 0
salgland$ranksum = 0
salgland_copy = salgland
salgland = salgland[order(salgland$log2FCtpm),]
for (i in 1:nrow(salgland)){
  salgland$sbrank[i] = i
}
salgland = salgland[order(salgland$UU.Count),]
for (i in 1:nrow(salgland)){
  salgland$uucountrank[i] = i
}
salgland = salgland[order(salgland$UU.Richness),]
for (i in 1:nrow(salgland)){
  salgland$uurichrank[i] = i
}
salgland = salgland[order(salgland$female_tpm),]
for (i in 1:nrow(salgland)){
  salgland$tpmrank[i] = i
}
salgland$ranksum = salgland$sbrank + salgland$uucountrank + salgland$uurichrank + salgland$tpmrank + salgland$ranksum

################################################################################################
### Liver
################################################################################################

brain <- edgeR_UU$Brain[myvars]
brain$sbrank = 0
brain$uucountrank = 0
brain$uurichrank = 0
brain$tpmrank = 0
brain$ranksum = 0
brain_copy = brain
brain = brain[order(brain$log2FCtpm),]
for (i in 1:nrow(brain)){
  brain$sbrank[i] = i
}
brain = brain[order(brain$UU.Count),]
for (i in 1:nrow(brain)){
  brain$uucountrank[i] = i
}
brain = brain[order(brain$UU.Richness),]
for (i in 1:nrow(brain)){
  brain$uurichrank[i] = i
}
brain = brain[order(brain$female_tpm),]
for (i in 1:nrow(brain)){
  brain$tpmrank[i] = i
}
brain$ranksum = brain$sbrank + brain$uucountrank + brain$uurichrank + brain$tpmrank + brain$ranksum

################################################################################################
### Adrenal Gland
################################################################################################

adrgland <- edgeR_UU$`Adrenal Gland`[myvars]
adrgland$sbrank = 0
adrgland$uucountrank = 0
adrgland$uurichrank = 0
adrgland$tpmrank = 0
adrgland$ranksum = 0
adrgland_copy = adrgland
adrgland = adrgland[order(adrgland$log2FCtpm),]
for (i in 1:nrow(adrgland)){
  adrgland$sbrank[i] = i
}
adrgland = adrgland[order(adrgland$UU.Count),]
for (i in 1:nrow(adrgland)){
  adrgland$uucountrank[i] = i
}
adrgland = adrgland[order(adrgland$UU.Richness),]
for (i in 1:nrow(adrgland)){
  adrgland$uurichrank[i] = i
}
adrgland = adrgland[order(adrgland$female_tpm),]
for (i in 1:nrow(adrgland)){
  adrgland$tpmrank[i] = i
}
adrgland$ranksum = adrgland$sbrank + adrgland$uucountrank + adrgland$uurichrank + adrgland$tpmrank + adrgland$ranksum

################################################################################################
### Thyroid
################################################################################################

thyroid <- edgeR_UU$Thyroid[myvars]
thyroid$sbrank = 0
thyroid$uucountrank = 0
thyroid$uurichrank = 0
thyroid$tpmrank = 0
thyroid$ranksum = 0
thyroid_copy = thyroid
thyroid = thyroid[order(thyroid$log2FCtpm),]
for (i in 1:nrow(thyroid)){
  thyroid$sbrank[i] = i
}
thyroid = thyroid[order(thyroid$UU.Count),]
for (i in 1:nrow(thyroid)){
  thyroid$uucountrank[i] = i
}
thyroid = thyroid[order(thyroid$UU.Richness),]
for (i in 1:nrow(thyroid)){
  thyroid$uurichrank[i] = i
}
thyroid = thyroid[order(thyroid$female_tpm),]
for (i in 1:nrow(thyroid)){
  thyroid$tpmrank[i] = i
}
thyroid$ranksum = thyroid$sbrank + thyroid$uucountrank + thyroid$uurichrank + thyroid$tpmrank + thyroid$ranksum

################################################################################################
### Lung
################################################################################################

lung <- edgeR_UU$Lung[myvars]
lung$sbrank = 0
lung$uucountrank = 0
lung$uurichrank = 0
lung$tpmrank = 0
lung$ranksum = 0
lung_copy = lung
lung = lung[order(lung$log2FCtpm),]
for (i in 1:nrow(lung)){
  lung$sbrank[i] = i
}
lung = lung[order(lung$UU.Count),]
for (i in 1:nrow(lung)){
  lung$uucountrank[i] = i
}
lung = lung[order(lung$UU.Richness),]
for (i in 1:nrow(lung)){
  lung$uurichrank[i] = i
}
lung = lung[order(lung$female_tpm),]
for (i in 1:nrow(lung)){
  lung$tpmrank[i] = i
}
lung$ranksum = lung$sbrank + lung$uucountrank + lung$uurichrank + lung$tpmrank + lung$ranksum

################################################################################################
### Pancreas
################################################################################################

pancreas <- edgeR_UU$Pancreas[myvars]
pancreas$sbrank = 0
pancreas$uucountrank = 0
pancreas$uurichrank = 0
pancreas$tpmrank = 0
pancreas$ranksum = 0
pancreas_copy = pancreas
pancreas = pancreas[order(pancreas$log2FCtpm),]
for (i in 1:nrow(pancreas)){
  pancreas$sbrank[i] = i
}
pancreas = pancreas[order(pancreas$UU.Count),]
for (i in 1:nrow(pancreas)){
  pancreas$uucountrank[i] = i
}
pancreas = pancreas[order(pancreas$UU.Richness),]
for (i in 1:nrow(pancreas)){
  pancreas$uurichrank[i] = i
}
pancreas = pancreas[order(pancreas$female_tpm),]
for (i in 1:nrow(pancreas)){
  pancreas$tpmrank[i] = i
}
pancreas$ranksum = pancreas$sbrank + pancreas$uucountrank + pancreas$uurichrank + pancreas$tpmrank + pancreas$ranksum

################################################################################################
### Liver
################################################################################################

esophagus <- edgeR_UU$Esophagus[myvars]
esophagus$sbrank = 0
esophagus$uucountrank = 0
esophagus$uurichrank = 0
esophagus$tpmrank = 0
esophagus$ranksum = 0
esophagus_copy = esophagus
esophagus = esophagus[order(esophagus$log2FCtpm),]
for (i in 1:nrow(esophagus)){
  esophagus$sbrank[i] = i
}
esophagus = esophagus[order(esophagus$UU.Count),]
for (i in 1:nrow(esophagus)){
  esophagus$uucountrank[i] = i
}
esophagus = esophagus[order(esophagus$UU.Richness),]
for (i in 1:nrow(esophagus)){
  esophagus$uurichrank[i] = i
}
esophagus = esophagus[order(esophagus$female_tpm),]
for (i in 1:nrow(esophagus)){
  esophagus$tpmrank[i] = i
}
esophagus$ranksum = esophagus$sbrank + esophagus$uucountrank + esophagus$uurichrank + esophagus$tpmrank + esophagus$ranksum

################################################################################################
### Stomach
################################################################################################

stomach <- edgeR_UU$Stomach[myvars]
stomach$sbrank = 0
stomach$uucountrank = 0
stomach$uurichrank = 0
stomach$tpmrank = 0
stomach$ranksum = 0
stomach_copy = stomach
stomach = stomach[order(stomach$log2FCtpm),]
for (i in 1:nrow(stomach)){
  stomach$sbrank[i] = i
}
stomach = stomach[order(stomach$UU.Count),]
for (i in 1:nrow(stomach)){
  stomach$uucountrank[i] = i
}
stomach = stomach[order(stomach$UU.Richness),]
for (i in 1:nrow(stomach)){
  stomach$uurichrank[i] = i
}
stomach = stomach[order(stomach$female_tpm),]
for (i in 1:nrow(stomach)){
  stomach$tpmrank[i] = i
}
stomach$ranksum = stomach$sbrank + stomach$uucountrank + stomach$uurichrank + stomach$tpmrank + stomach$ranksum

################################################################################################
### Colon
################################################################################################

colon <- edgeR_UU$Colon[myvars]
colon$sbrank = 0
colon$uucountrank = 0
colon$uurichrank = 0
colon$tpmrank = 0
colon$ranksum = 0
colon_copy = colon
colon = colon[order(colon$log2FCtpm),]
for (i in 1:nrow(colon)){
  colon$sbrank[i] = i
}
colon = colon[order(colon$UU.Count),]
for (i in 1:nrow(colon)){
  colon$uucountrank[i] = i
}
colon = colon[order(colon$UU.Richness),]
for (i in 1:nrow(colon)){
  colon$uurichrank[i] = i
}
colon = colon[order(colon$female_tpm),]
for (i in 1:nrow(colon)){
  colon$tpmrank[i] = i
}
colon$ranksum = colon$sbrank + colon$uucountrank + colon$uurichrank + colon$tpmrank + colon$ranksum

################################################################################################
### Small Intestine
################################################################################################

smallint <- edgeR_UU$`Small Intestine`[myvars]
smallint$sbrank = 0
smallint$uucountrank = 0
smallint$uurichrank = 0
smallint$tpmrank = 0
smallint$ranksum = 0
smallint_copy = smallint
smallint = smallint[order(smallint$log2FCtpm),]
for (i in 1:nrow(smallint)){
  smallint$sbrank[i] = i
}
smallint = smallint[order(smallint$UU.Count),]
for (i in 1:nrow(smallint)){
  smallint$uucountrank[i] = i
}
smallint = smallint[order(smallint$UU.Richness),]
for (i in 1:nrow(smallint)){
  smallint$uurichrank[i] = i
}
smallint = smallint[order(smallint$female_tpm),]
for (i in 1:nrow(smallint)){
  smallint$tpmrank[i] = i
}
smallint$ranksum = smallint$sbrank + smallint$uucountrank + smallint$uurichrank + smallint$tpmrank + smallint$ranksum

################################################################################################
### Liver
################################################################################################

nerve <- edgeR_UU$Nerve[myvars]
nerve$sbrank = 0
nerve$uucountrank = 0
nerve$uurichrank = 0
nerve$tpmrank = 0
nerve$ranksum = 0
nerve_copy = nerve
nerve = nerve[order(nerve$log2FCtpm),]
for (i in 1:nrow(nerve)){
  nerve$sbrank[i] = i
}
nerve = nerve[order(nerve$UU.Count),]
for (i in 1:nrow(nerve)){
  nerve$uucountrank[i] = i
}
nerve = nerve[order(nerve$UU.Richness),]
for (i in 1:nrow(nerve)){
  nerve$uurichrank[i] = i
}
nerve = nerve[order(nerve$female_tpm),]
for (i in 1:nrow(nerve)){
  nerve$tpmrank[i] = i
}
nerve$ranksum = nerve$sbrank + nerve$uucountrank + nerve$uurichrank + nerve$tpmrank + nerve$ranksum

################################################################################################
### Pituitary
################################################################################################

pituitary <- edgeR_UU$Pituitary[myvars]
pituitary$sbrank = 0
pituitary$uucountrank = 0
pituitary$uurichrank = 0
pituitary$tpmrank = 0
pituitary$ranksum = 0
pituitary_copy = pituitary
pituitary = pituitary[order(pituitary$log2FCtpm),]
for (i in 1:nrow(pituitary)){
  pituitary$sbrank[i] = i
}
pituitary = pituitary[order(pituitary$UU.Count),]
for (i in 1:nrow(pituitary)){
  pituitary$uucountrank[i] = i
}
pituitary = pituitary[order(pituitary$UU.Richness),]
for (i in 1:nrow(pituitary)){
  pituitary$uurichrank[i] = i
}
pituitary = pituitary[order(pituitary$female_tpm),]
for (i in 1:nrow(pituitary)){
  pituitary$tpmrank[i] = i
}
pituitary$ranksum = pituitary$sbrank + pituitary$uucountrank + pituitary$uurichrank + pituitary$tpmrank + pituitary$ranksum

################################################################################################
### Liver
################################################################################################

bladder <- edgeR_UU$Bladder[myvars]
bladder$sbrank = 0
bladder$uucountrank = 0
bladder$uurichrank = 0
bladder$tpmrank = 0
bladder$ranksum = 0
bladder_copy = bladder
bladder = bladder[order(bladder$log2FCtpm),]
for (i in 1:nrow(bladder)){
  bladder$sbrank[i] = i
}
bladder = bladder[order(bladder$UU.Count),]
for (i in 1:nrow(bladder)){
  bladder$uucountrank[i] = i
}
bladder = bladder[order(bladder$UU.Richness),]
for (i in 1:nrow(bladder)){
  bladder$uurichrank[i] = i
}
bladder = bladder[order(bladder$female_tpm),]
for (i in 1:nrow(bladder)){
  bladder$tpmrank[i] = i
}
bladder$ranksum = bladder$sbrank + bladder$uucountrank + bladder$uurichrank + bladder$tpmrank + bladder$ranksum




all_tissues = data.frame(gene=character(0), blood=numeric(0), spleen=numeric(0), kidney=numeric(0), liver=numeric(0), adipose=numeric(0), adrgland=numeric(0), bladder=numeric(0), brain=numeric(0), colon=numeric(0), heart=numeric(0), kidney=numeric(0), pituitary=numeric(0), salgland=numeric(0), skin=numeric(0), thyroid=numeric(0))
all_tissues[1,1] = 'XIST'
for (i in 1:nrow(blood)){
  gene = blood$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(kidney)){
  gene = kidney$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(liver)){
  gene = liver$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(adipose)){
  gene = adipose$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(adrgland)){
  gene = adrgland$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(bladder)){
  gene = bladder$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(brain)){
  gene = brain$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(colon)){
  gene = colon$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(heart)){
  gene = heart$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(muscle)){
  gene = muscle$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(pituitary)){
  gene = pituitary$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(salgland)){
  gene = salgland$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(skin)){
  gene = skin$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}
for (i in 1:nrow(thyroid)){
  gene = thyroid$gene[i]
  if (! gene %in% all_tissues$gene){
    all_tissues[nrow(all_tissues) + 1, 1] = gene
  }
}

############### Compile the ranks ################33

for (i in 1:nrow(all_tissues)){
  print(i)
  gene = all_tissues$gene[i]
  for (j in 1:nrow(blood)){
    if (gene == blood$gene[j]){
      all_tissues$blood[i] = blood$ranksum[j]
      break
    }
  }
  for (k in 1:nrow(spleen)){
    if (gene == spleen$gene[k]){
      all_tissues$spleen[i] = spleen$ranksum[k]
      break
    }
  }
  for (l in 1:nrow(kidney)){
    if (gene == kidney$gene[l]){
      all_tissues$kidney[i] = kidney$ranksum[l]
      break
    }
  }
  for (m in 1:nrow(liver)){
    if (gene == liver$gene[m]){
      all_tissues$liver[i] = liver$ranksum[m]
      break
    }
  }
  for (n in 1:nrow(adipose)){
    if (gene == adipose$gene[n]){
      all_tissues$adipose[i] = adipose$ranksum[n]
      break
    }
  }
  for (o in 1:nrow(adrgland)){
    if (gene == adrgland$gene[o]){
      all_tissues$adrgland[i] = adrgland$ranksum[o]
      break
    }
  }
  for (p in 1:nrow(bladder)){
    if (gene == bladder$gene[p]){
      all_tissues$bladder[i] = bladder$ranksum[p]
      break
    }
  }
  for (q in 1:nrow(brain)){
    if (gene == brain$gene[q]){
      all_tissues$brain[i] = brain$ranksum[q]
      break
    }
  }
  for (r in 1:nrow(colon)){
    if (gene == colon$gene[r]){
      all_tissues$colon[i] = colon$ranksum[r]
      break
    }
  }
  for (s in 1:nrow(heart)){
    if (gene == heart$gene[s]){
      all_tissues$heart[i] = heart$ranksum[s]
      break
    }
  }
  for (y in 1:nrow(muscle)){
    if (gene == heart$gene[y]){
      all_tissues$muscle[i] = muscle$ranksum[y]
      break
    }
  }
  for (t in 1:nrow(pituitary)){
    if (gene == pituitary$gene[t]){
      all_tissues$pituitary[i] = pituitary$ranksum[t]
      break
    }
  }
  for (u in 1:nrow(salgland)){
    if (gene == salgland$gene[u]){
      all_tissues$salgland[i] = salgland$ranksum[u]
      break
    }
  }
  for (v in 1:nrow(skin)){
    if (gene == skin$gene[v]){
      all_tissues$skin[i] = skin$ranksum[v]
      break
    }
  }
  for (w in 1:nrow(spleen)){
    if (gene == spleen$gene[w]){
      all_tissues$spleen[i] = spleen$ranksum[w]
      break
    }
  }
  for (x in 1:nrow(thyroid)){
    if (gene == thyroid$gene[x]){
      all_tissues$thyroid[i] = thyroid$ranksum[x]
      break
    }
  }
}

for (i in 1:19349){
  if (is.na(all_tissues$blood[i])){
    all_tissues$blood[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$adipose[i])){
    all_tissues$adipose[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$adrgland[i])){
    all_tissues$adrgland[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$bladder[i])){
    all_tissues$bladder[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$brain[i])){
    all_tissues$brain[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$colon[i])){
    all_tissues$colon[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$esophagus[i])){
    all_tissues$esophagus[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$heart[i])){
    all_tissues$heart[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$kidney[i])){
    all_tissues$kidney[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$liver[i])){
    all_tissues$liver[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$lung[i])){
    all_tissues$lung[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$muscle[i])){
    all_tissues$muscle[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$nerve[i])){
    all_tissues$nerve[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$pancreas[i])){
    all_tissues$pancreas[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$pituitary[i])){
    all_tissues$pituitary[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$salgland[i])){
    all_tissues$salgland[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$skin[i])){
    all_tissues$skin[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$smallint[i])){
    all_tissues$smallint[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$spleen[i])){
    all_tissues$spleen[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$stomach[i])){
    all_tissues$stomach[i] = 0
  }
}
for (i in 1:19349){
  if (is.na(all_tissues$thyroid[i])){
    all_tissues$thyroid[i] = 0
  }
}

path <- "C:/Users/13174/OneDrive - Johns Hopkins/Documents/Darrah Lab/XIST Paper Primary Data"
filename = paste0(path, '/all_tissues.xlsx')
write_xlsx(as.data.frame(all_tissues), path = filename)














all_tissues_percents = all_tissues
all_tissues_percents[is.na(all_tissues_percents)] = 0

all_tissues_percents$blood = all_tissues_percents$blood/max(all_tissues_percents$blood)
all_tissues_percents$blood = all_tissues_percents$blood/max(all_tissues_percents$blood)



