###############################################################################
### Figure 5A
###############################################################################

female_sle <- read.table(file = "clipboard", 
                         sep = "\t", header=FALSE) # Use Female SLE data from original authors 
cell_type_identifiers <- read.table(file = "clipboard", 
                                    sep = "\t", header=TRUE) # Use cluster_per_cell sheet from original authors
epithelial_cell_levels = c()
myeloid_cell_levels = c()
t_cell_levels = c()
b_cell_levels = c()
dividing_cell_levels = c()


for (i in 1:2838){
  cell_name = cell_type_identifiers[[i,1]]
  cell_type = cell_type_identifiers[[i,2]]
  for (j in 2:2364){
    cell_id = female_sle[[1,j]]
    if (cell_name == cell_id){
      xist_expression = female_sle[[42,j]]
      if (cell_type == "CE0"){
        epithelial_cell_levels = rbind(xist_expression, epithelial_cell_levels)
      }
      if (cell_type == "CD0"){
        dividing_cell_levels = rbind(xist_expression, dividing_cell_levels)
      }
      if (cell_type == "CM0" | cell_type == "CM2" | cell_type == "CM3" | cell_type == "CM1" | cell_type == "CM4"){
        myeloid_cell_levels = rbind(xist_expression, myeloid_cell_levels)
      }
      if (cell_type == "CT5b" | cell_type == "CT0a" | cell_type == "CT3b" | cell_type == "CT0b" | cell_type == "CT5a" | cell_type == "CT4" | cell_type == "CT1" | cell_type == "CT3a" | cell_type == "CT2" | cell_type == "CT6"){
        t_cell_levels = rbind(xist_expression, t_cell_levels)
      }
      if (cell_type == "CB2a" | cell_type == "CB0" | cell_type == "CB2b" | cell_type == "CB1" | cell_type == "CB3"){
        b_cell_levels =  rbind(xist_expression, b_cell_levels)
      }
    }
  }
}

###############################################################################
### Label which cells came from which patients
###############################################################################
key_with_disease_status = read.table(file.choose()) # Choose celseq_meta file
# Sort cells into individual patients: 
sixoeight = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
othreesix = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
othreeseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofouro = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofourone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofourtwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofourthree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofourfour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofourfive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofoursix = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofourseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofoureight = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofournine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofivetwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofivethree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofiveseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ofivenine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
osixone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
osixtwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oninefive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneofour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneofive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneoseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneoeight = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneonine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneoneone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneonethree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneonefour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneonefive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwofour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwosix = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwoseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threethreeone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threethreethree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threethreenine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefouro = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefourone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefourthree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefourfour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefoureight = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefournine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefiveo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threefivetwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threesixone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threesixtwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threesixfour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
threenineone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixoeight = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonine = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixoneo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixoneone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonetwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonethree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonefour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonefive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixonesix = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
sixoneseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
eightfourone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
eightsevenone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
eightseventwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
eightseventhree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
eightsevenfour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
eightsevenfive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ninesixone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
ninesixtwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneoneoneone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneoneonetwo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
oneoneonethree = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwoonefour = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwoonefive = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwoonesix = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
onetwooneseven = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
twoonineo = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
twooeightone = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
none = data.frame(cell_name=character(0), stringsAsFactors = FALSE)
### So what we want to do now is change this up so that we include all cells, or all B cells, or all T cells, or all M cells
for (i in 1:8297){
  if (key_with_disease_status[i,4] != "Epithelial"){
    
    if (key_with_disease_status[i,3] == "none"){
      none = rbind(none, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0608"){
      sixoeight = rbind(sixoeight, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0609"){
      sixonine = rbind(sixonine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0036"){
      othreesix = rbind(othreesix, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0040"){
      ofouro = rbind(ofouro, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0041"){
      ofourone = rbind(ofourone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0042"){
      ofourtwo = rbind(ofourtwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0043"){
      ofourthree = rbind(ofourthree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0044"){
      ofourfour = rbind(ofourfour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0045"){
      ofourfive = rbind(ofourfive, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0046"){
      ofoursix = rbind(ofoursix, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0047"){
      ofourseven = rbind(ofourseven, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0048"){
      ofoureight = rbind(ofoureight, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0049"){
      ofournine = rbind(ofournine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0052"){
      ofivetwo = rbind(ofivetwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0053"){
      ofivethree = rbind(ofivethree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0057"){
      ofiveseven = rbind(ofiveseven, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0059"){
      ofivenine = rbind(ofivenine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0061"){
      osixone = rbind(osixone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0062"){
      osixtwo = rbind(osixtwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0095"){
      oninefive = rbind(ofivenine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0104"){
      oneofour = rbind(oneofour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0105"){
      oneofive = rbind(oneofive, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0107"){
      oneoseven = rbind(oneoseven, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0108"){
      oneoeight = rbind(oneoeight, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0109"){
      oneonine = rbind(oneonine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0111"){
      oneoneone = rbind(oneoneone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0113"){
      oneonethree = rbind(oneonethree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0114"){
      oneonefour = rbind(oneonefour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0115"){
      oneonefive = rbind(oneonefive, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0124"){
      onetwofour = rbind(onetwofour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0126"){
      onetwosix = rbind(onetwosix, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0127"){
      onetwoseven = rbind(onetwoseven, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0331"){
      threethreeone = rbind(threethreeone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0333"){
      threethreethree = rbind(threethreethree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0339"){
      threethreenine = rbind(threethreenine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0340"){
      threefouro = rbind(threefouro, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0341"){
      threefourone = rbind(threefourone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0343"){
      threefourthree = rbind(threefourthree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0344"){
      threefourfour = rbind(threefourfour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0348"){
      threefoureight = rbind(threefoureight, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0349"){
      threefournine = rbind(threefournine, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0350"){
      threefiveo = rbind(threefiveo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0352"){
      threefivetwo = rbind(threefivetwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0361"){
      threesixone = rbind(threesixone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0362"){
      threesixtwo = rbind(threesixtwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0364"){
      threesixfour = rbind(threesixfour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0391"){
      threenineone = rbind(threenineone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0610"){
      sixoneo = rbind(sixoneo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0611"){
      sixoneone = rbind(sixoneone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0612"){
      sixonetwo = rbind(sixonetwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0613"){
      sixonethree = rbind(sixonethree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0614"){
      sixonefour = rbind(sixonefour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0615"){
      sixonefive = rbind(sixonefive, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0616"){
      sixonesix = rbind(sixonesix, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0617"){
      sixoneseven = rbind(sixoneseven, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0841"){
      eightfourone = rbind(eightfourone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0871"){
      eightsevenone = rbind(eightsevenone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0872"){
      eightseventwo = rbind(eightseventwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0873"){
      eightseventhree = rbind(eightseventhree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0874"){
      eightsevenfour = rbind(eightsevenfour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0875"){
      eightsevenfive = rbind(eightsevenfive, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0961"){
      ninesixone = rbind(ninesixone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-0962"){
      ninesixtwo = rbind(ninesixtwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1111"){
      oneoneoneone = rbind(oneoneoneone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1112"){
      oneoneonetwo = rbind(oneoneonetwo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1113"){
      oneoneonethree = rbind(oneoneonethree, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1214"){
      onetwoonefour = rbind(onetwoonefour, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1215"){
      onetwoonefive = rbind(onetwoonefive, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1216"){
      onetwoonesix = rbind(onetwoonesix, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-1217"){
      onetwooneseven = rbind(onetwooneseven, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-2090"){
      twoonineo = rbind(twoonineo, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
    if (key_with_disease_status[i,3] == "200-2081"){
      twooeightone = rbind(twooeightone, toString(key_with_disease_status[i,1]), stringsAsFactors = FALSE)
    }
  }
}
####
my_data <- read.table(file = "clipboard", 
                      sep = "\t", header=FALSE, stringsAsFactors = FALSE) # Copy in SLE Females tab of ampdata_ifngenes_xist
# Female lupus patients
organized_data_f_sle = matrix(0,45,2573)
ofouro_expression = matrix(0,45,1)
osixone_expression = matrix(0,45,12)
oninefive_expression = matrix(0,45,1)
oneofour_expression = matrix(0,45,1)
oneofive_expression = matrix(0,45,2)
oneoseven_expression = matrix(0,45,6)
oneonine_expression = matrix(0,45,94)
oneonefour_expression = matrix(0,45,12)
oneonefive_expression = matrix(0,45,9)
threefourone_expression = matrix(0,45,7)
threesixone_expression = matrix(0,45,292)
threesixtwo_expression = matrix(0,45,167)
threenineone_expression = matrix(0,45,211)
sixonetwo_expression = matrix(0,45,19)
eightfourone_expression = matrix(0,45,379)
eightseventhree_expression = matrix(0,45,210)
eightsevenfour_expression = matrix(0,45,279)
eightsevenfive_expression = matrix(0,45,182)
ninesixone_expression = matrix(0,45,230)
onetwoonesix_expression = matrix(0,45,149)
onetwooneseven_expression = matrix(0,45,72)
osixtwo_expression = matrix(0,45,5)
oneoeight_expression = matrix(0,45,81)
onetwofour_expression = matrix(0,45,170)
found_cell_counter = 0 
ofouro_counter = 0
osixone_counter = 0
oninefive_counter = 0
oneofour_counter = 0
oneofive_counter = 0
oneoseven_counter = 0
oneonine_counter = 0
oneonefour_counter = 0
oneonefive_counter = 0
threefourone_counter = 0
threesixone_counter = 0
threesixtwo_counter = 0
threenineone_counter = 0
sixonetwo_counter = 0
eightfourone_counter = 0
eightseventhree_counter = 0
eightsevenfour_counter = 0
eightsevenfive_counter = 0
ninesixone_counter = 0
onetwoonesix_counter = 0
onetwooneseven_counter = 0
osixtwo_counter = 0
oneoeight_counter = 0
onetwofour_counter = 0
for (i in 1:2840){
  print(i)
  if (grepl((my_data[1,i]), ofouro)){
    found_cell_counter = found_cell_counter + 1
    ofouro_counter = ofouro_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    ofouro_expression[,ofouro_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), osixone)){
    found_cell_counter = found_cell_counter + 1
    osixone_counter = osixone_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    osixone_expression[,osixone_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oninefive)){
    found_cell_counter = found_cell_counter + 1
    oninefive_counter = oninefive_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oninefive_expression[,oninefive_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneofour)){
    found_cell_counter = found_cell_counter + 1
    oneofour_counter = oneofour_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneofour_expression[,oneofour_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneofive)){
    found_cell_counter = found_cell_counter + 1
    oneofive_counter = oneofive_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneofive_expression[,oneofive_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneoseven)){
    found_cell_counter = found_cell_counter + 1
    oneoseven_counter = oneoseven_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneoseven_expression[,oneoseven_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneonine)){
    found_cell_counter = found_cell_counter + 1
    oneonine_counter = oneonine_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneonine_expression[,oneonine_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneonefour)){
    found_cell_counter = found_cell_counter + 1
    oneonefour_counter = oneonefour_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneonefour_expression[,oneonefour_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneonefive)){
    found_cell_counter = found_cell_counter + 1
    oneonefive_counter = oneonefive_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneonefive_expression[,oneonefive_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), threefourone)){
    found_cell_counter = found_cell_counter + 1
    threefourone_counter = threefourone_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    threefourone_expression[,threefourone_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), threesixone)){
    found_cell_counter = found_cell_counter + 1
    threesixone_counter = threesixone_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    threesixone_expression[,threesixone_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), threesixtwo)){
    found_cell_counter = found_cell_counter + 1
    threesixtwo_counter = threesixtwo_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    threesixtwo_expression[,threesixtwo_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), threenineone)){
    found_cell_counter = found_cell_counter + 1
    threenineone_counter = threenineone_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    threenineone_expression[,threenineone_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), sixonetwo)){
    found_cell_counter = found_cell_counter + 1
    sixonetwo_counter = sixonetwo_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    sixonetwo_expression[,sixonetwo_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), eightfourone)){
    found_cell_counter = found_cell_counter + 1
    eightfourone_counter = eightfourone_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    eightfourone_expression[,eightfourone_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), eightseventhree)){
    found_cell_counter = found_cell_counter + 1
    eightseventhree_counter = eightseventhree_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    eightseventhree_expression[,eightseventhree_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), eightsevenfour)){
    found_cell_counter = found_cell_counter + 1
    eightsevenfour_counter = eightsevenfour_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    eightsevenfour_expression[,eightsevenfour_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), eightsevenfive)){
    found_cell_counter = found_cell_counter + 1
    eightsevenfive_counter = eightsevenfive_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    eightsevenfive_expression[,eightsevenfive_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), ninesixone)){
    found_cell_counter = found_cell_counter + 1
    ninesixone_counter = ninesixone_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    ninesixone_expression[,ninesixone_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), onetwoonesix)){
    found_cell_counter = found_cell_counter + 1
    onetwoonesix_counter = onetwoonesix_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    onetwoonesix_expression[,onetwoonesix_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), onetwooneseven)){
    found_cell_counter = found_cell_counter + 1
    onetwooneseven_counter = onetwooneseven_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    onetwooneseven_expression[,onetwooneseven_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), osixtwo)){
    found_cell_counter = found_cell_counter + 1
    osixtwo_counter = osixtwo_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    osixtwo_expression[,osixtwo_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), oneoeight)){
    found_cell_counter = found_cell_counter + 1
    oneoeight_counter = oneoeight_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    oneoeight_expression[,oneoeight_counter] = my_data[,i]
  }
  if (grepl((my_data[1,i]), onetwofour)){
    found_cell_counter = found_cell_counter + 1
    onetwofour_counter = onetwofour_counter + 1
    organized_data_f_sle[,found_cell_counter] = my_data[,i]
    onetwofour_expression[,onetwofour_counter] = my_data[,i]
  }
}
XIST_in_onetwooneseven = onetwooneseven_expression[42,1:72]
XIST_in_onetwoonesix = onetwoonesix_expression[42,1:149]
XIST_in_ninesixone = ninesixone_expression[42,1:230]
XIST_in_eightsevenfive = eightsevenfive_expression[42,1:182]
XIST_in_eightsevenfour = eightsevenfour_expression[42,1:279]
XIST_in_eightseventhree = eightseventhree_expression[42,1:210]
XIST_in_eightfourone = eightfourone_expression[42,1:379]
XIST_in_sixonetwo = sixonetwo_expression[42,1:19]
XIST_in_threenineone = threenineone_expression[42,1:1] # only has 1
XIST_in_threesixtwo = threesixtwo_expression[42,1:167] 
XIST_in_threesixone = threesixone_expression[42,1:292]
XIST_in_threefourone = threefourone_expression[42,1:7]
XIST_in_oneonefive = oneonefive_expression[42,1:9]
XIST_in_oneonefour = oneonefour_expression[42,1:12]
XIST_in_oneonine = oneonine_expression[42,1:94]
XIST_in_oneoseven = oneoseven_expression[42,1:6]
XIST_in_oneofive = oneofive_expression[42,1:2] # only has 2
XIST_in_osixone = osixone_expression[42,1:12] 
XIST_in_ofouro = ofouro_expression[42,1:1]
XIST_in_osixtwo = osixtwo_expression[42,1:4]
XIST_in_oneoeight = oneoeight_expression[42,1:77]
XIST_in_onetwofour = onetwofour_expression[42,1:158]
ifn_sig_in_onetwooneseven = onetwooneseven_expression[2:41, 1:72]
ifn_sig_in_onetwoonesix = onetwoonesix_expression[2:41, 1:149]    # Should make me a matrix with all the ifn-responsive genes for all 1217 cells
ifn_sig_in_ninesixone = ninesixone_expression[2:41, 1:230]
ifn_sig_in_eightsevenfive = eightsevenfive_expression[2:41, 1:182]
ifn_sig_in_eightsevenfour = eightsevenfour_expression[2:41, 1:279]
ifn_sig_in_eightseventhree = eightseventhree_expression[2:41, 1:210]
ifn_sig_in_eightfourone = eightfourone_expression[2:41, 1:379]
ifn_sig_in_sixonetwo = sixonetwo_expression[2:41, 1:19]
ifn_sig_in_threenineone = threenineone_expression[2:41, 1:1]
ifn_sig_in_threesixtwo = threesixtwo_expression[2:41, 1:167]
ifn_sig_in_threesixone = threesixone_expression[2:41, 1:292]
ifn_sig_in_threefourone = threefourone_expression[2:41, 1:7]
ifn_sig_in_oneonefive = oneonefive_expression[2:41, 1:9]
ifn_sig_in_oneonefour = oneonefour_expression[2:41, 1:12]
ifn_sig_in_oneonine = oneonine_expression[2:41, 1:94]
ifn_sig_in_oneoseven = oneoseven_expression[2:41, 1:6]
ifn_sig_in_oneofive = oneofive_expression[2:41, 1:2]
ifn_sig_in_osixone = osixone_expression[2:41, 1:12]
ifn_sig_in_ofouro = ofouro_expression[2:41, 1:1]
ifn_sig_in_osixtwo = osixtwo_expression[2:41, 1:4]
ifn_sig_in_oneoeight = oneoeight_expression[2:41, 1:77]
ifn_sig_in_onetwofour = onetwofour_expression[2:41, 1:158]
#### Calculate IFN scores
ifn_score_in_onetwooneseven = sum(as.numeric(ifn_sig_in_onetwooneseven[,1:72]))/72
ifn_score_in_onetwoonesix = sum(as.numeric(ifn_sig_in_onetwoonesix[,1:149]))/149
ifn_score_in_ninesixone = sum(as.numeric(ifn_sig_in_ninesixone[,1:230]))/230
ifn_score_in_eightsevenfive = sum(as.numeric(ifn_sig_in_eightsevenfive[,1:182]))/182
ifn_score_in_eightsevenfour = sum(as.numeric(ifn_sig_in_eightsevenfour[,1:279]))/279
ifn_score_in_eightseventhree = sum(as.numeric(ifn_sig_in_eightseventhree[,1:210]))/210
ifn_score_in_eightfourone = sum(as.numeric(ifn_sig_in_eightfourone[,1:379]))/379
ifn_score_in_sixonetwo = sum(as.numeric(ifn_sig_in_sixonetwo[,1:19]))/19
ifn_score_in_threenineone = sum(as.numeric(ifn_sig_in_threenineone[,1:1]))/1
ifn_score_in_threesixtwo = sum(as.numeric(ifn_sig_in_threesixtwo[,1:167]))/167
ifn_score_in_threesixone = sum(as.numeric(ifn_sig_in_threesixone[,1:292]))/292
ifn_score_in_threefourone = sum(as.numeric(ifn_sig_in_threefourone[,1:7]))/7
ifn_score_in_oneonefive = sum(as.numeric(ifn_sig_in_oneonefive[,1:9]))/9
ifn_score_in_oneonefour = sum(as.numeric(ifn_sig_in_oneonefour[,1:12]))/12
ifn_score_in_oneonine = sum(as.numeric(ifn_sig_in_oneonine[,1:94]))/94
ifn_score_in_oneoseven = sum(as.numeric(ifn_sig_in_oneoseven[,1:6]))/6
ifn_score_in_oneofive = sum(as.numeric(ifn_sig_in_oneofive[,1:2]))/2
ifn_score_in_osixone = sum(as.numeric(ifn_sig_in_osixone[,1:12]))/12
ifn_score_in_ofouro = sum(as.numeric(ifn_sig_in_ofouro[,1:1]))/1
ifn_score_in_osixtwo = sum(as.numeric(ifn_sig_in_osixtwo[,1:4]))/4
ifn_score_in_oneoeight = sum(as.numeric(ifn_sig_in_oneoeight[,1:77]))/77
ifn_score_in_onetwofour = sum(as.numeric(ifn_sig_in_onetwofour[,1:158]))/158
### Calculate avg xist expression
avg_xist_in_onetwooneseven = mean(as.numeric(XIST_in_onetwooneseven))
avg_xist_in_onetwoonesix = mean(as.numeric(XIST_in_onetwoonesix))
avg_xist_in_ninesixone = mean(as.numeric(XIST_in_ninesixone))
avg_xist_in_eightsevenfive = mean(as.numeric(XIST_in_eightsevenfive))
avg_xist_in_eightsevenfour = mean(as.numeric(XIST_in_eightsevenfour))
avg_xist_in_eightseventhree = mean(as.numeric(XIST_in_eightseventhree))
avg_xist_in_eightfourone = mean(as.numeric(XIST_in_eightfourone))
avg_xist_in_sixonetwo = mean(as.numeric(XIST_in_sixonetwo))
avg_xist_in_threenineone = mean(as.numeric(XIST_in_threenineone))
avg_xist_in_threesixtwo = mean(as.numeric(XIST_in_threesixtwo))
avg_xist_in_threesixone = mean(as.numeric(XIST_in_threesixone))
avg_xist_in_threefourone = mean(as.numeric(XIST_in_threefourone))
avg_xist_in_oneonefive = mean(as.numeric(XIST_in_oneonefive))
avg_xist_in_oneonefour = mean(as.numeric(XIST_in_oneonefour))
avg_xist_in_oneonine = mean(as.numeric(XIST_in_oneonine))
avg_xist_in_oneoseven = mean(as.numeric(XIST_in_oneoseven))
avg_xist_in_oneofive = mean(as.numeric(XIST_in_oneofive))
avg_xist_in_osixone = mean(as.numeric(XIST_in_osixone))
avg_xist_in_ofouro = mean(as.numeric(XIST_in_ofouro))
avg_xist_in_osixtwo = mean(as.numeric(XIST_in_osixtwo))
avg_xist_in_oneoeight = mean(as.numeric(XIST_in_oneoeight))
avg_xist_in_onetwofour = mean(as.numeric(XIST_in_onetwofour))