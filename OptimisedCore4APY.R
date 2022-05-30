# "Optimisation of the core subset for the APY approximation of genomic relationships"
# Pocrnic, Lindgren, Tolhurst, Herring & Gorjanc (2022)
# DOI: 

# Clean the workspace
rm(list = ls())

# Load needed libraries
library("tidyverse")
library("AlphaSimR")
library("Matrix")

source("functions_opticore.R")

# Note: BLUPF90 has to be available (http://nce.ads.uga.edu/software/)
# Create BLUPF90 parameter files and bash scripts
prepare_par()

#### Simulation ####

founderPop = runMacs(nInd = 3000,
                     nChr = 10,
                     segSites = 1100,
                     species = "CATTLE")

SP = SimParam$new(founderPop)

SP$addTraitA(100)
SP$setVarE(h2 = 0.30)
SP$addSnpChip(1000)
SP$setSexes("yes_sys")

Parents = newPop(founderPop)

# Select randomly 1500 Females and 50 Males from the founder population
Sires = selectInd(Parents, 50, use = "rand", sex = "M")
Dams = selectInd(Parents, 1500, use = "rand", sex = "F")

# Record Data
Generation = 0
phenotypes = NULL
phenotypes = data_rec(phenotypes, c(Sires, Dams))

# List of Genotypes
rm(genoM)
genoM = list()

for(Generation in 1:20){
  # Generation = 1
  cat("Working on the Generation:", Generation, "\n")

  # Mate 50 Sires and 1500 Dams, and get 3000 Progeny
  pop = randCross2(males = Sires, females = Dams, nCrosses = 1500, nProgeny = 2)
  
  # Record Data
  phenotypes = data_rec(phenotypes, pop)
  
  # Export Genotypes for the last 5 generations
  if (Generation > 15) {
    genoM[[Generation]] = pullSnpGeno(pop)
  }
  
  # Select the best 50% (750/1500) females as replacement dams
  YoungDams = selectInd(pop, 750, use = "pheno", sex = "F")
  # Replacement rate for females is 50% (750/1500 removed)
  OldDams = selectInd(Dams, 750, use = "pheno", sex = "F")
  # Dams for the next generation
  Dams = c(OldDams, YoungDams)
  
  # Select the best 3% (45/1500) males as replacement sires
  YoungSires = selectInd(pop, 45, use = "gv", sex = "M")
  # Replacement rate for males is 90% (45/50 removed)
  OldSires = selectInd(Sires, 5, use = "gv", sex = "M")
  # Sires for the next generation
  Sires = c(OldSires, YoungSires)

}

# Remove phenotypes from the last generation (validation)
phenotypes = phenotypes %>%
  dplyr::mutate(Trait = replace(Trait, Generation == 20, NA))

# Separately add 500 individuals from the generation 21
# This individuals will be used for testing the updating algorithm
Generation = Generation + 1
pop = randCross2(males = Sires, females = Dams, nCrosses = 500, nProgeny = 1, balance = FALSE)
genoM_update = pullSnpGeno(pop)
# Create genotype matrix for updatinh
p_u = colMeans(genoM_update) / 2
P_u = matrix(p_u, nrow = nrow(genoM_update), ncol = ncol(genoM_update), byrow = TRUE)
Z_update = genoM_update - 2 * P_u
dim(Z_update)

# Create initial genotype matrix
M = do.call("rbind", genoM)
p_i = colMeans(M) / 2
P_i = matrix(p_i, nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
Z = M - 2 * P_i
dim(Z)

#### Singular Value Decomposition #### 

# Run SVD to find the number of eigenvalues that explain certain percentage of variance in Z
start_time = Sys.time()
  svd_Z = svd(Z)
end_time = Sys.time()
time_svd = end_time - start_time

variance.explained = prop.table(svd_Z$d^2)
variance_explained_table = tibble(eigen_number = seq(1:length(svd_Z$d)), eigen_value = svd_Z$d^2, var_explained = prop.table(svd_Z$d^2), cumvar = cumsum(var_explained) )

# Determining number of dimensions that explain ~n% of total variance
n_final = NULL
for (v in c(0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.98, 0.99)) {
  n_dims = 0
  total_variance = 0
  for (i in 1:length(variance_explained_table$var_explained)) {
    total_variance = total_variance + variance_explained_table$var_explained[i]
    ifelse(total_variance <= v, n_dims <- n_dims + 1, break)
  }
  n_final = cbind(n_final, n_dims)
}
colnames(n_final) = c("e10", "e30", "e50", "e70", "e90", "e95", "e98", "e99")

write.table(n_final, "coresize.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ", na = "0")


#### Initialize GRM ####

zz = tcrossprod(Z)
k = 2 * (sum(p_i * (1 - p_i)))
G = zz / k
dim(G)

rm(zz)

# Add to diagonal
nid = nrow(G)
# Could be done without the 0.99*G part 
fG = 0.99*G + (diag(0.01, nrow(G)))

# Create a reference file for genotyped animals and mark the core/noncore set:
geno_ref = tibble(geno_id = rownames(M), geno_order = seq(1:length(geno_id)), core_status = rep(NA, times = length(geno_id)))
geno_ref_backup = geno_ref


#### A) Direct (Full) Inverse ####
start_time = Sys.time()
 Ginv_full = solve(fG)
end_time = Sys.time()
time_inv_full = end_time - start_time

checksym = isSymmetric(Ginv_full)

if (checksym == "FALSE") {
  Ginv_full2 = forceSymmetric(Ginv_full)
  # isSymmetric(Ginv_full2)
  G_out1 = as(Ginv_full2, "dsTMatrix")
} else {
  G_out1 = as(Ginv_full, "dsTMatrix")
}

# G_out1 = as(Ginv_full2, "dgTMatrix") -> prints all trios with repeated elements 
G_out2 = summary(G_out1)
write.table(G_out2, file = "apyinv.txt", row.names = F, col.names = F)

# Clean intermediate files
rm(Ginv_full2, G_out1, G_out2)

# Export inverse created in R, to use as external file in blupf90
# use option "user_file" in blupf90

# Create datafile for BLUPF90
pheno_export = inner_join(phenotypes, geno_ref, by = c("Aid" = "geno_id")) %>%
  arrange(geno_order)

# Export data
write.table(cbind(pheno_export[,c("Aid", "Trait", "geno_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# Run blupf90
start_time = Sys.time()
 tmp2 = run_gblupf90(pheno_export, "Full", "Regular")
end_time = Sys.time()
time_gblup_full = end_time - start_time

final_data = inner_join(phenotypes, tmp2, by = "Aid")


#### B) APY with random core selection ####

# Set the number of core animals using the loop
rm(final_data2)
final_data2 = final_data
rm(geno_refAPY)
geno_refAPY = geno_ref[,-3]

start_time = Sys.time()
for(j in 1:(length(n_final))){

 varex = colnames(n_final)[j]
 ncore = n_final[j]
 
 # Run 5 replicates of random core sampling, for each number of core animals 
  for(apy_rep in 1:5){
    # Select random animals as core 
    index_core = sample_n(geno_ref[2], size = ncore)
  
    geno_ref2 = geno_ref %>%
      dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))
  
    # Run APY inverse function
    apyginvlist = APY_inverse(fG, geno_ref2)
  
    checksym = isSymmetric(apyginvlist$APY_Ginv)
  
    if (checksym == "FALSE") {
      Ginv_full2 = forceSymmetric(apyginvlist$APY_Ginv)
      # isSymmetric(Ginv_full2)
      G_out1 = as(Ginv_full2, "dsTMatrix")
    } else {
      G_out1 = as(apyginvlist$APY_Ginv, "dsTMatrix")
    }
  
    G_out2 = summary(G_out1)
    write.table(G_out2, file = "apyinv.txt", row.names = F, col.names = F)

    # Clean intermediate files
    rm(Ginv_full2, G_out1, G_out2)
   
    # Create datafile for BLUPF90
    pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
      arrange(apy_order)
    
    # Export data:
    write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
    # Run blupf90
    ebvcolname = paste("APY_", varex, "_Rep", apy_rep, sep="")
  
    tmp2 = run_gblupf90(pheno_export, ebvcolname, "APY")
  
    final_data2 = inner_join(final_data2, tmp2, by = "Aid")
    
    # Save "geno_ref2" from each replicate - list of core animals selected each time
    colnames(geno_ref2) = c("geno_id", "geno_order" , ebvcolname)
    geno_refAPY = inner_join(geno_refAPY, geno_ref2, by = "geno_id")
    geno_refAPY = dplyr::select(geno_refAPY, -starts_with("geno_order"))
    
  }

}
end_time = Sys.time()
time_APY_random = end_time - start_time

rm(apyginvlist)


#### C) APY with conditional core selection ####

rm(final_data3)
final_data3 = final_data2
rm(geno_refPrior)
geno_refPrior = geno_ref[,-3]

A <- Z

start_time = Sys.time()
for(j in 1:(length(n_final))){
  
  varex = colnames(n_final)[j]
  ncore = n_final[j]

knots_prior <- find_knots(A, p = ncore, method = "prior")

opticore <- knots_prior$data

index_core = NULL


index_core$geno_order = opticore$knots

geno_ref2 = geno_ref %>%
  dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))

# Run APY inverse function
apyginvlist = APY_inverse(fG, geno_ref2)

checksym = isSymmetric(apyginvlist$APY_Ginv)

if (checksym == "FALSE") {
  Ginv_full2 = forceSymmetric(apyginvlist$APY_Ginv)
  # isSymmetric(Ginv_full2)
  G_out1 = as(Ginv_full2, "dsTMatrix")
} else {
  G_out1 = as(apyginvlist$APY_Ginv, "dsTMatrix")
}

G_out2 = summary(G_out1)
write.table(G_out2, file = "apyinv.txt", row.names = F, col.names = F)

# Clean intermediate files
rm(Ginv_full2, G_out1, G_out2)

# Create datafile for BLUPF90
pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
  arrange(apy_order)

# Export data:
write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# Run blupf90

ebvcolname = paste("Prior_", varex, sep="")

tmp2 = run_gblupf90(pheno_export, ebvcolname, "APY")

final_data3 = inner_join(final_data3, tmp2, by = "Aid")

# Save "geno_ref2" from each replicate - list of core animals selected each time
colnames(geno_ref2) = c("geno_id", "geno_order" , ebvcolname)
geno_refPrior = inner_join(geno_refPrior, geno_ref2, by = "geno_id")
geno_refPrior = dplyr::select(geno_refPrior, -starts_with("geno_order"))

}
end_time = Sys.time()
time_prior = end_time - start_time

rm(apyginvlist, A)


#### D) APY with diagonal core selection ####

rm(final_data4)
final_data4 = final_data3
rm(geno_refNaive1)
geno_refNaive1 = geno_ref[,-3]

# Get diagonal of G (either form G or Z) 
geno_ref_diag = geno_ref
geno_ref_diag$v_prior <- rowSums(Z^2)

start_time = Sys.time()
for(j in 1:(length(n_final))){
  
  varex = colnames(n_final)[j]
  ncore = n_final[j]
  
  index_core = dplyr::top_n(geno_ref_diag, n = ncore, wt = geno_ref_diag[4])
  # Top n diagonals as core
  
  geno_ref2 = geno_ref %>%
    dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))
  
  # Run APY inverse function
  apyginvlist = APY_inverse(fG, geno_ref2)
  
  checksym = isSymmetric(apyginvlist$APY_Ginv)
  
  if (checksym == "FALSE") {
    Ginv_full2 = forceSymmetric(apyginvlist$APY_Ginv)
    # isSymmetric(Ginv_full2)
    G_out1 = as(Ginv_full2, "dsTMatrix")
  } else {
    G_out1 = as(apyginvlist$APY_Ginv, "dsTMatrix")
  }
  
  G_out2 = summary(G_out1)
  write.table(G_out2, file = "apyinv.txt", row.names = F, col.names = F)

  # Clean intermediate files
  rm(Ginv_full2, G_out1, G_out2)
  
  # Create datafile for BLUPF90
  pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
    arrange(apy_order)
  
  # Export data:
  write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
  # Run blupf90
  
  ebvcolname = paste("Diag_", varex, sep="")
  
  tmp2 = run_gblupf90(pheno_export, ebvcolname, "APY")
  
  final_data4 = inner_join(final_data4, tmp2, by = "Aid")
  
  # Save "geno_ref2" from each replicate - list of core animals selected each time
  colnames(geno_ref2) = c("geno_id", "geno_order" , ebvcolname)
  geno_refNaive1 = inner_join(geno_refNaive1, geno_ref2, by = "geno_id")
  geno_refNaive1 = dplyr::select(geno_refNaive1, -starts_with("geno_order"))
  
}
end_time = Sys.time()
time_naive1 = end_time - start_time

rm(apyginvlist)

  
#### E) APY with weighted diagonal core selection ####
  
rm(final_data5)
final_data5 = final_data4
rm(geno_refNaive2)
geno_refNaive2 = geno_ref[,-3]

start_time = Sys.time()
for(j in 1:(length(n_final))){
  
  varex = colnames(n_final)[j]
  ncore = n_final[j]
  
  # Run 5 replicates of random core sampling, for each number of core animals 
  for(apy_rep in 1:5){
    # Select random animals as core 
    index_core = sample_n(geno_ref_diag, size = ncore, weight = unlist(geno_ref_diag[4]))
    # Random n with diagonals as weights 
    # Weights; prop.table(geno_ref[4])
    # Alternative: index_core = sample(nrow(geno_ref), ncore, prob = unlist(geno_ref[4]))
    # Alternative: index_core = sample(geno_ref$geno_order, ncore, prob = unlist(geno_ref[4]))
  
    geno_ref2 = geno_ref %>%
      dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))
    
    # Run APY inverse function
    apyginvlist = APY_inverse(fG, geno_ref2)
    
    checksym = isSymmetric(apyginvlist$APY_Ginv)
    
    if (checksym == "FALSE") {
      Ginv_full2 = forceSymmetric(apyginvlist$APY_Ginv)
      # isSymmetric(Ginv_full2)
      G_out1 = as(Ginv_full2, "dsTMatrix")
    } else {
      G_out1 = as(apyginvlist$APY_Ginv, "dsTMatrix")
    }
    
    G_out2 = summary(G_out1)
    write.table(G_out2, file = "apyinv.txt", row.names = F, col.names = F)
    
    # Clean intermediate files
    rm(Ginv_full2, G_out1, G_out2)
    
    # Create datafile for BLUPF90
    pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
      arrange(apy_order)
    
    # Export data:
    write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
    
    # Run blupf90
    
    ebvcolname = paste("Weight_", varex, "_Rep", apy_rep, sep="")
    
    tmp2 = run_gblupf90(pheno_export, ebvcolname, "APY")
    
    final_data5 = inner_join(final_data5, tmp2, by = "Aid")
    
    # Save "geno_ref2" from each replicate - list of core animals selected each time
    colnames(geno_ref2) = c("geno_id", "geno_order" , ebvcolname)
    geno_refNaive2 = inner_join(geno_refNaive2, geno_ref2, by = "geno_id")
    geno_refNaive2 = dplyr::select(geno_refNaive2, -starts_with("geno_order"))
    
  }
  
}
end_time = Sys.time()
time_naive2 = end_time - start_time

rm(apyginvlist)


#### F) APY with conditional core selection - using reduced rank matrix ####

rm(final_data6)
final_data6 = final_data5
rm(geno_refPriorSVD)
geno_refPriorSVD = geno_ref[,-3]

# Example 99% variance = 3130
# Note; before A = Z

A <- svd_Z$u[,1:3130] * rep(svd_Z$d[1:3130], each = nrow(svd_Z$u))
dim(A)

start_time = Sys.time()
for(j in 1:(length(n_final))){
  
  varex = colnames(n_final)[j]
  ncore = n_final[j]
  
  knots_prior <- find_knots(A, p = ncore, method = "prior")
  
  opticore <- knots_prior$data
  
  index_core = NULL
  
  index_core$geno_order = opticore$knots

geno_ref2 = geno_ref %>%
  dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))

# Run APY inverse function
apyginvlist = APY_inverse(fG, geno_ref2)

checksym = isSymmetric(apyginvlist$APY_Ginv)

if (checksym == "FALSE") {
  Ginv_full2 = forceSymmetric(apyginvlist$APY_Ginv)
  # isSymmetric(Ginv_full2)
  G_out1 = as(Ginv_full2, "dsTMatrix")
} else {
  G_out1 = as(apyginvlist$APY_Ginv, "dsTMatrix")
}

G_out2 = summary(G_out1)
write.table(G_out2, file = "apyinv.txt", row.names = F, col.names = F)

# Clean intermediate files
rm(Ginv_full2, G_out1, G_out2)

# Create datafile for BLUPF90
pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
  arrange(apy_order)

# Export data:
write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# Run blupf90

ebvcolname = paste("SVDprior_", varex, sep="")

tmp2 = run_gblupf90(pheno_export, ebvcolname, "APY")

final_data6 = inner_join(final_data6, tmp2, by = "Aid")

# Save "geno_ref2" from each replicate - list of core animals selected each time
colnames(geno_ref2) = c("geno_id", "geno_order" , ebvcolname)
geno_refPriorSVD = inner_join(geno_refPriorSVD, geno_ref2, by = "geno_id")
geno_refPriorSVD = dplyr::select(geno_refPriorSVD, -starts_with("geno_order"))

}
end_time = Sys.time()
time_priorSVD = end_time - start_time

rm(apyginvlist, A)


#### G) Updating with arrival of new data ####

# First set of core individuals (10% var = 10 animals), prior method
A1 = Z
A_new = Z_update

knots_prior <- find_knots(A1, p = 10, method = "prior")

# Find another 5 knots in rbind(A1, A_new):
info2 <- find_knots(
     A = rbind(A1, A_new), # This A won't actually be used
     knots = knots_prior$data$knots,
     B = rbind(knots_prior$B,
               expand_B(A1, knots_prior$data$knots, A_new)),
     p = 5,
     method = "prior")

knots_prior_15 <- find_knots(rbind(A1, A_new), p = 15, method = "prior")

knots_prior$data$knots
info2$data$knots
knots_prior_15$data$knots

table(info2$data$knots %in% knots_prior$data$knots)
table(info2$data$knots %in% knots_prior_15$data$knots)

  