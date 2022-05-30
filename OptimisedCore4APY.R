# "Optimisation of the core subset for the APY approximation of genomic relationships"
# Pocrnic, Lindgren, Tolhust, Herring & Gorjanc (2022)
# DOI: 

rm(list = ls())

# Load R libraries
library("tidyverse")
library("AlphaSimR")
library("Matrix")
library("umap")
library("gridExtra")

setwd("~/Desktop/test_optiapy/")
getwd()

source("functions_opticore.R")

# Note: BLUPF90 has to be available (http://nce.ads.uga.edu/software/)

# Create BLUPF90 parameter files and bash scripts
prepare_par()

#### Simulation ####

founderPop = runMacs(nInd = 3000,
                     nChr = 10,
                     segSites = 1100,
                     species = "CATTLE")

# save.image("founders.RData")
# load("founders.RData")

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

#### SVD and PCA #### 
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


# Run PCA (Redundant since SVD is already saved) 

start_time = Sys.time()
 # Note: Default in prcomp is center = TRUE, but our Z matrix is already centered
 res.pca = prcomp(Z, center = FALSE)
end_time = Sys.time()
time_pca = end_time - start_time

ind.coord <- res.pca$x
# head(ind.coord[, 1:2])
# tail(ind.coord[, 1:2])
# aa = summary(res.pca)
# First 10 PC's (Var, Cumulative)
# aa$importance[2, 1:10]
# aa$importance[3, 1:10]
# aa$importance[3, 2505]
# pc_1 = round(aa$importance[2, 1:1]*100, 2)
# pc_2 = round(aa$importance[2, 2:2]*100, 2)
# rm(aa)

# Plot PC1 vs. PC2 
gen_no = phenotypes %>%
  filter(Generation > 15) 

ind.coord_2 = data.frame(ind.coord[, 1:2], factor(gen_no$Generation))

p_pca = ggplot(ind.coord_2, aes(x = ind.coord_2[, 1], y = ind.coord_2[, 2], group = ind.coord_2[, 3])) + geom_point(aes(color = ind.coord_2[, 3], alpha = ind.coord_2[, 3])) + 
  theme_classic() + labs(y = "PC 2 (1.25%)", x = "PC 1 (1.18%)", colour = "Generation", alpha = "Generation", tag = "A") + scale_color_manual(values = color_roslin)

# Raw PCA plot (no alpha)
p_pca2 = ggplot(ind.coord_2, aes(x = ind.coord_2[, 1], y = ind.coord_2[, 2], group = ind.coord_2[, 3])) + geom_point(aes(color = ind.coord_2[, 3])) + 
  labs(x = paste("PC 1 (1.3 %)"), y = paste("PC 2 (1.2 %)"), colour = "Generation") +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  scale_color_manual(values = color_roslin) + theme(legend.position = "top") 

#### UMAP Visualization ####

custom.config = umap.defaults
# custom.config$n_neighbors = 5
# custom.configs$min_dist = 0.5
# NN and MD are two main parameters to play with in UMAP

start_time = Sys.time()
 res.umap = umap(d = Z, config = custom.config)
end_time = Sys.time()
time_umap = end_time - start_time

# head(res.umap$layout, 5)
ind.coord_3 = data.frame(res.umap$layout, factor(gen_no$Generation), row.names(res.umap$layout))
# With outliers 
p_umap = ggplot(ind.coord_3, aes(x = ind.coord_3[, 1], y = ind.coord_3[, 2], group = ind.coord_3[, 3])) + 
  geom_point(aes(color = ind.coord_3[, 3], alpha = ind.coord_3[, 3])) + 
  theme_classic() + labs(y = "", x = "", colour = "Generation", alpha = "Generation", tag = "B") + scale_color_manual(values = color_roslin) + theme(legend.position = "top")

ggsave(plot = p_umap + PaperTheme +  
         theme(
           axis.text.x = element_blank(), 
           axis.text.y = element_blank(), 
           axis.ticks = element_blank()) , 
       filename = "UMAP_outliers.png", height = PaperSize, width = PaperSize * 1.5, unit = "cm")

# Take note:
# Issue are the outliers in the first genotyped generation, coming from two sires only
# GG: It is likely some odd/rare combination of base population haplotypes in some base population individuals. 
# Here is a nice example, how base pop diversity can be soon lost.

# Remove outliers 
ind.coord_33 = inner_join(ind.coord_3, phenotypes, by = c("row.names.res.umap.layout." = "Aid") )
# 
ind.coord_333 = ind.coord_33 %>%
  filter(Sid!=46859 & Sid!=47751)

length(unique(ind.coord_333$Sid))

p_umap = ggplot(ind.coord_333, aes(x = ind.coord_333[, 1], y = ind.coord_333[, 2], group = ind.coord_333[, 3])) + 
  geom_point(aes(color = ind.coord_333[, 3], alpha = ind.coord_333[, 3])) + 
  theme_classic() + labs(y = "", x = "", colour = "Generation", alpha = "Generation", tag = "B") + scale_color_manual(values = color_roslin) + theme(legend.position = "top")

ggsave(plot = p_umap + PaperTheme +  
         theme(
           axis.text.x = element_blank(), 
           axis.text.y = element_blank(), 
           axis.ticks = element_blank()) , 
       filename = "UMAP.png", height = PaperSize, width = PaperSize * 1.5, unit = "cm")


# Without alpha
p_umap2 = ggplot(ind.coord_333, aes(x = ind.coord_333[, 1], y = ind.coord_333[, 2])) +
  geom_point(aes(color = as.factor(ind.coord_333[, 3]), group = as.factor(ind.coord_333[, 3]), shape = as.factor(ind.coord_333[, 3]))) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + labs(y = "", x = "", colour = "Generation", shape = "Generation") + 
  scale_color_manual(values = color_roslin) + 
  scale_shape_manual(values = c(16, 16, 16, 16, 16)) +
  theme(legend.position = "top") 

# Combine PCA and UMAP

##### PCA & UMAP panel plot ##### 
p1a = p_pca2 + labs(tag = "a)") 

p2a = p_umap2 + labs(y = "", x = "", tag = "b)") 

get_legend = function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend = get_legend(p1a)

panel_2 = grid.arrange(p1a + theme(legend.position="none"), 
                       p2a + theme(legend.position="none"),
                       legend,
                       nrow = 2, ncol = 2,
                       layout_matrix = rbind(c(3,3), c(1,2)),
                       heights = c(1, 10),
                       widths = c(10, 10) )

ggsave(file = "Panel_2_sim.png", panel_2, height = PaperSize, width = PaperSize * 1.5, unit = "cm")

# UMAP per generation (Manuscript additional file)
# Start, change for each generation; i = 16..20
# Do the loop in final version

min_x = min(ind.coord_333[, 1])
max_x = max(ind.coord_333[, 1])
min_y = min(ind.coord_333[, 2])
max_y = max(ind.coord_333[, 2])

boja = 0
for (i in 16:20) {
  ind.coord_4 = ind.coord_333 %>%
    filter(ind.coord_333[,3]==i)
  
  boja = boja + 1
  
  pp = ggplot(ind.coord_4, aes(x = ind.coord_4[, 1], y = ind.coord_4[, 2], group = ind.coord_4[, 3])) + 
    geom_point(aes(color = ind.coord_4[, 3]), alpha = ind.coord_4[, 3]) + 
    theme_classic() + labs(y = "", x = "", colour = "Generation", alpha = "Generation") + scale_color_manual(values = color_roslin[boja]) +
    ylim(min_y, max_y) + xlim(min_x, max_x) 
  
  ggsave(plot = pp + PaperTheme +
           theme(
             axis.text.x = element_blank(), 
             axis.text.y = element_blank(), 
             axis.ticks = element_blank()) ,
         filename = paste("UMAP_gen", i, ".png", sep=""),
         height = PaperSize, width = PaperSize * 1.5, unit = "cm")
  
  # pp_gen = paste("UMAP_gen", i, sep="")
  # assign(pp_gen, pp)
}

# Panel per generation
panel_umapgen = grid.arrange(UMAP_gen16 + PaperTheme,
                             UMAP_gen17 + PaperTheme, 
                             UMAP_gen18 + PaperTheme, 
                             UMAP_gen19 + PaperTheme, 
                             UMAP_gen20 + PaperTheme, 
                             nrow = 3, ncol=2,
                             layout_matrix = rbind(c(1,2), c(3,4), c(5,6)),
                             heights = c(10, 10, 10),
                             widths = c(10, 10)
)

ggsave(file = "UMAP_GEN.png", panel_umapgen, height = PaperSize, width = PaperSize * 1.5, unit = "cm")



# UMAP is clearly grouping individuals into Paternal Half-Sibs Clusters
ind.coord_4 = ind.coord_333 %>%
  filter(ind.coord_333[,3]==20)

length(unique(ind.coord_4$Sid))
table(ind.coord_4$Sid)

# Some code to generate many distinctive colors
# From: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# library("RColorBrewer")
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# Alternative
# library("viridis")
# col_vector = viridis_pal(option = "D")(50) 

pp = ggplot(ind.coord_4, aes(x = ind.coord_4[, 1], y = ind.coord_4[, 2], shape = ind.coord_4[, 5])) + geom_point(aes(color = ind.coord_4[, 5])) + 
  theme_classic() + labs(y = "", x = "", colour = "Sire", shape = "Sire") + scale_color_manual(values = col_vector) + ylim(min_y, max_y) + xlim(min_x, max_x) + 
  scale_shape_manual(values=c(0:25,48:57,65:78))
pp + theme(legend.position="none")

ggsave(plot = pp + PaperTheme +
         theme(legend.position="none", 
               axis.text.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks = element_blank()) ,
       filename = "UMAP_gen20_HS.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm")



# Clean the things
rm(ind.coord, ind.coord_2, ind.coord_3, ind.coord_33, ind.coord_333, ind.coord_4)
gc()
# Save after SVD, PCA, and UMAP
save.image("after_svd_umap.RData")
# load("after_svd_umap.RData")


#### Initialize GRM ####
# This is "G" - Genomic (Covariance) Relationship Matrix:
zz = tcrossprod(Z)
# zz = Z %*% t(Z)
k = 2 * (sum(p_i * (1 - p_i)))
G = zz / k
dim(G)
format(object.size(G), units = "GB", standard = "SI")
object_size(G)

rm(zz)

# Add to diagonal:
nid = nrow(G)
fG = 0.99*G + (diag(0.01, nrow(G)))

# Create a reference file for genotyped animals and mark the core/noncore set:
geno_ref = tibble(geno_id = rownames(M), geno_order = seq(1:length(geno_id)), core_status = rep(NA, times = length(geno_id)))
geno_ref_backup = geno_ref

# Calculate true PA and MS
phenotypes_backup = phenotypes

phenotypes$TBV_sire = 0
phenotypes$TBV_dam = 0
phenotypes$TMS = 0
for (i in 1:nrow(phenotypes)){
  if(phenotypes$Sid[i] > 0 | phenotypes$Did[i] > 0)
  { 
    phenotypes$TBV_sire[i] = phenotypes$TBV[phenotypes$Aid==phenotypes$Sid[i] ]
    phenotypes$TBV_dam[i] = phenotypes$TBV[phenotypes$Aid==phenotypes$Did[i] ]
    phenotypes$TMS[i] = phenotypes$TBV[i] - 0.5 * (phenotypes$TBV_sire[i] + phenotypes$TBV_dam[i])
  } 
  
}


#### A) Direct Inverse ####
start_time = Sys.time()
 Ginv_full = solve(fG)
end_time = Sys.time()
time_inv_full = end_time - start_time
# ~ 1.5h on Mac
format(object.size(Ginv_full), units = "GB", standard = "SI")
object_size(Ginv_full)

# tol = 100 * .Machine$double.eps
# isSymmetric(Ginv_full, tol = 1e-12)

# Potentially faster alternative (and will be symmetric)
# Ginv_full_alt = chol2inv(chol(fG))
# max(abs(Ginv_full - Ginv_full_alt))

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

save.image("after_inv_full.RData")
# load("after_inv_full.RData")


# Create datafile for BLUPF90

pheno_export = inner_join(phenotypes, geno_ref, by = c("Aid" = "geno_id")) %>%
  arrange(geno_order)

# Export data:
write.table(cbind(pheno_export[,c("Aid", "Trait", "geno_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# It is simulated data, so we don't run renumf90
# Run blupf90
# The question here is do we really need blupf90 at all, or should we instead move solving MME to R
# I talked with Chris, and unfortunately AlphaMME does not support pre-created inverse as input

start_time = Sys.time()
 tmp2 = run_gblupf90(pheno_export, "Full", "Regular")
# Warning message (once per session) - run again, check details
end_time = Sys.time()
time_gblup_full = end_time - start_time
# ~ 5min on Mac

final_data = inner_join(phenotypes, tmp2, by = "Aid")


#### B) APY with random core selection ####

# Set the number of core animals using the loop

rm(final_data2)
final_data2 = final_data
rm(geno_refAPY)
geno_refAPY = geno_ref[,-3]

# n_final = read.table(file = "coresize.txt", header = TRUE)

start_time = Sys.time()
for(j in 1:(length(n_final))){

 varex = colnames(n_final)[j]
 ncore = n_final[j]
 
 # Run 5 replicates of random core sampling, for each number of core animals 
  for(apy_rep in 1:5){
    # Select random animals as core 
    # set.seed(12345)
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
  
    # Export the APY inverse created in R, to use as external file in blupf90
    # use option "user_file" in blupf90
    # This could be avoided, but APY option is not available in the public version
    
    # Clean intermediate files
    rm(Ginv_full2, G_out1, G_out2)
   
    # Create datafile for BLUPF90
    pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
      arrange(apy_order)
  
    # Export data:
    write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
    # It is simulated data, so we don't run renumf90
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
# ~ 5.5h on Mac

rm(apyginvlist)

save.image("after_apy.RData")
# load("after_apy.RData")

#### C1) Optimized core selection - Prior ####

rm(final_data3)
final_data3 = final_data2
rm(geno_refPrior)
geno_refPrior = geno_ref[,-3]

A <- Z

# n_final = read.table(file = "coresize.txt", header = TRUE)
# j = 8
# varex = colnames(n_final)[j]
# ncore = as.numeric(n_final[j])

start_time = Sys.time()
for(j in 1:(length(n_final))){
  
  varex = colnames(n_final)[j]
  ncore = n_final[j]

knots_prior <- find_knots(A, p = ncore, method = "prior")

opticore <- knots_prior$data

index_core = NULL

# This was wrong? has to be <= nrow(A)
# index_core$geno_order = as.numeric(row.names(opticore))
# row.names(opticore) -> original ID of animals, original row names in Z

index_core$geno_order = opticore$knots

geno_ref2 = geno_ref %>%
  dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))


  # dplyr::mutate(order_selected = rep(NA, times = length(geno_id))) %>%

  # dplyr::mutate(order_selected = if_else(geno_order %in% opticore$knots, true = opticore$k, false = 0))

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

# Export the APY inverse created in R, to use as external file in blupf90
# use option "user_file" in blupf90
# This could be avoided, but APY option is not available in the public version

# Clean intermediate files
rm(Ginv_full2, G_out1, G_out2)

# Create datafile for BLUPF90
pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
  arrange(apy_order)

# Export data:
write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# It is simulated data, so we don't run renumf90
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
# ~ 5h on Mac

rm(apyginvlist, A)

save.image("after_prior.RData")
# load("after_prior.RData")

# DataSet1 <- new.env() 
# load("after_apy.RData", envir = DataSet1)
# ls(DataSet1)
# head(DataSet1$final_data2)
# tmp2APY = DataSet1$final_data2 %>%
#  dplyr::select(starts_with("APY_e99"), starts_with("SE_APY_e99"), "Aid")
# final_data3 = inner_join(final_data3, tmp2APY, by = "Aid")
# geno_refAPY = DataSet1$geno_refAPY
# rm(DataSet1)



#### SKIP C2) Optimized core selection - Max or Min Posterior ####

# Not working!!! 
# Computational issues (look into test_functions.R)

# Note: Max algorithm is not useful anyways 
# profvis({
#   knots_max <- find_knots(A, p = p, method = "posterior_max")
# })

### Knots average 
# Warning: Run on small example to compare core selection vs. prior method

knots_mean <- find_knots(A, p = ncore, method = "posterior_avg")

opticore_mean <- knots_mean$data

index_core = NULL

index_core$geno_order = opticore_mean$knots

geno_ref_compare = geno_ref2 %>%
  dplyr::mutate(core_status_mean = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))

geno_ref_compare_2 = geno_ref_compare %>%
  dplyr::mutate(core_status_together = NA) %>%
  dplyr::mutate(core_status_together = replace(core_status_together, core_status == 1 & core_status_mean == 0, "PRIOR")) %>%
  dplyr::mutate(core_status_together = replace(core_status_together, core_status == 0 & core_status_mean == 1, "POST")) %>%
  dplyr::mutate(core_status_together = replace(core_status_together, core_status == 1 & core_status_mean == 1, "BOTH")) %>%
  dplyr::mutate(core_status_together = replace(core_status_together, core_status == 0 & core_status_mean == 0, "NONCORE"))
 
# Plot PCA 
res.pca <- prcomp(Z, center = FALSE)
ind.coord <- res.pca$x
# head(ind.coord[, 1:2])
# aa = summary(res.pca)
# First 10 PC's (Var, Cumulative)
# aa$importance[2, 1:10]
# aa$importance[2, 1:10]

ind.coord_2 = data.frame(ind.coord[, 1:2], geno_ref_compare_2$core_status_together)

p = ggplot(ind.coord_2, aes(x = ind.coord_2[, 1], y = ind.coord_2[, 2], group = ind.coord_2[, 3])) + geom_point(aes(shape = ind.coord_2[, 3], color = ind.coord_2[, 3], alpha = ind.coord_2[, 3])) + 
  theme_classic() + labs(y = "PC 2 (4.3%)", x = "PC 1 (5.1%)", colour = "Knots", shape = "Knots", alpha = "Knots") +
  scale_color_manual(values = c("darkgreen", "lightgreen", "darkblue", "brown" )) +
  scale_alpha_manual(values = c(1, 0.2, 1, 1))


#### D1) Naive Method 1: Highest Diagonal Value ####

# n_final = read.table(file = "coresize.txt", header = TRUE)
# j = 8
# varex = colnames(n_final)[j]
# ncore = as.numeric(n_final[j])

rm(final_data4)
final_data4 = final_data3
rm(geno_refNaive1)
geno_refNaive1 = geno_ref[,-3]

# Get diagonal of G (either form G or Z) 
geno_ref_diag = geno_ref
geno_ref_diag$v_prior <- rowSums(Z^2)
# When working with Z, we actually get rowSums(Z^2)[i] = G[i,i]*k
# geno_ref_diag$v_prior <- diag(G)

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
  
  # Export the APY inverse created in R, to use as external file in blupf90
  # use option "user_file" in blupf90
  # This could be avoided, but APY option is not available in the public version
  
  # Clean intermediate files
  rm(Ginv_full2, G_out1, G_out2)
  
  # Create datafile for BLUPF90
  pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
    arrange(apy_order)
  
  # Export data:
  write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
  # It is simulated data, so we don't run renumf90
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
# ~ 1h on Mac

rm(apyginvlist)

save.image("after_naive1.RData")
# load("after_naive1.RData")

# DataSet1 <- new.env() 
# load("after_prior.RData", envir = DataSet1)
# ls(DataSet1)
# head(DataSet1$final_data3)
# tmp2APY = DataSet1$final_data3 %>%
#  dplyr::select(starts_with("APY_e99"), starts_with("SE_APY_e99"), starts_with("Prior_e99"), starts_with("SE_Prior_e99"), "Aid")
# final_data4 = inner_join(final_data4, tmp2APY, by = "Aid")
# geno_refPrior = DataSet1$geno_refPrior 
# rm(DataSet1)


  
  
#### D2) Naive Method 2: Random Sample with Weights (Diagonal Values) ####
  
rm(final_data5)
final_data5 = final_data4
rm(geno_refNaive2)
geno_refNaive2 = geno_ref[,-3]

# n_final = read.table(file = "coresize.txt", header = TRUE)
# j = 8
# varex = colnames(n_final)[j]
# ncore = as.numeric(n_final[j])

# Same as in the D1) Naive Method 1 
# Get diagonal of G (either form G or Z) 
# geno_ref_diag = geno_ref
# geno_ref_diag$v_prior <- rowSums(Z^2)
# When working with Z, we actually get rowSums(Z^2)[i] = G[i,i]*k
# geno_ref_diag$v_prior <- diag(G)

start_time = Sys.time()
for(j in 1:(length(n_final))){
  
  varex = colnames(n_final)[j]
  ncore = n_final[j]
  
  # Run 5 replicates of random core sampling, for each number of core animals 
  for(apy_rep in 1:5){
    # Select random animals as core 
    # set.seed(12345)
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
    
    # Export the APY inverse created in R, to use as external file in blupf90
    # use option "user_file" in blupf90
    # This could be avoided, but APY option is not available in the public version
    
    # Clean intermediate files
    rm(Ginv_full2, G_out1, G_out2)
    
    # Create datafile for BLUPF90
    pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
      arrange(apy_order)
    
    # Export data:
    write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
    
    # It is simulated data, so we don't run renumf90
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
# ~ 5.5h on Mac

rm(apyginvlist)

save.image("after_naive2.RData")
# load("after_naive2.RData")

# For the last version load this as SVD was not re-tested: 4.3.3022.

# rm(tmp2APY)
# tmp2APY = final_data5 %>%
#  dplyr::select(starts_with("Weight_"), starts_with("SE_Weight_"), "Aid")
# final_data5_backup = final_data5
# rm(final_data5)
# final_data5 = inner_join(final_data4, tmp2APY, by = "Aid")

#### E) Optimized core selection - Prior with SVD of Z ####

rm(final_data6)
final_data6 = final_data5
rm(geno_refPriorSVD)
geno_refPriorSVD = geno_ref[,-3]

# Here is another moment where we chose how many dimensions or PCs we select
# Example 99% variance = 3130
# Note; before A = Z

A <- svd_Z$u[,1:3130] * rep(svd_Z$d[1:3130], each = nrow(svd_Z$u))
dim(A)
# Equivalent. What is the proper term: Reduced Rank PCA?
# aa = svd_Z$u[,1] * rep(svd_Z$d[1], each = nrow(svd_Z$u))
# aa = res.pca$x[,1]

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

# Export the APY inverse created in R, to use as external file in blupf90
# use option "user_file" in blupf90
# This could be avoided, but APY option is not available in the public version

# Clean intermediate files
rm(Ginv_full2, G_out1, G_out2)

# Create datafile for BLUPF90
pheno_export = inner_join(phenotypes, apyginvlist$index_file[, c(1,4)], by = c("Aid" = "geno_id")) %>%
  arrange(apy_order)

# Export data:
write.table(cbind(pheno_export[,c("Aid", "Trait", "apy_order")], rep(1, nrow(pheno_export))), "Blupf90.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# It is simulated data, so we don't run renumf90
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
# ~2h on Mac

rm(apyginvlist, A)

save.image("after_priorSVD.RData")
# load("after_priorSVD.RData")


#### F) Updating with arrival of new data ####

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


  