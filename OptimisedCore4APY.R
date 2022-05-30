# "Optimal core definitions for the APY model"
# Pocrnic, Lindgren, Tolhust, Herring, Gorjanc (2022)
# DOI: 


# CLEAR FOR OPEN VERSION ON GITHUB

rm(list = ls())

library("tidyverse")
library("AlphaSimR")
library("Matrix")
library("umap")
library("gridExtra")
library("VennDiagram")


library("scoringRules")
library("verification")
# library("profvis")

library("pryr")
# mem_used()


setwd("~/Desktop/test_optiapy/")
getwd()

# source("functions_opticore.R")

# Note: blupf90 has to be available

# Create blupf90 parameter files and bash scripts
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

# Way (1) - Empirical: Sequential Update vs. Refreshed
  # See if those new 5 are different than if i would select 15 combined
  # At what point is worth updating the new core
# Way (2) - Find some metric
# Somehow find a-priory measure if updating is needed 
# e.g. if the variance of new Z is X-times higher than max old one; should be updated withj 

# Consideration
# Number of new core
# Are we just adding or should we remove from old core

  
#### Analyses #### 

#### Calculate estimated PA and estimated MS ####
final_data6_MS = final_data6[ ,1:82] %>%
  dplyr::select(-starts_with("SE")) %>%
  dplyr::select(-ends_with("_Rep2") & -ends_with("_Rep3") & -ends_with("_Rep4") & -ends_with("_Rep5"))

# Calculate MS for Direct Inverse
# Make a function here, don't be lazy!!! 

final_data6_MS$Full_sire = 0
final_data6_MS$Full_sire_e10 = 0
final_data6_MS$Full_sire_e30 = 0
final_data6_MS$Full_sire_e50 = 0
final_data6_MS$Full_sire_e70 = 0
final_data6_MS$Full_sire_e90 = 0
final_data6_MS$Full_sire_e95 = 0
final_data6_MS$Full_sire_e98 = 0

final_data6_MS$Full_dam = 0
final_data6_MS$Full_dam_e10 = 0
final_data6_MS$Full_dam_e30 = 0
final_data6_MS$Full_dam_e50 = 0
final_data6_MS$Full_dam_e70 = 0
final_data6_MS$Full_dam_e90 = 0
final_data6_MS$Full_dam_e95 = 0
final_data6_MS$Full_dam_e98 = 0

final_data6_MS$Full_EMS = 0
final_data6_MS$Full_EMS_e10 = 0
final_data6_MS$Full_EMS_e30 = 0
final_data6_MS$Full_EMS_e50 = 0
final_data6_MS$Full_EMS_e70 = 0
final_data6_MS$Full_EMS_e90 = 0
final_data6_MS$Full_EMS_e95 = 0
final_data6_MS$Full_EMS_e98 = 0

for (i in 1:nrow(final_data6_MS)){
  if(final_data6_MS$Sid[i] > final_data6_MS$Aid[1] )
  { 
    final_data6_MS$Full_sire[i] = final_data6_MS$Full[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e10[i] = final_data6_MS$APY_e10_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e30[i] = final_data6_MS$APY_e30_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e50[i] = final_data6_MS$APY_e50_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e70[i] = final_data6_MS$APY_e70_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e90[i] = final_data6_MS$APY_e90_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e95[i] = final_data6_MS$APY_e95_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
    final_data6_MS$Full_sire_e98[i] = final_data6_MS$APY_e98_Rep1[final_data6_MS$Aid==final_data6_MS$Sid[i] ]
  } 
  if( final_data6_MS$Did[i] > final_data6_MS$Aid[1])
  {
    final_data6_MS$Full_dam[i] = final_data6_MS$Full[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e10[i] = final_data6_MS$APY_e10_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e30[i] = final_data6_MS$APY_e30_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e50[i] = final_data6_MS$APY_e50_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e70[i] = final_data6_MS$APY_e70_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e90[i] = final_data6_MS$APY_e90_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e95[i] = final_data6_MS$APY_e95_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
    final_data6_MS$Full_dam_e98[i] = final_data6_MS$APY_e98_Rep1[final_data6_MS$Aid==final_data6_MS$Did[i] ]
  }
  
  final_data6_MS$Full_EMS[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire[i] + final_data6_MS$Full_dam[i])
  final_data6_MS$Full_EMS_e10[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e10[i] + final_data6_MS$Full_dam_e10[i]) 
  final_data6_MS$Full_EMS_e30[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e30[i] + final_data6_MS$Full_dam_e30[i]) 
  final_data6_MS$Full_EMS_e50[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e50[i] + final_data6_MS$Full_dam_e50[i]) 
  final_data6_MS$Full_EMS_e70[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e70[i] + final_data6_MS$Full_dam_e70[i]) 
  final_data6_MS$Full_EMS_e90[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e90[i] + final_data6_MS$Full_dam_e90[i]) 
  final_data6_MS$Full_EMS_e95[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e95[i] + final_data6_MS$Full_dam_e95[i]) 
  final_data6_MS$Full_EMS_e98[i] = final_data6_MS$Full[i] - 0.5 * (final_data6_MS$Full_sire_e98[i] + final_data6_MS$Full_dam_e98[i]) 
}



final_data6_MS_G20 = final_data6_MS %>%
  filter(Generation == 20)

rez_MS = as_tibble(colnames(n_final))
rez_MS$nvar = as.numeric(str_sub(rez_MS$value, 2))  
rez_MS$ncore = n_final[1,]
rez_MS$Full = rep(round(cor(final_data6_MS_G20$TBV, final_data6_MS_G20$Full, use = "pairwise.complete.obs"), digits = 2), each = length(n_final))
rez_MS$MS_Full = rep(round(cor(final_data6_MS_G20$TMS, final_data6_MS_G20$Full_EMS, use = "pairwise.complete.obs"), digits = 2), each = length(n_final))


# Calculate ACC for APY runs
cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98")
rez_apy = data.frame(matrix(ncol = 1, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = final_data6_MS_G20 %>%
    dplyr::select(starts_with(cnm[x]), "TBV")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = round(cor(tmpX$TBV, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  rescor_mean = mean(rescor)
  
  rez_apy = rbind(rez_apy, c(rescor_mean))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("APY")

rez_MS = bind_cols(rez_MS, rez_apy)


# Calculate MS for APY runs
cnm = c("Full_EMS_e10", "Full_EMS_e30", "Full_EMS_e50", "Full_EMS_e70", "Full_EMS_e90", "Full_EMS_e95", "Full_EMS_e98")
rez_apy = data.frame(matrix(ncol = 1, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = final_data6_MS_G20 %>%
    dplyr::select(starts_with(cnm[x]), "TMS")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = round(cor(tmpX$TMS, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  rescor_mean = mean(rescor)
  
  rez_apy = rbind(rez_apy, c(rescor_mean))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("APY_EMS")

rez_MS = bind_cols(rez_MS, rez_apy)

colnames(rez_MS) = c("nDim", "nVar", "nCore", "AccFull", "MS_AccFull", "AccAPY", "MS_AccAPY")  

#### MS Accuracy Graphs ####
p = ggplot(rez_MS, aes(x = nCore, y = AccAPY)) + theme_classic() + geom_line(aes(colour = "APY")) + 
  geom_line(aes(x = nCore, y = AccFull, colour = "Full")) + 
  geom_line(aes(x = nCore, y = MS_AccFull, colour = "MS_Full")) + 
  geom_line(aes(x = nCore, y = MS_AccAPY, colour = "MS_APY")) +
  xlab("Number of core animals") + ylab("Correlation with TBV or TMS") + labs(colour = "Accuracy") +
  theme(legend.position = "top") + scale_color_manual(values = color_roslin) + ylim(0.1, 0.8)

ggsave(plot = p + PaperTheme, filename = "Accuracy_MST.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm")






#### Accuracy ####

final_data6_backup = final_data6
# final_data6 = final_data6_backup

# For the last version: 4.3.3022.
final_data6 = final_data5

# Distribution of EBV vs. TBV
ggplot(final_data6, aes(x = TBV)) + 
  geom_histogram(binwidth = 1, fill = "red", alpha = 0.2) +     
  geom_histogram(aes(x = Full), binwidth = 1, fill = "blue", alpha = 0.2)

mean(final_data6$TBV)
sd(final_data6$TBV)
mean(final_data6$Full)
sd(final_data6$Full)

# Standardize TBV 
TBV_Mean = mean(final_data6$TBV, na.rm = TRUE) 
TBV_SD = sd(final_data6$TBV, na.rm = TRUE) 
final_data6$TBVSTD = (final_data6$TBV - TBV_Mean) / TBV_SD

ggplot(final_data6, aes(x = Full, y = TBVSTD)) + geom_point(size = 1, alpha = 0.2) + 
  geom_smooth(method="lm", formula = y ~ x, col = "red", size = 1)

# While EBV will always be “centered”, with expectation “0", the TBV from ASR can have gv_mu much different than zero - especially as we start calculating EBV after N past generations of selection we do not “record”

# Standardize all GEBV to units of SD(TBV)
for(k in 11:(length(final_data6)-1)){
  final_data6[k] = final_data6[k] / TBV_SD
}


##### Accuracy for validation animals only (Generation 20) #####
tmp_ac = final_data6 %>%
  filter(Generation == 20)

rezultati = NULL

rezultati = as_tibble(colnames(n_final))
rezultati$nvar = as.numeric(str_sub(rezultati$value, 2))  
rezultati$ncore = n_final[1,]
# Calculate accuracy for Direct inverse
rezultati$Full = rep(cor(tmp_ac$TBV, tmp_ac$Full, use = "pairwise.complete.obs"), each = length(n_final))

# Note: Create a function for calculating this accuracies - current code is messy!!!
# For correlations(GEBV_FULL, GEBV_APY) change "TBV" to "Full" in the code 

# Calculate accuracy for random scenarios
cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98", "APY_e99")
rez_apy = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = tmp_ac %>%
    dplyr::select(starts_with(cnm[x]), "TBV")

  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = cor(tmpX$TBV, tmpX[k], use = "pairwise.complete.obs")
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }

  rescor_mean = mean(rescor)
  rescor_L = quantile(rescor, na.rm=TRUE, probs=0.025)
  rescor_U = quantile(rescor, na.rm=TRUE, probs=0.975)
  
  rez_apy = rbind(rez_apy, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("APY", "APY_L", "APY_U")

rezultati = bind_cols(rezultati, rez_apy)


# round(cor(as.matrix(cbind(tmp_ac$APY_e98_Rep1, tmp_ac$APY_e98_Rep2, tmp_ac$APY_e98_Rep3, tmp_ac$APY_e98_Rep4, tmp_ac$APY_e98_Rep5)), use = "pairwise.complete.obs"), digits = 2)


# Calculate accuracy for optimal scenarios
tmpX = tmp_ac %>%
  dplyr::select(starts_with("Prior"), "TBV")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  rescor_new = cor(tmpX$TBV, tmpX[k], use = "pairwise.complete.obs")
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Prior"
rezultati = bind_cols(rezultati, rescor)

# Calculate accuracy for diagonal scenarios
tmpX = tmp_ac %>%
  dplyr::select(starts_with("Diag"), "TBV")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  rescor_new = cor(tmpX$TBV, tmpX[k], use = "pairwise.complete.obs")
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Diag"
rezultati = bind_cols(rezultati, rescor)

# Calculate accuracy for weighted scenarios
cnm = c("Weight_e10", "Weight_e30", "Weight_e50", "Weight_e70", "Weight_e90", "Weight_e95", "Weight_e98", "Weight_e99")
rez_wgt = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = tmp_ac %>%
    dplyr::select(starts_with(cnm[x]), "TBV")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = cor(tmpX$TBV, tmpX[k], use = "pairwise.complete.obs")
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  rescor_mean = mean(rescor)
  rescor_L = quantile(rescor, na.rm=TRUE, probs=0.025)
  rescor_U = quantile(rescor, na.rm=TRUE, probs=0.975)
  rez_wgt = rbind(rez_wgt, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_wgt)[x] = cnm[x]
  
}
colnames(rez_wgt) = c("WGT", "WGT_L", "WGT_U")

rezultati = bind_cols(rezultati, rez_wgt)

# round(cor(as.matrix(cbind(tmp_ac$Weight_e98_Rep1, tmp_ac$Weight_e98_Rep2, tmp_ac$Weight_e98_Rep3, tmp_ac$Weight_e98_Rep4, tmp_ac$Weight_e98_Rep5)), use = "pairwise.complete.obs"), digits = 2)

# Make a backup
bu_rezultati_prediction = rezultati 

# Calculate accuracy for SVD Prior runs
# Skip putting this to graph for the paper, only table needed (14.03.2022.)
tmpX = tmp_ac %>%
  dplyr::select(starts_with("SVD"), "TBV")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  rescor_new = round(cor(tmpX$TBV, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "SVDPrior"
rezultati = bind_cols(rezultati, rescor)


##### Accuracy Graphs #####
# Accuracy for paper (all scenarios)
p = ggplot(rezultati, aes(x = nvar, y = APY)) + theme_classic(base_size = 12, base_family = "Times New Roman") + 
  geom_line(aes(colour = "Random", linetype = "Random")) + 
  geom_ribbon(aes(x = nvar, ymin = APY_L, ymax = APY_U), alpha = .1, colour = NA) +
  geom_line(aes(x = nvar, y = Prior, colour = "Optimised", linetype = "Optimised")) + 
  geom_line(aes(x = nvar, y = Full, colour = "Full", linetype = "Full")) +
  geom_line(aes(x = nvar, y = Diag, colour = "Diagonal", linetype = "Diagonal")) + 
  geom_line(aes(x = nvar, y = WGT, colour = "Weighted", linetype = "Weighted")) + 
  geom_ribbon(aes(x = nvar, ymin = WGT_L, ymax = WGT_U), alpha = .1, colour = NA) +
  xlab("% Variance explained") + ylab("Cor(GEBV, TBV)") + labs(colour = "Core selection", linetype = "Core selection") +
  theme(legend.position = "top") + scale_color_manual(breaks = c("Random", "Diagonal", "Weighted", "Optimised", "Full"), values = c(color_roslin[1], color_roslin[2], color_roslin[3], color_roslin[5], color_roslin[4] )) +
  scale_linetype_manual(breaks = c("Random", "Diagonal", "Weighted", "Optimised", "Full"), values = c("dashed", "dotdash", "dotted", "solid", "solid") ) +
  ylim(0.0, 0.80) + xlim(10, 100)


ggsave(file = "accuracy_simulation.png", p, height = PaperSize, width = PaperSize * 1.5, unit = "cm")


#### CRPS ####
# Gregor: 
## We get posterior distribution for each breeding value. For each breeding value we also have the true value.
## We need to calculate CRPS for each breeding value - it tells us how close and concentrated posteriors distribution is around the true value.

rez_crps = as_tibble(colnames(n_final))
rez_crps$nvar = as.numeric(str_sub(rez_crps$value, 2))  
rez_crps$ncore = n_final[1,]
rez_crps$Full = rep(round(mean( scoringRules::crps(as.vector(tmp_ac$TBVSTD), family = "normal",  location = tmp_ac$Full, scale = tmp_ac$SE_Full) )
, digits = 2), each = length(n_final))


# Calculate CRPS for APY runs
cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98")
cnmse = c("SE_APY_e10", "SE_APY_e30", "SE_APY_e50", "SE_APY_e70", "SE_APY_e90", "SE_APY_e95", "SE_APY_e98")
rez_apy_crps = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  # x = 1
  tmpX = tmp_ac %>%
    dplyr::select(starts_with(cnm[x]), starts_with(cnmse[x]), "TBVSTD")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:5){
    # k = 1
    rescor_new = round(mean( scoringRules::crps(as.vector(tmpX$TBVSTD), family = "normal",  location = dplyr::pull(tmpX[k]), scale = dplyr::pull(tmpX[k+5])) ), digits = 2)
    rescor_new = as.data.frame(rescor_new)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  rescor_mean = round(mean(rescor$rescor_new), 2)
  rescor_L = round(quantile(rescor$rescor_new, na.rm=TRUE, probs=0.025), 2)
  rescor_U = round(quantile(rescor$rescor_new, na.rm=TRUE, probs=0.975), 2)
  
  rez_apy_crps = rbind(rez_apy_crps, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_apy_crps)[x] = cnm[x]
  
}
colnames(rez_apy_crps) = c("APY", "APY_L", "APY_U")

rez_crps = bind_cols(rez_crps, rez_apy_crps)

# Calculate CRPS for Prior runs
tmpX = tmp_ac %>%
  dplyr::select(starts_with("Prior"), starts_with("SE_Prior"), "TBVSTD")

rescor = NULL
rescor_new = NULL
for(k in 1:7){
  rescor_new = round(mean( scoringRules::crps(as.vector(tmpX$TBVSTD), family = "normal",  location = dplyr::pull(tmpX[k]), scale = dplyr::pull(tmpX[k+7])) ), digits = 2)
  rescor_new = as.data.frame(rescor_new)
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}

rescor = as.data.frame(rescor)
colnames(rescor) = "Prior"
rez_crps = bind_cols(rez_crps, rescor)


# Calculate CRPS for Diag runs
tmpX = tmp_ac %>%
  dplyr::select(starts_with("Diag"), starts_with("SE_Diag"), "TBVSTD")

rescor = NULL
rescor_new = NULL
for(k in 1:7){
  rescor_new = round(mean( scoringRules::crps(as.vector(tmpX$TBVSTD), family = "normal",  location = dplyr::pull(tmpX[k]), scale = dplyr::pull(tmpX[k+7])) ), digits = 2)
  rescor_new = as.data.frame(rescor_new)
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Diag"
rez_crps = bind_cols(rez_crps, rescor)

# Calculate CRPS for Weighted Diag runs
cnm = c("Weight_e10", "Weight_e30", "Weight_e50", "Weight_e70", "Weight_e90", "Weight_e95", "Weight_e98")
cnmse = c("SE_Weight_e10", "SE_Weight_e30", "SE_Weight_e50", "SE_Weight_e70", "SE_Weight_e90", "SE_Weight_e95", "SE_Weight_e98")
rez_wgt_crps = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = tmp_ac %>%
    dplyr::select(starts_with(cnm[x]), starts_with(cnmse[x]), "TBVSTD")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:5){
    rescor_new = round(mean( scoringRules::crps(as.vector(tmpX$TBVSTD), family = "normal",  location = dplyr::pull(tmpX[k]), scale = dplyr::pull(tmpX[k+5])) ), digits = 2)
    rescor_new = as.data.frame(rescor_new)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  rescor_mean = round(mean(rescor$rescor_new), 2)
  rescor_L = round(quantile(rescor$rescor_new, na.rm=TRUE, probs=0.025), 2)
  rescor_U = round(quantile(rescor$rescor_new, na.rm=TRUE, probs=0.975), 2)
  
  rez_wgt_crps = rbind(rez_wgt_crps, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_wgt_crps)[x] = cnm[x]
  
}
colnames(rez_wgt_crps) = c("WGT", "WGT_L", "WGT_U")

rez_crps = bind_cols(rez_crps, rez_wgt_crps)



# Calculate CRPS for SVD Prior runs
tmpX = tmp_ac %>%
  dplyr::select(starts_with("SVDPrior"), starts_with("SE_SVDPrior"), "TBVSTD")

rescor = NULL
rescor_new = NULL
for(k in 1:7){
  rescor_new = round(mean( scoringRules::crps(as.vector(tmpX$TBVSTD), family = "normal",  location = dplyr::pull(tmpX[k]), scale = dplyr::pull(tmpX[k+7])) ), digits = 2)
  rescor_new = as.data.frame(rescor_new)
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}

rescor = as.data.frame(rescor)
colnames(rescor) = "SVDPrior"
rez_crps = bind_cols(rez_crps, rescor)




p = ggplot(rez_crps, aes(x = ncore, y = APY)) + theme_classic() + geom_line(aes(colour = "Random")) + 
  geom_ribbon(aes(x = ncore, ymin = APY_L, ymax = APY_U), alpha = .2, colour = NA) +
  geom_line(aes(x = ncore, y = Full, colour = "Full")) + 
  geom_line(aes(x = ncore, y = Prior, colour = "Optimal")) + 
  xlab("Number of core animals") + ylab("CRPS") + labs(colour = "GBLUP scenario") +
  theme(legend.position = "top") + scale_color_manual(values = color_roslin) +
  ylim(0.2, 0.8) + xlim(0.0, 2500)


ggsave(plot = p + PaperTheme, filename = "crps.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "crps_pres.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm")

#### Correlation Graphs - All animals ####  
rezultati = NULL

solutions_all = final_data6 

# Calculate corr for optimal scenarios
tmpX = solutions_all %>%
  dplyr::select(starts_with("Prior"), "Full")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
  rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Prior"
rezultati = bind_cols(rezultati, rescor)


# Calculate corr for random scenarios
cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98", "APY_e99")
rez_apy = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = solutions_all %>%
    dplyr::select(starts_with(cnm[x]), "Full")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
    rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  # rescor_mean = round(mean(rescor), 2)
  # rescor_L = round(quantile(rescor, na.rm=TRUE, probs=0.025), 2)
  # rescor_U = round(quantile(rescor, na.rm=TRUE, probs=0.975), 2)
  
  rescor_mean = mean(rescor)
  rescor_L = quantile(rescor, na.rm=TRUE, probs=0.025)
  rescor_U = quantile(rescor, na.rm=TRUE, probs=0.975)
  
  rez_apy = rbind(rez_apy, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("APY", "APY_L", "APY_U")

rezultati = bind_cols(rezultati, rez_apy)


# Calculate corr for diagonal scenarios
tmpX = solutions_all %>%
  dplyr::select(starts_with("Diag"), "Full")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
  rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Diagonal"

rezultati = bind_cols(rezultati, rescor)


# Calculate corr for weighted scenarios
cnm = c("Weight_e10", "Weight_e30", "Weight_e50", "Weight_e70", "Weight_e90", "Weight_e95", "Weight_e98", "Weight_e99")
rez_apy = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = solutions_all %>%
    dplyr::select(starts_with(cnm[x]), "Full")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
    rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  # rescor_mean = round(mean(rescor), 2)
  # rescor_L = round(quantile(rescor, na.rm=TRUE, probs=0.025), 2)
  # rescor_U = round(quantile(rescor, na.rm=TRUE, probs=0.975), 2)
  
  rescor_mean = mean(rescor)
  rescor_L = quantile(rescor, na.rm=TRUE, probs=0.025)
  rescor_U = quantile(rescor, na.rm=TRUE, probs=0.975)
  
  rez_apy = rbind(rez_apy, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("Weight", "Weight_L", "Weight_U")

rezultati = bind_cols(rezultati, rez_apy)

#### Finalize the corr table and plot - All animals #### 
rezultati$ncore = t(n_final[1,])
rezultati$nvar = as.numeric(str_sub(rownames(rezultati), 8)) 

# Make a backup
bu_rezultati_all = rezultati 


# Graph for paper (all scenarios)
p1_all = ggplot(rezultati, aes(x = nvar, y = APY)) + theme_classic(base_size = 12, base_family = "Times New Roman") + 
  geom_line(aes(colour = "Random", linetype = "Random")) + 
  geom_ribbon(aes(x = nvar, ymin = APY_L, ymax = APY_U), alpha = .2, colour = NA) +
  geom_line(aes(x = nvar, y = Prior, colour = "Optimised", linetype = "Optimised")) + 
  geom_line(aes(x = nvar, y = Diagonal, colour = "Diagonal", linetype = "Diagonal")) + 
  geom_line(aes(x = nvar, y = Weight, colour = "Weighted", linetype = "Weighted")) + 
  geom_ribbon(aes(x = nvar, ymin = Weight_L, ymax = Weight_U), alpha = .2, colour = NA) +
  xlab("% Variance explained") + ylab("Cor(Full, APY)") + labs(colour = "Core selection", linetype = "Core selection") +
  theme(legend.position = "top") + 
  scale_color_manual(breaks = c("Random", "Diagonal", "Weighted", "Optimised"), values = c(color_roslin[1], color_roslin[2], color_roslin[3], color_roslin[5] )) +
  scale_linetype_manual(breaks = c("Random", "Diagonal", "Weighted", "Optimised"), values = c("dashed", "dotdash", "dotted", "solid") ) +
  ylim(0.0, 1.0) + xlim(10, 100)


#### Correlation Graphs - Validation animals ####

accuracy_optimal = final_data6 %>%
  filter(Generation == 20)

rezultati = NULL

# Calculate corr for optimal scenarios
tmpX = accuracy_optimal %>%
  dplyr::select(starts_with("Prior"), "Full")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
  rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Prior"
rezultati = bind_cols(rezultati, rescor)


# Calculate corr for random scenarios
cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98", "APY_e99")
rez_apy = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = accuracy_optimal %>%
    dplyr::select(starts_with(cnm[x]), "Full")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
    rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  # rescor_mean = round(mean(rescor), 2)
  # rescor_L = round(quantile(rescor, na.rm=TRUE, probs=0.025), 2)
  # rescor_U = round(quantile(rescor, na.rm=TRUE, probs=0.975), 2)
  
  rescor_mean = mean(rescor)
  rescor_L = quantile(rescor, na.rm=TRUE, probs=0.025)
  rescor_U = quantile(rescor, na.rm=TRUE, probs=0.975)
  
  rez_apy = rbind(rez_apy, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("APY", "APY_L", "APY_U")

rezultati = bind_cols(rezultati, rez_apy)


# Calculate corr for diagonal scenarios
tmpX = accuracy_optimal %>%
  dplyr::select(starts_with("Diag"), "Full")

rescor = NULL
rescor_new = NULL
for(k in 1:(length(tmpX)-1)){
  # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
  rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
  rownames(rescor_new) = colnames(tmpX)[k]
  rescor = rbind(rescor, rescor_new)
}
rescor = as.data.frame(rescor)
colnames(rescor) = "Diagonal"

rezultati = bind_cols(rezultati, rescor)


# Calculate corr for weighted scenarios
cnm = c("Weight_e10", "Weight_e30", "Weight_e50", "Weight_e70", "Weight_e90", "Weight_e95", "Weight_e98", "Weight_e99")
rez_apy = data.frame(matrix(ncol = 3, nrow = 0))
for(x in 1:length(cnm)){
  
  tmpX = accuracy_optimal %>%
    dplyr::select(starts_with(cnm[x]), "Full")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    # rescor_new = round(cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs"), digits = 2)
    rescor_new = cor(tmpX$Full, tmpX[k], use = "pairwise.complete.obs")
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  
  # rescor_mean = round(mean(rescor), 2)
  # rescor_L = round(quantile(rescor, na.rm=TRUE, probs=0.025), 2)
  # rescor_U = round(quantile(rescor, na.rm=TRUE, probs=0.975), 2)
  
  rescor_mean = mean(rescor)
  rescor_L = quantile(rescor, na.rm=TRUE, probs=0.025)
  rescor_U = quantile(rescor, na.rm=TRUE, probs=0.975)
  
  rez_apy = rbind(rez_apy, c(rescor_mean, rescor_L, rescor_U))
  rownames(rez_apy)[x] = cnm[x]
  
}
colnames(rez_apy) = c("Weight", "Weight_L", "Weight_U")

rezultati = bind_cols(rezultati, rez_apy)

#### Finalize the corr table and plot - Validation animals #### 
rezultati$ncore = t(n_final[1,])
rezultati$nvar = as.numeric(str_sub(rownames(rezultati), 8)) 

# Make a backup
bu_rezultati_validation = rezultati 

# Graph for paper (all scenarios)
p2_all = ggplot(rezultati, aes(x = nvar, y = APY)) + theme_classic(base_size = 12, base_family = "Times New Roman") + 
  geom_line(aes(colour = "Random", linetype = "Random")) + 
  geom_ribbon(aes(x = nvar, ymin = APY_L, ymax = APY_U), alpha = .2, colour = NA) +
  geom_line(aes(x = nvar, y = Prior, colour = "Optimised", linetype = "Optimised")) + 
  geom_line(aes(x = nvar, y = Diagonal, colour = "Diagonal", linetype = "Diagonal")) + 
  geom_line(aes(x = nvar, y = Weight, colour = "Weighted", linetype = "Weighted")) + 
  geom_ribbon(aes(x = nvar, ymin = Weight_L, ymax = Weight_U), alpha = .2, colour = NA) +
  xlab("% Variance explained") + ylab("Cor(Full, APY)") + labs(colour = "Core selection", linetype = "Core selection") +
  theme(legend.position = "top") + 
  scale_color_manual(breaks = c("Random", "Diagonal", "Weighted", "Optimised"), values = c(color_roslin[1], color_roslin[2], color_roslin[3], color_roslin[5] )) +
  scale_linetype_manual(breaks = c("Random", "Diagonal", "Weighted", "Optimised"), values = c("dashed", "dotdash", "dotted", "solid") ) +
  ylim(0.0, 1.0) + xlim(10, 100)


####  Plot Correlations together ####

get_legend = function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# For paper (all scenarios)
p1a = p1_all + labs(tag = "a)") 
p2a = p2_all + labs(y = "", x = "", tag = "b)") 
legend = get_legend(p2a)

panel_1 = grid.arrange(p1a + theme(legend.position="none"), 
                       p2a + theme(legend.position="none"),
                       legend,
                       nrow = 2, ncol = 2,
                       layout_matrix = rbind(c(3,3), c(1,2)),
                       heights = c(1, 10),
                       widths = c(10, 10) )

ggsave(file = "Panel_1_sim.png", panel_1, height = PaperSize, width = PaperSize * 1.5, unit = "cm")


#### Testing CRPS ####

# All of these give the same results: 

cr1 = verification::crps(obs = final_data6$TBVSTD, pred = cbind(final_data6$Full, final_data6$SE_Full) )$CRPS
cr2 = verification::crps(obs = final_data6$TBVSTD, pred = cbind(final_data6$Full, final_data6$SE_Full) )

cr3 = mean( scoringRules::crps(as.vector(final_data6$TBVSTD), family = "normal",  location = final_data6$Full, scale = final_data6$SE_Full) )
cr4 = mean( scoringRules::crps_norm(final_data6$TBVSTD, mean = final_data6$Full, sd = final_data6$SE_Full) )

mean(cr2$crps) == cr1 
mean(cr2$crps) == cr3 
mean(cr2$crps) == cr4

# mean( scoringRules::crps(as.vector(final_data$TBV), family = "normal",  mean = mean(final_data$Full), sd = mean(final_data$SE_Full)) )
# mean( scoringRules::crps_norm(final_data$TBV, mean = mean(final_data$Full), sd = mean(final_data$SE_Full)) )


#### Panel plot: UMAP with core selections ####

# ind.coord_3 = data.frame(ind.coord[, 1:2], geno_refPrior$Prior_e10, geno_refPrior$geno_id)
# ind.coord_3 = ind.coord_3[order(ind.coord_3[, 3]),]

ind.coord_p1 = left_join(ind.coord_333, geno_refPrior, by = c("row.names.res.umap.layout." = "geno_id"))
ind.coord_p1 = ind.coord_p1[order(ind.coord_p1$Prior_e10),]


p1 = ggplot(ind.coord_p1, aes(x = ind.coord_p1[, 1], y = ind.coord_p1[, 2], color = factor(Prior_e10), alpha = factor(Prior_e10) ) ) + 
  geom_point(aes(shape = factor(Prior_e10), color = factor(Prior_e10), alpha = factor(Prior_e10) ) ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + scale_alpha_manual(limits = c("0", "1"), values = c(0.2, 1)) + scale_color_manual(values = c("aquamarine", "darkblue")) +
  labs(y = "", x = "", colour = "Core Status", alpha = "Core Status", shape = "Core Status", tag = "d)") +
  theme(legend.position = "top")


ind.coord_p2 = left_join(ind.coord_333, geno_refAPY, by = c("row.names.res.umap.layout." = "geno_id"))
ind.coord_p2 = ind.coord_p2[order(ind.coord_p2$APY_e10_Rep1, ind.coord_p2$APY_e10_Rep2, ind.coord_p2$APY_e10_Rep3, ind.coord_p2$APY_e10_Rep4, ind.coord_p2$APY_e10_Rep5),]

ind.coord_p2 = ind.coord_p2 %>%
  dplyr::mutate(APY_e10_Rep2 = ifelse(APY_e10_Rep2 == 1, 2, NA))  %>%
  dplyr::mutate(APY_e10_Rep3 = ifelse(APY_e10_Rep3 == 1, 3, NA))  %>%
  dplyr::mutate(APY_e10_Rep4 = ifelse(APY_e10_Rep4 == 1, 4, NA))  %>%
  dplyr::mutate(APY_e10_Rep5 = ifelse(APY_e10_Rep5 == 1, 5, NA))
  
p2 = ggplot(ind.coord_p2, aes(x = ind.coord_p2[, 1], y = ind.coord_p2[, 2] ) ) + 
  geom_point(aes(shape = factor(APY_e10_Rep1), color = factor(APY_e10_Rep1), alpha = factor(APY_e10_Rep1) ) ) +
  geom_point(aes(shape = factor(APY_e10_Rep2), color = factor(APY_e10_Rep2), alpha = factor(APY_e10_Rep2) ) ) +
  geom_point(aes(shape = factor(APY_e10_Rep3), color = factor(APY_e10_Rep3), alpha = factor(APY_e10_Rep3) ) ) +
  geom_point(aes(shape = factor(APY_e10_Rep4), color = factor(APY_e10_Rep4), alpha = factor(APY_e10_Rep4) ) ) +
  geom_point(aes(shape = factor(APY_e10_Rep5), color = factor(APY_e10_Rep5), alpha = factor(APY_e10_Rep5) ) ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + scale_alpha_manual(limits = c("0", "1", "2", "3", "4", "5"), values = c(0.2, 1, 1, 1, 1, 1)) + scale_color_manual(values = c("aquamarine", "darkblue", "darkgoldenrod", "darkred", "deeppink4", "black")) +
  labs(y = "", x = "", colour = "Core Status (Replicate)", alpha = "Core Status (Replicate)", shape = "Core Status (Replicate)", tag = "a)") +
  theme(legend.position = "top")


ind.coord_p3 = left_join(ind.coord_333, geno_refNaive1, by = c("row.names.res.umap.layout." = "geno_id"))
ind.coord_p3 = ind.coord_p3[order(ind.coord_p3$Diag_e10),]

p3 = ggplot(ind.coord_p3, aes(x = ind.coord_p3[, 1], y = ind.coord_p3[, 2], color = factor(Diag_e10), alpha = factor(Diag_e10) ) ) + 
  geom_point(aes(shape = factor(Diag_e10), color = factor(Diag_e10), alpha = factor(Diag_e10) ) ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + scale_alpha_manual(limits = c("0", "1"), values = c(0.2, 1)) + scale_color_manual(values = c("aquamarine", "darkblue")) +
  labs(y = "", x = "", colour = "Core Status", alpha = "Core Status", shape = "Core Status", tag = "b)") +
  theme(legend.position = "top")


ind.coord_p4 = left_join(ind.coord_333, geno_refNaive2, by = c("row.names.res.umap.layout." = "geno_id"))
ind.coord_p4 = ind.coord_p4[order(ind.coord_p4$Weight_e10_Rep1, ind.coord_p4$Weight_e10_Rep2, ind.coord_p4$Weight_e10_Rep3, ind.coord_p4$Weight_e10_Rep4, ind.coord_p4$Weight_e10_Rep5),]

ind.coord_p4 = ind.coord_p4 %>%
  dplyr::mutate(Weight_e10_Rep2 = ifelse(Weight_e10_Rep2 == 1, 2, NA))  %>%
  dplyr::mutate(Weight_e10_Rep3 = ifelse(Weight_e10_Rep3 == 1, 3, NA))  %>%
  dplyr::mutate(Weight_e10_Rep4 = ifelse(Weight_e10_Rep4 == 1, 4, NA))  %>%
  dplyr::mutate(Weight_e10_Rep5 = ifelse(Weight_e10_Rep5 == 1, 5, NA))

p4 = ggplot(ind.coord_p4, aes(x = ind.coord_p4[, 1], y = ind.coord_p4[, 2] ) ) + 
  geom_point(aes(shape = factor(Weight_e10_Rep1), color = factor(Weight_e10_Rep1), alpha = factor(Weight_e10_Rep1) ) ) +
  geom_point(aes(shape = factor(Weight_e10_Rep2), color = factor(Weight_e10_Rep2), alpha = factor(Weight_e10_Rep2) ) ) +
  geom_point(aes(shape = factor(Weight_e10_Rep3), color = factor(Weight_e10_Rep3), alpha = factor(Weight_e10_Rep3) ) ) +
  geom_point(aes(shape = factor(Weight_e10_Rep4), color = factor(Weight_e10_Rep4), alpha = factor(Weight_e10_Rep4) ) ) +
  geom_point(aes(shape = factor(Weight_e10_Rep5), color = factor(Weight_e10_Rep5), alpha = factor(Weight_e10_Rep5) ) ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + scale_alpha_manual(limits = c("0", "1", "2", "3", "4", "5"), values = c(0.2, 1, 1, 1, 1, 1)) + scale_color_manual(values = c("aquamarine", "darkblue", "darkgoldenrod", "darkred", "deeppink4", "black")) +
  labs(y = "", x = "", colour = "Core Status (Replicate)", alpha = "Core Status (Replicate)", shape = "Core Status (Replicate)", tag = "c)") +
  theme(legend.position = "top")


legend <- get_legend(p2)

panel_1 = grid.arrange(p2 + theme(legend.position="none"), 
                       p3 + theme(legend.position="none"), 
                       p4 + theme(legend.position="none"), 
                       p1 + theme(legend.position="none"),
                       legend,
                       nrow = 3, ncol=2,
                       layout_matrix = rbind(c(5,5), c(1,2), c(3,4)),
                       heights = c(3, 10, 10),
                       widths = c(10, 10)
                       )

ggsave(file = "Panel_umap_core.png", panel_1, height = PaperSize, width = PaperSize * 1.5, unit = "cm")



#### Venn Diagram: Core selection overlap between the methods ####

# geno_ref_Total = data.frame(geno_refPrior, geno_refAPY, geno_refNaive1, geno_refNaive2)

# table(geno_refPrior$Prior_e10 %in% geno_refAPY$APY_e10_Rep1) 
# geno_refPrior$Prior_e10 == geno_refAPY$APY_e10_Rep1

# pr98 = geno_refPrior %>% filter(Prior_e98==1) %>% dplyr::select(geno_id) %>% unlist()
# ap98 = geno_refAPY %>% filter(APY_e98_Rep1==1) %>% dplyr::select(geno_id) %>% unlist()
# table(pr98 %in% ap98)

# Compare Prior vs SVD Prior
pr98 = geno_refPrior %>% filter(Prior_e98==1) %>% dplyr::select(geno_id) %>% unlist()
svdpr98 = geno_refPriorSVD %>% filter(SVDprior_e98==1) %>% dplyr::select(geno_id) %>% unlist()
table(pr98 %in% svdpr98)

# Compare Prior vs Diagonal
diag98 = geno_refNaive1 %>% filter(Diag_e98==1) %>% dplyr::select(geno_id) %>% unlist()
table(pr98 %in% diag98)

# Compare APY replicates
ap198 = geno_refAPY %>% filter(APY_e98_Rep4==1) %>% dplyr::select(geno_id) %>% unlist()
ap298 = geno_refAPY %>% filter(APY_e98_Rep5==1) %>% dplyr::select(geno_id) %>% unlist()
table(ap198 %in% ap298)

# Compare APY replicates
wd198 = geno_refNaive2 %>% filter(Weight_e98_Rep2==1) %>% dplyr::select(geno_id) %>% unlist()
wd298 = geno_refNaive2 %>% filter(Weight_e98_Rep5==1) %>% dplyr::select(geno_id) %>% unlist()
table(wd198 %in% wd298)


venn.diagram(
  x = list(
    geno_refAPY %>% filter(APY_e98_Rep1==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refPrior %>% filter(Prior_e98==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refNaive1 %>% filter(Diag_e98==1) %>% dplyr::select(geno_id)  %>% unlist(),  
    geno_refNaive2 %>% filter(Weight_e98_Rep1==1) %>% dplyr::select(geno_id)  %>% unlist()
    ),
    category.names = c("Random", "Optimised", "Max diagonal", "Random diag weight"),
    filename = "Venn_e98.png",
    output = TRUE,
    cex = c(rep(1.5, 5), 3, rep(1.5, 9)),
    cat.cex = 2,
    cat.pos = c(-9, 4, 0, 0),
    lwd = 2,
    lty = 'blank',
    fill = color_roslin[1:4],
  imagetype="png" ,
  height = 20 , 
  width = 20 , 
  units = "cm",
  resolution = 600,
  compression = "lzw")

venn.diagram(
  x = list(
    geno_refAPY %>% filter(APY_e98_Rep1==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refAPY %>% filter(APY_e98_Rep2==1) %>% dplyr::select(geno_id)  %>% unlist(),  
    geno_refAPY %>% filter(APY_e98_Rep3==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refAPY %>% filter(APY_e98_Rep4==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refAPY %>% filter(APY_e98_Rep5==1) %>% dplyr::select(geno_id)  %>% unlist()
  ),
  category.names = c("1", "2", "3", "4", "5"),
  filename = "Venn_rnd_e98.png",
  output = TRUE,
  cex = 1.5,
  cat.cex = 2,
  lwd = 2,
  lty = 'blank',
  fill = color_roslin[1:5],
  imagetype="png" ,
  height = 20 , 
  width = 20 , 
  units = "cm",
  resolution = 600,
  compression = "lzw")
# print.mode = "percent"  

venn.diagram(
  x = list(
    geno_refNaive2 %>% filter(Weight_e98_Rep1==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refNaive2 %>% filter(Weight_e98_Rep2==1) %>% dplyr::select(geno_id)  %>% unlist(),  
    geno_refNaive2 %>% filter(Weight_e98_Rep3==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refNaive2 %>% filter(Weight_e98_Rep4==1) %>% dplyr::select(geno_id)  %>% unlist(),
    geno_refNaive2 %>% filter(Weight_e98_Rep5==1) %>% dplyr::select(geno_id)  %>% unlist()
  ),
  category.names = c("", "", "", "", ""),
  filename = "Venn_rnd_diag_e98.png",
  output = TRUE,
  lwd = 2,
  lty = 'blank',
  fill = color_roslin[1:5],
  imagetype="png" ,
  height = 20 , 
  width = 20 , 
  units = "cm",
  resolution = 600,
  compression = "lzw")

# print.mode = "percent"  


#### Test if the Prior method is selecting from each HS-family ####

# Take 3000 genotypes from population 20 that showed nice clustering with UMAP
# There are 50 clusters in gen 20, that equal to 50 paternal HS families
# Apply prior core selection with nCore = 50, and see if they will be spread nicely across the 50 clusters

gen_no_20 = phenotypes %>%
  filter(Generation == 20) 

ind.coord_2 = data.frame(ind.coord[, 1:2], factor(gen_no$Generation))

Z_20 = Z[row.names(Z) %in% gen_no_20$Aid,]

rm(geno_ref20)
geno_ref20 = geno_ref[,-3]
geno_ref20 = geno_ref20[geno_ref20$geno_id==gen_no_20$Aid, ]
  
knots_prior <- find_knots(Z_20, p = 50, method = "prior")
  
opticore <- knots_prior$data
index_core = NULL
# Had to change this, as now I have different dimesnsion of Z, so the order is not the 
# So, instead of doing it per order, do it per ID
# index_core$geno_order = opticore$knots
index_core$geno = row.names(opticore)
 
geno_ref2 = geno_ref %>%
  dplyr::mutate(core_status = if_else(geno_id %in% index_core$geno, true = 1, false = 0)) 

geno_ref20 = left_join(geno_ref20, geno_ref2, by = "geno_id")

# Select UMAP results for gen 20 only
ind.coord_g20 = ind.coord_333 %>%
  filter(ind.coord_333[,3]==20)

ind.coord_p1 = left_join(ind.coord_g20, geno_ref20, by = c("row.names.res.umap.layout." = "geno_id"))
ind.coord_p1 = ind.coord_p1[order(ind.coord_p1$core_status),]



pp = ggplot(ind.coord_4, aes(x = ind.coord_p1[, 1], y = ind.coord_p1[, 2], group = ind.coord_p1[, 3])) + 
  geom_point(aes(color = ind.coord_p1[, 3]), alpha = ind.coord_p1[, 3]) + 
  theme_classic() + labs(y = "", x = "", colour = "Generation", alpha = "Generation") + scale_color_manual(values = color_roslin[5]) +
  theme(legend.position = "top") + ylim(min_y, max_y) + xlim(min_x, max_x) 

p1 = ggplot(ind.coord_p1, aes(x = ind.coord_p1[, 1], y = ind.coord_p1[, 2], color = factor(core_status), alpha = factor(core_status) ) ) + 
  geom_point(aes(shape = factor(core_status), color = factor(core_status), alpha = factor(core_status) ) ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + scale_alpha_manual(limits = c("0", "1"), values = c(0.5, 1)) + scale_color_manual(values = c("aquamarine", "darkblue")) +
  labs(y = "", x = "", colour = "Core Status", alpha = "Core Status", shape = "Core Status") +
  theme(legend.position = "top") + ylim(min_y, max_y) + xlim(min_x, max_x) 

ggsave(plot = p1 + PaperTheme +
         theme(
           axis.text.x = element_blank(), 
           axis.text.y = element_blank(), 
           axis.ticks = element_blank()) ,
       filename = "UMAP_gen20_core.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm")

# Version for manuscript
ind.coord_4 = ind.coord_333 %>%
  filter(ind.coord_333[,3]==20)

p_umap20 = ggplot(ind.coord_4, aes(x = ind.coord_4[, 1], y = ind.coord_4[, 2])) +
  geom_point(aes(color = as.factor(ind.coord_4[, 3]), group = as.factor(ind.coord_4[, 3]), shape = as.factor(ind.coord_4[, 3]))) +
  theme_classic(base_size = 12, base_family = "Times New Roman") + labs(y = "", x = "", colour = "Generation", shape = "Generation") + 
  # scale_color_manual(values = "lightgrey") + 
  # scale_shape_manual(values = c(16)) +
  theme(legend.position = "top")

ind.coord_p2 = ind.coord_p1 %>%
  filter(core_status == 1)

p_50sim = p_umap20 + geom_point(data = ind.coord_p2, mapping = aes(x = ind.coord_p2[, 1], y = ind.coord_p2[, 2], shape = "Optimised", color = "Optimised"), size = c(2.5) ) + 
  guides(size = FALSE) + 
  scale_color_manual(values = c("lightgrey", color_roslin[5])) + 
  scale_shape_manual(values = c(16, 4))

ggsave(file = "p_50sim.png", p_50sim, height = PaperSize, width = PaperSize * 1.5, unit = "cm")



#### Evaluate re-ranking from multiple evaluation rounds ####

rm(rank_res)
tmpX = tmp_ac %>%
  dplyr::select("Full", "Aid")

rank_top = tmpX %>%
  dplyr::slice_max(order_by = TBV, prop = 0.01, with_ties = FALSE) %>%
  dplyr::select("Aid") 
colnames(rank_top) ="top_EBV" 

rank_min = tmpX %>%
  dplyr::slice_min(order_by = TBV, prop = 0.01, with_ties = FALSE) %>%
  dplyr::select("Aid") 
colnames(rank_top) ="top_EBV" 

rank_res = cbind(rank_top, rank_min)



# Calculate re-ranking for Prior runs
tmpX = tmp_ac %>%
  dplyr::select(starts_with("Prior"), "Aid")


for(k in 1:(length(tmpX)-1)){
  rank_top = tmpX %>%
    dplyr::slice_max(order_by = tmpX[k], prop = 0.01, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  # dplyr::rename( (colnames(tmpX[k]) ) = Aid)
  colnames(rank_top) = paste0("top_", colnames(tmpX[k]) ) 
  
  rank_min = tmpX %>%
    dplyr::slice_min(order_by = tmpX[k], prop = 0.01, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_min) = paste0("min_", colnames(tmpX[k]) ) 
  
  rank_res_new = cbind(rank_top, rank_min)    
  rank_res = cbind(rank_res, rank_res_new)
}




# Calculate re-ranking for APY runs
cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98")
# for(x in 1:length(cnm)){
  
  tmpX = tmp_ac %>%
    dplyr::select(starts_with(cnm[7]), "Aid")
  
  for(k in 1:(length(tmpX)-1)){
    rank_top = tmpX %>%
      dplyr::slice_max(order_by = tmpX[k], prop = 0.01, with_ties = FALSE) %>%
      dplyr::select("Aid") 
      # dplyr::rename( (colnames(tmpX[k]) ) = Aid)
      colnames(rank_top) = paste0("top_", colnames(tmpX[k]) ) 

    rank_min = tmpX %>%
      dplyr::slice_min(order_by = tmpX[k], prop = 0.01, with_ties = FALSE) %>%
      dplyr::select("Aid") 
      colnames(rank_min) = paste0("min_", colnames(tmpX[k]) ) 
    
    rank_res_new = cbind(rank_top, rank_min)    
    rank_res = cbind(rank_res, rank_res_new)
  }

  
  
  
  table(rank_res$top_EBV %in% rank_res$top_Prior_e98)
  table(rank_res$top_EBV %in% rank_res$top_APY_e98_Rep1)
  table(rank_res$top_EBV %in% rank_res$top_APY_e98_Rep2)
  table(rank_res$top_EBV %in% rank_res$top_APY_e98_Rep3)
  table(rank_res$top_EBV %in% rank_res$top_APY_e98_Rep4)
  table(rank_res$top_EBV %in% rank_res$top_APY_e98_Rep5)
  
  table(rank_res$top_Prior_e98 %in% rank_res$top_APY_e98_Rep1)
  table(rank_res$top_Prior_e98 %in% rank_res$top_APY_e98_Rep2)
  table(rank_res$top_Prior_e98 %in% rank_res$top_APY_e98_Rep3)
  table(rank_res$top_Prior_e98 %in% rank_res$top_APY_e98_Rep4)
  table(rank_res$top_Prior_e98 %in% rank_res$top_APY_e98_Rep5)
  
  table(rank_res$top_APY_e98_Rep1 %in% rank_res$top_APY_e98_Rep2)
  table(rank_res$top_APY_e98_Rep1 %in% rank_res$top_APY_e98_Rep3)
  table(rank_res$top_APY_e98_Rep1 %in% rank_res$top_APY_e98_Rep4)
  table(rank_res$top_APY_e98_Rep1 %in% rank_res$top_APY_e98_Rep5)
  table(rank_res$top_APY_e98_Rep2 %in% rank_res$top_APY_e98_Rep3)
  table(rank_res$top_APY_e98_Rep2 %in% rank_res$top_APY_e98_Rep4)
  table(rank_res$top_APY_e98_Rep2 %in% rank_res$top_APY_e98_Rep5)
  table(rank_res$top_APY_e98_Rep3 %in% rank_res$top_APY_e98_Rep4)
  table(rank_res$top_APY_e98_Rep3 %in% rank_res$top_APY_e98_Rep5)
  table(rank_res$top_APY_e98_Rep4 %in% rank_res$top_APY_e98_Rep5)
  
  

  
  
  venn.diagram(
    x = list(
      rank_res$top_APY_e98_Rep1,
      rank_res$top_APY_e98_Rep2,
      rank_res$top_APY_e98_Rep3,
      rank_res$top_APY_e98_Rep4,
      rank_res$top_APY_e98_Rep5
    ),
    category.names = c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5"),
    filename = "Venn_reranking_e10.png",
    output = TRUE,
    lwd = 2,
    lty = 'blank',
    fill = color_roslin[1:5],
    imagetype="png" ,
    height = 20 , 
    width = 20 , 
    units = "cm",
    resolution = 600,
    compression = "lzw")
  


  
  #### Re-Ranking Spearman's rank correlations ####
  
  # Do it for:
  ## Full data 
  ## Validation animals only (Generation 20)
  ## And repeat for correlations with TBV, and with Full
  tmp_ac = final_data6 %>%
    filter(Generation == 20)
  
  rezultati_sp = as_tibble(colnames(n_final))
  rezultati_sp$nvar = as.numeric(str_sub(rezultati_sp$value, 2))  
  rezultati_sp$ncore = n_final[1,]
  # Calculate accuracy for Direct inverse
  rezultati_sp$Full = rep(round(cor(tmp_ac$TBV, tmp_ac$Full, method = "pearson"), digits = 2), each = length(n_final))
  
  # Note: Create a function for calculating this accuracies - current code is messy!!!
  # For correlations(GEBV_FULL, GEBV_APY) change "TBV" to "Full" in the code 
  
  # Calculate accuracy for APY runs
  cnm = c("APY_e10", "APY_e30", "APY_e50", "APY_e70", "APY_e90", "APY_e95", "APY_e98")
  rez_apy = data.frame(matrix(ncol = 3, nrow = 0))
  for(x in 1:length(cnm)){
    
    tmpX = tmp_ac %>%
      dplyr::select(starts_with(cnm[x]), "TBV")
    
    rescor = NULL
    rescor_new = NULL
    for(k in 1:(length(tmpX)-1)){
      rescor_new = round(cor(tmpX$TBV, tmpX[k], method = "spearman"), digits = 2)
      rownames(rescor_new) = colnames(tmpX)[k]
      rescor = rbind(rescor, rescor_new)
    }
    
    rescor_mean = round(mean(rescor), 2)
    rescor_L = round(quantile(rescor, na.rm=TRUE, probs=0.025), 2)
    rescor_U = round(quantile(rescor, na.rm=TRUE, probs=0.975), 2)
    
    rez_apy = rbind(rez_apy, c(rescor_mean, rescor_L, rescor_U))
    rownames(rez_apy)[x] = cnm[x]
    
  }
  colnames(rez_apy) = c("APY", "APY_L", "APY_U")
  
  rezultati_sp = bind_cols(rezultati_sp, rez_apy)
  
  
  # round(cor(as.matrix(cbind(tmp_ac$APY_e98_Rep1, tmp_ac$APY_e98_Rep2, tmp_ac$APY_e98_Rep3, tmp_ac$APY_e98_Rep4, tmp_ac$APY_e98_Rep5)), use = "pairwise.complete.obs"), digits = 2)
  
  
  # Calculate accuracy for Prior runs
  tmpX = tmp_ac %>%
    dplyr::select(starts_with("Prior"), "TBV")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = round(cor(tmpX$TBV, tmpX[k], method = "spearman"), digits = 2)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  rescor = as.data.frame(rescor)
  colnames(rescor) = "Prior"
  rezultati_sp = bind_cols(rezultati_sp, rescor)
  
  # Calculate accuracy for Diag runs
  tmpX = tmp_ac %>%
    dplyr::select(starts_with("Diag"), "TBV")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = round(cor(tmpX$TBV, tmpX[k], method = "spearman"), digits = 2)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  rescor = as.data.frame(rescor)
  colnames(rescor) = "Diag"
  rezultati_sp = bind_cols(rezultati_sp, rescor)
  
  # Calculate accuracy for Weighted Diag runs
  cnm = c("Weight_e10", "Weight_e30", "Weight_e50", "Weight_e70", "Weight_e90", "Weight_e95", "Weight_e98")
  rez_wgt = data.frame(matrix(ncol = 3, nrow = 0))
  for(x in 1:length(cnm)){
    
    tmpX = tmp_ac %>%
      dplyr::select(starts_with(cnm[x]), "TBV")
    
    rescor = NULL
    rescor_new = NULL
    for(k in 1:(length(tmpX)-1)){
      rescor_new = round(cor(tmpX$TBV, tmpX[k], method = "spearman"), digits = 2)
      rownames(rescor_new) = colnames(tmpX)[k]
      rescor = rbind(rescor, rescor_new)
    }
    
    rescor_mean = round(mean(rescor), 2)
    rescor_L = round(quantile(rescor, na.rm=TRUE, probs=0.025), 2)
    rescor_U = round(quantile(rescor, na.rm=TRUE, probs=0.975), 2)
    rez_wgt = rbind(rez_wgt, c(rescor_mean, rescor_L, rescor_U))
    rownames(rez_wgt)[x] = cnm[x]
    
  }
  colnames(rez_wgt) = c("WGT", "WGT_L", "WGT_U")
  
  rezultati_sp = bind_cols(rezultati_sp, rez_wgt)
  
  # round(cor(as.matrix(cbind(tmp_ac$Weight_e98_Rep1, tmp_ac$Weight_e98_Rep2, tmp_ac$Weight_e98_Rep3, tmp_ac$Weight_e98_Rep4, tmp_ac$Weight_e98_Rep5)), use = "pairwise.complete.obs"), digits = 2)
  
  
  
  # Calculate accuracy for SVD Prior runs
  tmpX = tmp_ac %>%
    dplyr::select(starts_with("SVD"), "TBV")
  
  rescor = NULL
  rescor_new = NULL
  for(k in 1:(length(tmpX)-1)){
    rescor_new = round(cor(tmpX$TBV, tmpX[k], use = method = "spearman"), digits = 2)
    rownames(rescor_new) = colnames(tmpX)[k]
    rescor = rbind(rescor, rescor_new)
  }
  rescor = as.data.frame(rescor)
  colnames(rescor) = "SVDPrior"
  rezultati_sp = bind_cols(rezultati_sp, rescor)
  
  
 #### Evaluate sire re-ranking ####
  

tmp_sires = final_data6 %>% 
    filter(Generation == 20 & Sex == "M")
    
  rank_top = tmp_sires %>%
    dplyr::slice_max(order_by = TBVSTD, n = 45, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_top) ="top_TBV" 
  
  rank_full = tmp_sires %>%
    dplyr::slice_max(order_by = Full, n = 45, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_full) ="top_Full" 
  
  rank_rnd = tmp_sires %>%
    dplyr::slice_max(order_by = APY_e98_Rep5, n = 45, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_rnd) ="top_Rnd" 
  
  rank_prior = tmp_sires %>%
    dplyr::slice_max(order_by = Prior_e98, n = 45, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_prior) ="top_Prior" 
  
  
  rank_diag = tmp_sires %>%
    dplyr::slice_max(order_by = Diag_e98, n = 45, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_diag) ="top_Diag" 
  
  rank_weight = tmp_sires %>%
    dplyr::slice_max(order_by = Weight_e98_Rep5, n = 45, with_ties = FALSE) %>%
    dplyr::select("Aid") 
  colnames(rank_weight) ="top_Weight" 
  
  
  rank_sires = cbind(rank_top, rank_full, rank_rnd, rank_prior, rank_diag, rank_weight)
  
  table(rank_sires$top_TBV %in% rank_sires$top_Diag)
  table(rank_sires$top_TBV %in% rank_sires$top_Weight)
  table(rank_sires$top_TBV %in% rank_sires$top_Prior)
  table(rank_sires$top_TBV %in% rank_sires$top_Rnd)
  table(rank_sires$top_TBV %in% rank_sires$top_Full)
  
  rank_sires = cbind(rank_full, rank_rnd, rank_prior, rank_diag, rank_weight)
  
  table(rank_sires$top_Full %in% rank_sires$top_Diag)
  table(rank_sires$top_Full %in% rank_sires$top_Weight)
  table(rank_sires$top_Full %in% rank_sires$top_Prior)
  table(rank_sires$top_Full %in% rank_sires$top_Rnd)

#### Absolute difference and MSE ####
  
tmp_ac = final_data6 %>%
    filter(Generation == 20) %>%
    dplyr::select(Aid, Sid, Did, Sex, Generation, TBV, TBVSTD, Full, 
           starts_with("APY"), starts_with("Prior"), starts_with("Diag"), starts_with("Weight") )
  
  # Manually repeat for "Full" instead of "TBVSTD"
  tmp_diff = tmp_ac    
  for(k in 8:(length(tmp_ac))){
    absdiff = abs(tmp_ac$TBVSTD - tmp_ac[k])
    colnames(absdiff) = paste0("AbsDiff_", colnames(absdiff))
    tmp_diff = cbind(tmp_diff, absdiff)
  }
  
  # Manually repeat for different core sizes, e.g., e_10, ...
  round(mean(tmp_diff$AbsDiff_Full), digits = 2)
  round(mean(tmp_diff$AbsDiff_APY_e98_Rep1), digits = 2)
  round(mean(tmp_diff$AbsDiff_Weight_e98_Rep1), digits = 2)
  round(mean(tmp_diff$AbsDiff_Diag_e98), digits = 2)
  round(mean(tmp_diff$AbsDiff_Prior_e98), digits = 2)
  
  round(max(tmp_diff$AbsDiff_Full), digits = 4)
  round(max(tmp_diff$AbsDiff_APY_e98_Rep1), digits = 4)
  round(max(tmp_diff$AbsDiff_Weight_e98_Rep1), digits = 4)
  round(max(tmp_diff$AbsDiff_Diag_e98), digits = 4)
  round(max(tmp_diff$AbsDiff_Prior_e98), digits = 4)
  
  ggplot(tmp_diff, aes(x = "Prior", y = AbsDiff_Prior_e98)) + geom_point(size = 1, alpha = 0.2) + 
    geom_point(aes(x = "Prior", y = mean(AbsDiff_Prior_e98)), size = 2, alpha = 1, colour = "red") +
    geom_point(aes(x = "Prior", y = max(AbsDiff_Prior_e98)), size = 2, alpha = 1, colour = "blue") +
    geom_point(aes(x = "Random", y = AbsDiff_APY_e98_Rep1), size = 1, alpha = 0.2) +
    geom_point(aes(x = "Random", y = mean(AbsDiff_APY_e98_Rep1)), size = 2, alpha = 1, colour = "red") +
    geom_point(aes(x = "Random", y = max(AbsDiff_APY_e98_Rep1)), size = 2, alpha = 1, colour = "blue") +
    geom_point(aes(x = "Weight", y = AbsDiff_Weight_e98_Rep1), size = 1, alpha = 0.2) +
    geom_point(aes(x = "Weight", y = mean(AbsDiff_Weight_e98_Rep1)), size = 2, alpha = 1, colour = "red") +
    geom_point(aes(x = "Weight", y = max(AbsDiff_Weight_e98_Rep1)), size = 2, alpha = 1, colour = "blue") +
    geom_point(aes(x = "Diag", y = AbsDiff_Diag_e98), size = 1, alpha = 0.2) +
    geom_point(aes(x = "Diag", y = mean(AbsDiff_Diag_e98)), size = 2, alpha = 1, colour = "red") +
    geom_point(aes(x = "Diag", y = max(AbsDiff_Diag_e98)), size = 2, alpha = 1, colour = "blue") +
    theme_classic() +
    labs(y = "Absolute difference", x = "Method") 
    
  # MSE  
  # Manually repeat for "Full" instead of "TBVSTD"
  tmp_diff = tmp_ac    
  for(k in 8:(length(tmp_ac))){
    absdiff = (tmp_ac$Full- tmp_ac[k])^2
    colnames(absdiff) = paste0("AbsDiff_", colnames(absdiff))
    tmp_diff = cbind(tmp_diff, absdiff)
  }
  
  # Manually repeat for different core sizes, e.g., e_10, ... 
  round(mean(tmp_diff$AbsDiff_Full), digits = 4)
  round(mean(tmp_diff$AbsDiff_APY_e98_Rep1), digits = 4)
  round(mean(tmp_diff$AbsDiff_Weight_e98_Rep1), digits = 4)
  round(mean(tmp_diff$AbsDiff_Diag_e98), digits = 4)
  round(mean(tmp_diff$AbsDiff_Prior_e98), digits = 4)

  