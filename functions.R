# Functions for "Optimisation of the core subset for the APY approximation of genomic relationships"
# Pocrnic, Lindgren, Tolhurst, Herring & Gorjanc (2022) 
# DOI: 

#### Functions ####

# A: AA' is the covariance
# knots: already selected knots
# B: BB' is the conditional covariance given already selected knots
# p: The number of new knots to choose
# method:
#   prior: Maximises prior variance for the knot in each step
#   posterior_max: Minimises the maximum posterior variance in each step
#   posterior_avg: Minimises the average posterior variance in each step
#     (to avoid issues when the maximum isn't improved)
find_knots <- function(A,
                       knots = NULL,
                       B = NULL,
                       p,
                       method = c("prior", "posterior_max", "posterior_avg")) {
  method <- match.arg(method)
  if (is.null(knots)) {
    knots <- c()
    B <- A
  }
  if (is.null(B)) {
    B <- A
    for (i in knots) {
      b <- B[i, , drop = FALSE]
      B <- B - (B %*% t(b)) %*% (b / sum(b^2))
    }
  }
  target <- c()
  total_var <- c()
  for (k in seq_len(p)) {
    v_prior <- rowSums(B^2)
    if (method == "prior") {
      i <- which.max(v_prior)
      target <- c(target, v_prior[i])
    } else if (method == "posterior_max") {
      ok <- v_prior > 0
      v_post <- rep(Inf, nrow(A))
      v_post[ok] <-
        vapply(which(ok), function(i) {
          b <- B[i, , drop = FALSE]
          max(rowSums((B - (B %*% t(b)) %*% (b / sum(b^2)))^2))
        }, 1.0)
      # Make sure we don't pick an existing knot by accident
      v_post[knots] <- Inf
      i <- which.min(v_post)
      target <- c(target, v_post[i])
    } else if (method == "posterior_avg") {
      ok <- v_prior > 0
      v_post <- rep(Inf, nrow(A))
      v_post[ok] <-
        vapply(which(ok), function(i) {
          b <- B[i, , drop = FALSE]
          mean(rowSums((B - (B %*% t(b)) %*% (b / sum(b^2)))^2))
        }, 1.0)
      # Make sure we don't pick an existing knot by accident
      v_post[knots] <- Inf
      i <- which.min(v_post)
      target <- c(target, v_post[i])
    }
    knots <- c(knots, i)
    b <- B[i, , drop = FALSE]
    B <- B - (B %*% t(b)) %*% (b / sum(b^2))
    # B <- B - (tcrossprod(B, b)) %*% (b / sum(b^2))
    total_var <- c(total_var, sum(B^2)) # sum of the remaining variances
  }
  list(data = data.frame(k = seq_len(p),
                         knots = knots,
                         total_var = total_var,
                         target = target,
                         method = method),
       B = B)
}

# Aold: Aold'Aold is the covariance between existing points
# knots: Existing knots in Aold
# Anew: Rows to be added to an existing A-matrix
# Returns Bnew, the rows of the expanded B-matrix such that with
#   B = rbind(Bold, Bnew), the matrix BB' is the covariance
#   conditionally on the knots.
#
# Example:
#   # Find 5 knots in A1:
#   info1 <- find_knots(A = A1, p = 5)
#   # Find another 5 knots in rbind(A1, A2):
#   info2 <- find_knots(
#     A = rbind(A1, A2), # This A won't actually be used
#     knots = info$data$knots,
#     B = rbind(info$B,
#               expand_B(A1, info$data$knots, A2)),
#     p = 5)
expand_B <- function(Aold, knots, Anew) {
  B <- rbind(Aold[knots, , drop = FALSE], Anew)
  for (i in 1:length(knots)) {
    b <- B[i, , drop = FALSE]
    B <- B - (B %*% t(b)) %*% (b / sum(b^2))
  }
  Bnew <- B[length(knots) + seq_len(nrow(Anew)), , drop = FALSE]
  Bnew
}


# APY inverse (result is "Full" APY G inverse matrix)
# Inputs are GRM and vector indicating rows/columns in GRM corresponding to the core animals
# For test: fG and geno_ref2

APY_inverse <- function(GRM, corelist) {
  # Index core and non-core
  core = corelist$core_status == 1
  noncore = corelist$core_status == 0
  
  ncore = sum(core)
  nnoncore = sum(noncore)
  
  # Partition G to core and non-core
  Gcc = GRM[core, core]
  Gcn = GRM[core, noncore]
  Gnc = GRM[noncore, core]
  Gnn = GRM[noncore, noncore]
  
  # Gpart = rbind(cbind(Gcc, Gcn), cbind(Gnc, Gnn)) 
  
  # APY inverse based on above Gpart:
  Gcc_inv = solve(Gcc)
  
  Mnn_inv = matrix(data = 0, nrow = nnoncore, ncol = nnoncore)
  for(i in 1:nnoncore) 
  { 
    Mnn_inv[i,i] <- 1 / (Gnn[i,i] - Gnc[i,] %*% Gcc_inv %*% Gcn[,i]) 
  } 
  
  APY11 = Gcc_inv + Gcc_inv %*% Gcn %*% Mnn_inv %*% Gnc %*% Gcc_inv 
  APY12 = -1 * Gcc_inv %*% Gcn %*% Mnn_inv 
  APY21 = -1 * Mnn_inv %*% Gnc %*% Gcc_inv 
  # APY21 = t(APY12)
  APY22 = Mnn_inv 
  
  Ginv_APY = rbind(cbind(APY11, APY12), 
                   cbind(APY21, APY22))
  
  # Create index of ID's
  # Note - order is not the same as in the original G
  apy_ref = tibble(geno_id = c(rownames(Gcc), rownames(Gnn)), apy_order = seq(1:nrow(Ginv_APY)))
  orig_ref = full_join(corelist, apy_ref, by = "geno_id")
  
  list(APY_Ginv = Ginv_APY, index_file = orig_ref)
}


# Function to record data:
data_rec <- function(datafile, popname) {
  datafile = rbind(datafile,
                   tibble(Aid        = popname@id,
                          Sid        = popname@father,
                          Did        = popname@mother,
                          Trait      = popname@pheno,
                          Sex        = popname@sex,
                          TBV        = popname@gv,
                          Generation = Generation))
}


# Function for running GBLUP via BLUPF90 and reading the solutions with given name (e.g., EBV_colname = "Full")
run_gblupf90 <- function(datafile, EBV_colname, method = c("Regular", "APY")) { 
  # Run:
  system(command = "echo blupf90.par | $HOME/bin/blupf90 | tee blup.log")
  
  # Import blupf90 solutions to R: 
  sol1 <- read_table2("solutions", col_names = FALSE, skip = 1,
                      col_types = cols(.default = col_double(),
                                       X1 = col_double(),
                                       X2 = col_double(),
                                       X3 = col_double(),
                                       X4 = col_double(),
                                       X5 = col_double()))
  colnames(sol1) = c("Trait", "Effect", "Level", EBV_colname, paste0("SE_", EBV_colname))
  
  # Extract EBV and SE from the file (animal effect was effect #2)
  tmp = sol1 %>%
    filter(Trait == 1 & Effect == 2) %>%
    dplyr::select("Level", EBV_colname, paste0("SE_", EBV_colname))
  
  # Connect EBV with original animal ID
  
  if (method == "Regular") {
    tmp2 = dplyr::full_join(datafile, tmp, by = c("geno_order" = "Level")) %>%
      dplyr::select("Aid", EBV_colname, paste0("SE_", EBV_colname))  
  } else if (method == "APY") {
    tmp2 = full_join(datafile, tmp, by = c("apy_order" = "Level")) %>%
      dplyr::select("Aid", EBV_colname, paste0("SE_", EBV_colname)) 
  }
  
  return(tmp2)
}


prepare_par <- function() {
sink("blupf90.par", type="output")
writeLines("#blupf90 parametar file
DATAFILE
 Blupf90.dat
NUMBER_OF_TRAITS
 1
NUMBER_OF_EFFECTS
 2
OBSERVATION(S)
 2
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
 4      1 cross 
 3      15000 cross 
RANDOM_RESIDUAL VALUES
 2.3333    
RANDOM_GROUP
 2
RANDOM_TYPE
 user_file
FILE
 apyinv.txt                                                                   
(CO)VARIANCES
 1.0000
OPTION msg 1 
OPTION sol se
OPTION use_yams
")
sink()
}




