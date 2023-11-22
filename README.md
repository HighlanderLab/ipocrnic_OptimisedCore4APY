# Optimisation of the core subset for the APY approximation of genomic relationships
## Pocrnic, Lindgren, Tolhurst, Herring & Gorjanc (2022)

Published in GSE: [Optimisation of the core subset for the APY approximation of genomic relationships](https://doi.org/10.1186/s12711-022-00767-x)

- `functions.R` includes (1) `find_knots` function to run conditional sampling core selection algorithms presented in the manuscript; (2) `expand_B` extension function; (3) `APY_inverse` brute force APY function (in the case you use free version of `BLUPF90` that doesn't have `OPTION apy ...` available); (4) `data_rec` function to collect data from `AlphaSimR` simulation; (5) `run_gblupf90` function to run BLUPF90 and read solutions into `R`; and (6) `prepare_par` function to prepare parameter file for `BLUPF90`. 

- You might need to change the following line of code in `functions.R` based on your OS and location of the `BLUPF90` binary;
```
system(command = "echo blupf90.par | $HOME/bin/blupf90 | tee blup.log"))
```

- `OptimisedCore4APY.R` includes the code to simulate the data using `AlphaSimR` and to generate core selection scenarios as in the manuscript. 

- Please note that `BLUPF90` program (Ignacy Misztal et al., University of Georgia, USA) has to be installed on your system. For download and installation details, please consult [BLUPF90 website](http://nce.ads.uga.edu/software/). 
