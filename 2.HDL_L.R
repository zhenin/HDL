# 2. Run HDL-L
library(dplyr)

# Define file paths
base_path <- "/Your/Path"
sumstat_path <- file.path(base_path, "sumstat")
result_path <- file.path(base_path, "hdll")
bedp_path <- file.path(base_path, "bfile")
LD_path <- base_path  # Assuming LD matrices are in base_path, adjust as needed

# Create result directory if it doesn't exist
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

# Load required data
load(file.path(base_path, "HDLL_LOC_snps.RData"))  # Loads 'NEWLOC'

# Filter to include only chromosome 22 (adjust as needed)
NEWLOC <- NEWLOC %>% filter(CHR == 22)

# Define HDL-L analysis function
run_hdl_l <- function(sim) {
  # Load sumstats
  gwas1.df <- fread(file.path(sumstat_path, sprintf("y1.sim.%d.sum.stats.txt", sim)))
  gwas2.df <- fread(file.path(sumstat_path, sprintf("y2.sim.%d.sum.stats.txt", sim)))
  res_HDL_list <- list()
  for (i in 1:nrow(NEWLOC)) {
    chr <- NEWLOC$CHR[i]
    piece <- NEWLOC$piece[i]
    try({
      res <- HDL.L(
        gwas1 = gwas1.df,
        gwas2 = gwas2.df,
        Trait1name = "gwas1",
        Trait2name = "gwas2",
        LD.path = LD_path,
        bim.path = bedp_path,
        chr = chr,
        piece = piece,
        N0 = 0
      )
      res_HDL_list[[length(res_HDL_list) + 1]] <- res
    }, silent = TRUE)
  }
  res_HDL <- do.call(rbind, res_HDL_list)
  save(res_HDL, file = file.path(result_path, sprintf("result_HDL.Lsim%d.RData", sim)))
}

# Run HDL-L for each simulation
num_simulations <- 100  # Adjust if needed
for (sim in 1:num_simulations) {
  cat("Running HDL-L for simulation", sim, "\n")
  run_hdl_l(sim)
}