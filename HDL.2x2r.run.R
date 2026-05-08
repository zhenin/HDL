args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("HDL.2x2r demo: BMI (21001) x Asthma-diagnosed (6152_8)\n")
  cat("  Trait X : BMI (UKB field 21001)\n")
  cat("  Trait Y : Asthma-diagnosed (UKB field 6152_8)\n")
  cat("  Sample 1: UKB males\n")
  cat("  Sample 2: UKB females\n\n")
  cat("Usage:\n")
  cat("  Rscript HDL.2x2r.run.R LD.path=<path> [output.file=<file>] [Nref=<n>] [N0.11=<n>] [N0.22=<n>]\n\n")
  cat("Required:\n")
  cat("  LD.path=     Path to the LD reference panel directory\n")
  cat("               (UKB_imputed_SVD_eigen99_extraction; see https://github.com/zhenin/HDL/wiki/Reference-panels)\n\n")
  cat("Optional:\n")
  cat("  output.file= Path for the result file (default: HDL.2x2r_BMI_Asthma6152_result.txt)\n")
  cat("  Nref=        LD reference panel sample size (default: 335265)\n")
  cat("  N0.11=       Sample overlap for gwasX1 vs gwasY1 (default: min(N_X1, N_Y1))\n")
  cat("  N0.22=       Sample overlap for gwasX2 vs gwasY2 (default: min(N_X2, N_Y2))\n\n")
  cat("Example:\n")
  cat("  Rscript HDL.2x2r.run.R \\\n")
  cat("    LD.path=/path/to/UKB_imputed_SVD_eigen99_extraction \\\n")
  cat("    output.file=HDL.2x2r_BMI_Asthma6152_result.txt\n")
  q()
}

## ---- parse arguments -------------------------------------------------------- ##
.arg <- function(pattern, args) {
  x <- gsub(pattern, "", args[grep(pattern, args)])
  if (length(x) == 0) NULL else x
}

LD.path     <- .arg("LD.path=", args)
output.file <- .arg("output.file=", args)
Nref        <- .arg("Nref=", args)
N0.11       <- .arg("N0.11=", args)
N0.22       <- .arg("N0.22=", args)

if (is.null(LD.path)) stop("LD.path is required. Run without arguments to see usage.")

output.file <- ifelse(is.null(output.file), "HDL.2x2r_BMI_Asthma6152_result.txt", output.file)
Nref        <- ifelse(is.null(Nref), 335265, as.numeric(Nref))
N0.11       <- if (is.null(N0.11)) NULL else as.numeric(N0.11)
N0.22       <- if (is.null(N0.22)) NULL else as.numeric(N0.22)

args.print <- paste("Function arguments:", paste(args, collapse = "\n"), sep = "\n")
cat(args.print, "\n\n")

if (output.file != "" && file.exists(output.file)) file.remove(output.file)

## ---- locate bundled example data -------------------------------------------- ##
script.path <- normalizePath(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
data.dir    <- file.path(dirname(script.path), "HDL", "data")

## ---- load example summary statistics --------------------------------------- ##
cat("Loading example summary statistics...\n")
gwasX1 <- readRDS(file.path(data.dir, "21001_male.hdl.rds"))
gwasY1 <- readRDS(file.path(data.dir, "6152_8_male.hdl.rds"))
gwasX2 <- readRDS(file.path(data.dir, "21001_female.hdl.rds"))
gwasY2 <- readRDS(file.path(data.dir, "6152_8_female.hdl.rds"))

cat(sprintf("  gwasX1 (BMI,    male)   : N = %d, SNPs = %d\n", unique(gwasX1$N), nrow(gwasX1)))
cat(sprintf("  gwasY1 (Asthma, male)   : N = %d, SNPs = %d\n", unique(gwasY1$N), nrow(gwasY1)))
cat(sprintf("  gwasX2 (BMI,    female) : N = %d, SNPs = %d\n", unique(gwasX2$N), nrow(gwasX2)))
cat(sprintf("  gwasY2 (Asthma, female) : N = %d, SNPs = %d\n", unique(gwasY2$N), nrow(gwasY2)))

## ---- run HDL.2x2r ----------------------------------------------------------- ##
library(HDL)

res <- HDL.2x2r(
  gwasX1.df   = gwasX1,
  gwasY1.df   = gwasY1,
  gwasX2.df   = gwasX2,
  gwasY2.df   = gwasY2,
  LD.path     = LD.path,
  Nref        = Nref,
  N0.11       = N0.11,
  N0.22       = N0.22,
  output.file = output.file
)

## ---- prepend arguments to output file --------------------------------------- ##
if (output.file != "" && file.exists(output.file)) {
  lines <- readLines(output.file)
  writeLines(c(args.print, "", lines), con = output.file)
}
