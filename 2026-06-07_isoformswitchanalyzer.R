###############################################
## IsoformSwitchAnalyzeR Full Pipeline
## Annotated version for Jeff
###############################################

library(IsoformSwitchAnalyzeR)
library(tximport)
library(readr)
library(dplyr)

###############################################
## 1. Define input files
###############################################

# Directory containing Salmon quant directories
salmon_dir <- "/path/to/salmon_quant/"

# Sample metadata (must contain: sampleID, condition)
samples <- read.csv("/path/to/sample_metadata.csv")

# Transcript-to-gene mapping (from GTF or txdb)
gtf_file <- "/path/to/annotation.gtf"

###############################################
## 2. Import quantifications with tximport
###############################################

# Build vector of quant.sf paths
files <- file.path(salmon_dir, samples$sampleID, "quant.sf")
names(files) <- samples$sampleID

txi <- tximport(files, type="salmon", txOut=TRUE)

###############################################
## 3. Create switchAnalyzeRlist object
###############################################

switchList <- importRdata(
  isoformCountMatrix   = txi$counts,
  isoformRepExpression = txi$abundance,
  designMatrix         = samples,
  isoformExonAnno      = gtf_file,
  isoformNtFasta       = "/path/to/transcripts.fa",
  isoformProteinFasta  = "/path/to/proteins.fa",
  showProgress         = TRUE
)

###############################################
## 4. Prefilter low-expression isoforms
###############################################

switchList <- preFilter(
  switchList,
  geneExpressionCutoff     = 3,
  isoformExpressionCutoff  = 3,
  removeSingleIsoformGenes = TRUE
)

###############################################
## 5. Test for differential isoform usage
###############################################

switchList <- isoformSwitchTestDEXSeq(
  switchList,
  reduceToSwitchingGenes = TRUE,
  alpha = 0.05
)

###############################################
## 6. Test for differential expression (optional)
###############################################

switchList <- isoformSwitchTestDRIMSeq(
  switchList,
  reduceToSwitchingGenes = TRUE
)

###############################################
## 7. Annotate ORFs
###############################################

switchList <- analyzeORF(
  switchList,
  pathToOutput = "ORF_output/",
  overwrite = TRUE
)

###############################################
## 8. Annotate protein domains (Pfam)
###############################################

switchList <- analyzePFAM(
  switchList,
  pathToPFAMresultFile = "/path/to/pfam_scan_results.txt",
  showProgress = TRUE
)

###############################################
## 9. Annotate signal peptides
###############################################

switchList <- analyzeSignalP(
  switchList,
  pathToSignalP = "/usr/bin/signalp"
)

###############################################
## 10. Annotate coding potential
###############################################

switchList <- analyzeCPC2(
  switchList,
  pathToCPC2resultFile = "/path/to/cpc2_results.txt"
)

###############################################
## 11. Integrate all analyses
###############################################

switchList <- analyzeSwitchConsequences(
  switchList,
  consequencesToAnalyze = c(
    "intron_retention",
    "coding_potential",
    "NMD_status",
    "ORF_seq_similarity",
    "domain_gain_or_loss",
    "signal_peptide_gain_or_loss"
  )
)

###############################################
## 12. Export results
###############################################

switchSummary <- extractSwitchSummary(switchList)
write.csv(switchSummary, "isoform_switch_summary.csv", row.names = FALSE)

################################