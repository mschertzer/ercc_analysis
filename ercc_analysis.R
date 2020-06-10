#!/usr/bin/env Rscript

# Megan D Schertzer
# June 3, 2020
# Absolute RNA quantification using ERCC spike-ins

# usage: ./ercc_analysis.r -d 100 -e 2 -p 30 -r 1 <ERCC_FEATURECOUNTS_FILE> <GENE_FEATURECOUNTS_FILE> <BAM_READ_COUNT>
#        ./ercc_analysis.r --dilution 100 --ercc_microliter 2 --rna_per_cell_picogram 30 --rna_microgram 1 <ERCC_FEATURECOUNTS_FILE> <GENE_FEATURECOUNTS_FILE> <BAM_READ_COUNT>

# load packages
cat("Installing and loading required packages... \n")
cat("\n")
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(readr)) install.packages('readr')
library(readr)

if (!require(optparse)) install.packages('optparse')
library(optparse)

# command line arguements
cat("\n")
cat("Parsing command line arguments... \n")
cat("\n")
arguments <- parse_args(OptionParser(usage = "%prog [options] <ERCC_FEATURECOUNTS_FILE> <GENE_FEATURECOUNTS_FILE> <BAM_READ_COUNT>",
                                     option_list=list(
                                       make_option(c("-d","--dilution"), default = 100, type = "double", help = "The dilution of ercc spike-ins prepared prior to library prep. [Default %default] is same as used in paper"),
                                       make_option(c("-e","--ercc_microliter"), default = 2, type = "double", help = "Microliter amount of ercc spike-in dilution added to rna before library prep. [Default %default] is same as used in paper."),
                                       make_option(c("-p", "--rna_per_cell_picogram"), default = 30, type = "double", help = "This is the amount of total RNA per cell for your cell type. This is cell type specific and must be calculated for your cell type as described in the paper."),
                                       make_option(c("-r","--rna_microgram"), default = 1, type = "double", help = "Micrograms of starting total cellular RNA used in library prep. [Default %default] is same as used in paper."))), 
                        positional_arguments = 3)

# assign variables
opt = arguments$opt
ercc_name = arguments$args[1]
gene_name = arguments$args[2]
read_count =  as.numeric(arguments$args[3])
out_name = strsplit(arguments$args[2], "\\.")[[1]][1]

# Make sure the information is correct
cat("Here is the information that this script will use to analyze your data. Check that it is correct:\n")
cat("1. The featureCounts for ERCC spike-ins is in a file called", ercc_name, "\n")
cat("2. The featureCounts for all genes of interest is in a file called", gene_name, "\n")
cat("3. The total read counts from your sam/bam alignment file is", read_count, "\n")
cat("4. The dilution of ERCC spike-ins used in your experiment was 1 /", opt$dilution, "\n")
cat("5. The microliter amount of ERCC spike-ins used in your experiment was", opt$ercc_microliter, "\n")
cat("6. The picogram amount of total RNA/cell for your cell type is", opt$rna_per_cell_picogram, "\n")
cat("7. The microgram amount of total cellular RNA that you started the library with was", opt$rna_microgram, "\n")
cat("\n")

# Read in ercc file and go through checkpoints
cat("Starting analysis of", ercc_name,"... \n")
cat("\n")
if (!file.exists(ercc_name)) stop("File ",ercc_name," does not exist. Check that the featureCounts file is in the proper directory. \n")
ercc_fc <- read_tsv(ercc_name, col_names = TRUE, skip = 1)
if (nrow(ercc_fc) > 93) stop("ERCC featureCounts file has more rows than expected. Check the file order that you entered on the command line. The ERCC file should be first. \n")
ercc_fc <- ercc_fc %>% 
  select(1,6,7) %>% 
  rename(ERCCID=Geneid, Counts=3)

# Combine this ERCC concentration file (https://www.thermofisher.com/order/catalog/product/4456740#/4456740) with the ERCC featureCounts file
# ERCC concentrations from Mix 1
if (!file.exists("ERCC92_conc.txt")) stop("File ERCC92_conc.txt does not exist. Check that the file is in the proper directory. You can download this file from https://www.thermofisher.com/order/catalog/product/4456740#/4456740. \n")
ercc_conc <- read_tsv(file="ERCC92_conc.txt", col_names = TRUE)
ercc_conc <- ercc_conc %>% 
  select(1:4) %>% 
  rename(ERCCID = `ERCC ID`, Conc.Mix1 = 'concentration in Mix 1 (attomoles/ul)')

# Check the structure of the two dataframes- note the common column is called 'ERCCID'
#glimpse(ercc_conc)
#glimpse(ercc_fc)

# Combine this ERCC concentration file with the ERCC featureCounts file based on the ERCCID column
ercc_all <- ercc_conc %>% 
  full_join(., ercc_fc, by = "ERCCID")

# Perform calculations to compare ERCC RPKM versus ERCC copies
# transform Mix 1 concentrations= values to reflect dilutions in experiments
# 1 x 10^18 attomoles = 1 mole
# avogadros number = 6.022 x 10^23 molecules/mole
ercc_calculations <- ercc_all %>% 
  mutate(Attomoles.added = (Conc.Mix1/opt$dilution)*opt$ercc_microliter) %>% 
  mutate(Moles.added = Attomoles.added/(10^18)) %>% 
  mutate(Molecules.added = Moles.added * (6.022 * 10^23)) %>%
  mutate(RPM = Counts/(read_count/1000000)) %>% 
  mutate(RPKM = RPM/(Length/1000)) %>%
  mutate(LogMolecules = log2(Molecules.added)) %>%
  mutate(LogRPKM = log2(RPKM))

# Check out the structure of the new table
#glimpse(ercc_calculations)

# filter out ERCC transcripts that had zero reads align
ercc_calculations <- ercc_calculations %>% 
  arrange(., desc(LogRPKM)) %>% 
  filter(., LogRPKM != "-Inf")

# Plot ERCC standard curve, fit a regression line
plot_name = paste(out_name, "ercc_plot.pdf", sep="_")

pdf(file = plot_name, width=6, height=5)
plot(ercc_calculations$LogMolecules, ercc_calculations$LogRPKM, pch=16, xlab="log2(Molecules)", ylab="log2(Rpkm)", main=paste("ERCC Standards for", out_name, sep=" "))
quant <- lm(LogRPKM ~ LogMolecules, data=ercc_calculations)
abline(quant)
cf <- round(coef(quant), 3)
eq <- paste0("y = ",cf[1], ifelse(sign(cf[2])==1, "+", "-"), abs(cf[2]), "x")
mtext(eq, 3, line= -2)
r <- paste0("R2 = ", format(summary(quant)$r.squared, digits=3))
mtext(r, 3, line=-3)
dev.off()

# Load and process gene featureCount files
cat("\n")
cat("Starting analysis of", gene_name,"... \n")
cat("\n")

if (!file.exists(gene_name)) stop("File ",gene_name," does not exist. Check that the featureCounts file is in the proper directory. \n")
fc <- read_tsv(gene_name, col_names = TRUE, skip = 1)
genes_fc <- fc %>% rename(Counts = 7)
cat("The featureCounts file has", nrow(genes_fc), "genes. \n")
cat("\n")

# Use y=mx+b equation from ERCC Standards fitted line to calculate molecules per cell
m <- coef(quant)[2]
b <- coef(quant)[1]
cat("The slope of the ERCC standard curve is", m, "\n")
cat("The y-intercept of the ERCC standard curve is", b, "\n")
cat("\n")

# Calculate molecules per cell for each RNA transcript
# The last step in this calculation is from molecules/Xug RNA to molecules/cell which requires you to know the amount of RNA/per cell in your specific cell type and the starting amount of total RNA for the library prep (the latter is set as a parameter in the second chunk of this notebook)
# In TSCs, we measured 30 picograms/cell
mpc <- genes_fc %>% 
  mutate(RPM = Counts/(read_count/1000000)) %>% 
  mutate(RPKM = RPM/(Length/1000)) %>%
  mutate(LogRPKM = log2(RPKM)) %>%
  mutate(LogMolecules = (LogRPKM +abs(b))/m) %>%
  mutate(Molecules = 2^LogMolecules) %>%
  mutate(MoleculesPer.ug = Molecules*(opt$rna_per_cell_picogram/1000000)) %>%
  mutate(MPC = MoleculesPer.ug/opt$rna_microgram)

# Check out the structure of the new table
#glimpse(mpc) 

# Write a table as output
write_tsv(mpc, path = paste(out_name, "mpc.txt", sep="_"))

cat("Analysis is finished. There should be one output plot and one output table. \n")
cat("\n")
