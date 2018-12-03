# R Package - CloneHunter
This is a wrapped package of the above workflow with additional checks on the Celltag Library sequences. This package have a dependency on R version (R >= 3.5.1). This can be used as an alternative approach for this pipeline.

(You might need to install devtools to be able to install from Github first)
```r
install.packages("devtools")
```
Install the package from GitHub.
```r
library("devtools")
devtools::install_github("morris-lab/CloneHunter", subdir = "CloneHunter")
```
## Assessment of CellTag Library Complexity via Sequencing
In the first section, we would like to evaluate the CellTag library complexity using sequencing. Following is an example using the sequencing data we generated in lab for pooled CellTag library V2. 
### 1. Read in the fastq sequencing data and extract the CellTags
#### Note: This function is still under construction for extension to accept BAM file as input file. Currently it only works with FASTQ files.
```r
# Read in data file that come with the package
fpath <- system.file("extdata", "V2-1_S2_L001_R1_001.fastq", package = "CloneHunter")
output.path <- "./celltag_extracted_v2-1_r1.txt"
# Extract the CellTags
extracted.cell.tags <- CellTagExtraction(fastq.bam.input = fpath, celltag.version = "v2", extraction.output.filename = output.path, save.fullTag = FALSE, save.onlyTag = FALSE)
```
The extracted CellTag - `extracted.cell.tags` variable - contains a list of two vectors as following.

|First Vector-`extracted.cell.tags[[1]]`|Second Vector-`extracted.cell.tags[[1]]` |
|:-------------------------------------:|:---------------------------------------:|
|Full CellTag with conserved regions    |8nt CellTag region                       |

### 2. Count the CellTags and sort based on occurrence of each CellTag
```r
# Count the occurrence of each CellTag
cell.tag.count <- as.data.frame(table(extracted.cell.tags[[2]]), stringsAsFactors = F)
# Sort the CellTags in descending order of occurrence
cell.tag.count.sort <- cell.tag.count[order(-cell.tag.count$Freq), ]
colnames(cell.tag.count.sort) <- c("CellTag", "Count")
```

### 3. Generation of whitelist for this CellTag library
Here are are generating the whitelist for this CellTag library - CellTag V2. This will remove the CellTags with an occurrence number below the threshold. The threshold (using 90th percentile as an example) is determined via - Cutoff = floor[(90th quantile)/10]. The percentile can be changed while calling the function. Occurrence scatter plots are saved under the `output.dir`, which could be further used to determine the percentile for each different CellTag library.
```r
whitelisted.cell.tag <- CellTagWhitelistFiltering(count.sorted.table = cell.tag.count.sort, percentile = 0.9, output.dir="./", output.file = "my_favourite_v1.csv", save.output = TRUE)
```
The generated whitelist for each library can be used to filter and clean the single-cell CellTag UMI matrices.

## Clone Calling
In this section, we are presenting an alternative approach that utilizes this package that we established to carry out clone calling with single-cell CellTag UMI count matrices. In this pipeline below, we are using a subset dataset generated from the full data (Full data can be found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99915). Briefly, in our lab, we reprogram mouse embryonic fibroblasts (MEFs) to induced endoderm progenitors (iEPs). This dataset is a single-cell dataset that contains cells collected from different time points during the process. This subset is a part of the first replicate of the data. It contains cells collected at Day 28 with three different CellTag libraries - V1, V2 & V3. 

### 1. Read in the single-cell CellTag UMI count matrix
These matrices can be obtained from the first part of this tutorial. In this pipeline below, the matrices were saved as .Rds files. It will need to be changed when saving the matrices into different data types.
```r
# Read the RDS file
dt.mtx.path <- system.file("extdata", "hf1.d28.prefiltered.Rds", package = "CloneHunter")
sc.celltag <- readRDS(dt.mtx.path)
# Change the rownames
rnm <- sc.celltag$Cell.BC
sc.celltag <- sc.celltag[,-1]
rownames(sc.celltag) <- rnm
```

### 2. Binarize the single-cell CellTag UMI count matrix
Here we would like to binarize the count matrix to contain 0 or 1, where 0 indicates no such CellTag found in a single cell and 1 suggests the existence of such CellTag. The suggested cutoff that marks existence or absence is at least 2 counts per CellTag per Cell. For details, please refer to the paper - <<<<<<<<<<<<<<<<<<<<<<<<<<<LINK>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
```r
# Calling binarization
binary.sc.celltag <- SingleCellDataBinarization(sc.celltag, 2)
```

### 3. Metric plots to facilitate for additional filtering
We then generate scatter plots for the number of total celltag counts in each cell and the number each tag across all cells. These plots could help us further in filtering and cleaning the data.
```r
metric.p <- MetricPlots(celltag.data = binary.sc.celltag)
print(paste0("Mean CellTags Per Cell: ", metric.p[[1]]))
print(paste0("Mean CellTag Frequency across Cells: ", metric.p[[2]]))
```

### 4. Apply the whitelisted CellTags generated from assessment
##### Note: This filters the single-cell data based on the whitelist of CellTags one by one. By mean of that, if three CellTag libraries were used, the following commands need to be executed for 3 times and result matrices can be further joined (Example provided).

Based on the whitelist generated earlier, we filter the UMI count matrix to contain only whitelisted CelTags.
```r
whitelist.sc.data.v2 <- SingleCellDataWhitelist(binary.sc.celltag, whitels.cell.tag.file = "./my_favourite_v2_1.csv")
##########
# Only run if this sample has been tagged with more than 1 CellTags libraries
##########
## NOT RUN
whitelist.sc.data.v1 <- SingleCellDataWhitelist(binary.sc.celltag, whitels.cell.tag.file = "./my_favourite_v1.csv")
whitelist.sc.data.v3 <- SingleCellDataWhitelist(binary.sc.celltag, whitels.cell.tag.file = "./my_favourite_v3.csv")
# Concatenate the colnames with celltag versions
colnames(whitelist.sc.data.v1) <- paste("V1_", colnames(whitelist.sc.data.v1), sep = "")
colnames(whitelist.sc.data.v3) <- paste("V3_", colnames(whitelist.sc.data.v3), sep = "")
```
