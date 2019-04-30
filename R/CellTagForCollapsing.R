#' CellTag Starcode Prior Collapsing
#'
#' This function generate the .txt file that will be fed into starcode - https://github.com/gui11aume/starcode - to collapse similar CellTags.
#' @param ctm.after.whitelist CellTag single-cell matrix after whitelist filtering 
#' @param umi.matrix CellTag raw UMI count matrix
#' @param output.file The filepath and name to save the table for collapsing
#' @return The data frame for collapsing
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagDataForCollapsing(celltags.whitelisted.3, combo.umi.mat.10, "my_favoriate_collapse.csv")
#' 
CellTagDataForCollapsing <- function(celltag.obj, output.file) {
  # Get the data out from the CellTag object
  umi.matrix <- as.matrix(celltag.obj@raw.count)
  for.collapse <- t(umi.matrix)
  # Melt the matrix
  for.collapse <- melt(for.collapse)
  # Subset the matrix to only contain tags with positive UMI numbers
  for.collapse <- subset(for.collapse, value > 0)
  for.collapse$X1 <- as.character(for.collapse$X1)
  for.collapse$X2 <- as.character(for.collapse$X2)
  # Create the contatenation column
  for.collapse$concat <- paste0(for.collapse$X1, unlist(lapply(strsplit(for.collapse$X2, "-"), function(x) x[1])))
  write.table(for.collapse$concat, output.file, sep = "\t", row.names = F, quote = F, col.names = F)
  # Set CellTag object
  celltag.obj@pre.starcode <- for.collapse
  # Print the path saved
  cat("The file for collapsing is stored at: ", output.file, "\n")
  return(celltag.obj)
}

#' CellTag Starcode Post Collapsing
#'
#' This function processes the result generated from starcode - https://github.com/gui11aume/starcode.
#' @param ctm.after.whitelist CellTag single-cell matrix after whitelist filtering 
#' @param collapsed.rslt.file File path to the collapsed result file
#' @param collapsed.csv.file File path to the data frame file generated for collapsing
#' @param output.file The RDS file path and name to save the resulting UMI matrix
#' @return The collapsed and processed UMI matrices
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagDataPostCollapsing(celltags.whitelisted.3, "collapsed_test.txt", "collapsed.csv", "collapsed_data_matrix.Rds")
#' 
CellTagDataPostCollapsing <- function(celltag.obj, collapsed.rslt.file) {
  # Read in the collpased result
  collapsed <- read.table(collapsed.rslt.file, sep = "\t", header = F, stringsAsFactors = F)
  # Read in the file for collapsing
  collapsing <- celltag.obj@pre.starcode
  colnames(collapsing)[c(1:2)] <- c("CellTag", "Cell.Barcode")
  new.collapsing.df <- collapsing
  # Process the collapsing data file
  for (i in 1:nrow(collapsed)) {
    curr.row <- collapsed[i,]
    curr.centroid <- curr.row$V1
    curr.count <- curr.row$V2
    curr.ct <- substring(curr.centroid, 1, 8)
    if (curr.count > 1) {
      curr.collapse.set <- strsplit(curr.row$V3, ",")[[1]]
      curr.to.collapse <- setdiff(curr.collapse.set, curr.centroid)
      for (j in 1:length(curr.to.collapse)) {
        curr.for.c <- curr.to.collapse[j]
        curr.for.c.ct <- substring(curr.for.c, 1, 8)
        if (curr.for.c.ct != curr.ct) {
          ind <- which(collapsing$concat == curr.to.collapse[j])
          ind.cent <- which(collapsing$concat == curr.centroid)
          new.collapsing.df[ind, "concat"] <- curr.centroid
          new.collapsing.df[ind, "CellTag"] <- collapsing[ind.cent[1], "CellTag"]
          new.collapsing.df[ind, "Cell.Barcode"] <- collapsing[ind.cent[1], "Cell.Barcode"]
        }
      }
    }
  }
  new.collapsing.df <- setDT(new.collapsing.df)
  # Regenerate the new matrix
  new.matrix <- dcast(new.collapsing.df, Cell.Barcode~CellTag)
  # Give the matrix rownames
  cell.rnm <- new.matrix$Cell.Barcode
  cnms <- colnames(new.matrix)[2:ncol(new.matrix)]
  new.matrix <- as.matrix(new.matrix[, ..cnms])
  rownames(new.matrix) <- cell.rnm
  # Save the new matrix to the object
  celltag.obj@collapsed.count <- as(new.matrix, "dgCMatrix")
  return(celltag.obj)
}
