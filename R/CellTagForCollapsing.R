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
CellTagDataForCollapsing <- function(ctm.after.whitelist, umi.matrix, output.file) {
  # Generate filenames for collapse txt
  dirnm <- dirname(output.file)
  bnm <- basename(output.file)
  fnm.parts <- strsplit(bnm, "[.]")[[1]]
  fnm <- fnm.parts[c(1:(length(fnm.parts) - 1))]
  fn.txt <- paste0(dirnm, "/", paste(c(fnm, "collapse.txt"), collapse = "_"))
  # Get the filtered CellTags and Cell Barcodes
  filtered.names <- rownames(ctm.after.whitelist)
  for.collapse <- as.matrix(umi.matrix[, filtered.names])
  # Melt the matrix
  for.collapse <- melt(for.collapse)
  # Subset the matrix to only contain tags with positive UMI numbers
  for.collapse <- subset(for.collapse, value > 0)
  # Create the contatenation column
  for.collapse$concat <- paste0(for.collapse$X1, unlist(lapply(strsplit(for.collapse$X2, "-"), function(x) x[1])))
  write.csv(for.collapse, output.file, row.names = F, quote = F)
  print(fn.txt)
  write.table(for.collapse$concat, fn.txt, sep = "\t", row.names = F, quote = F, col.names = F)
  return(for.collapse)
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
CellTagDataPostCollapsing <- function(ctm.after.whitelist, collapsed.rslt.file, collapsed.csv.file, output.file) {
  # Read in the collpased result
  collapsed <- read.table(collapsed.rslt.file, sep = "\t", header = F, stringsAsFactors = F)
  # Read in the file for collapsing
  collapsing <- read.csv(collapsed.csv.file, stringsAsFactors = F, header = T)
  colnames(collapsing)[c(1:2)] <- c("CellTag", "Cell.Barcode")
  new.collapsing.df <- collapsing
  # Process the collapsing data file
  for (i in 1:nrow(collapsed)) {
    curr.row <- collapsed[i,]
    curr.centroid <- curr.row$V1
    curr.count <- curr.row$V2
    if (curr.count > 1) {
      curr.collapse.set <- strsplit(curr.row$V3, ",")[[1]]
      curr.to.collapse <- setdiff(curr.collapse.set, curr.centroid)
      for (j in 1:length(curr.to.collapse)) {
        ind <- which(collapsing$concat == curr.to.collapse[j])
        ind.cent <- which(collapsing$concat == curr.centroid)
        new.collapsing.df[ind, "concat"] <- curr.centroid
        new.collapsing.df[ind, "CellTag"] <- collapsing[ind.cent[1], "CellTag"]
        new.collapsing.df[ind, "Cell.Barcode"] <- collapsing[ind.cent[1], "Cell.Barcode"]
      }
    }
  }
  # Regenerate the new matrix
  new.matrix <- cast(new.collapsing.df, CellTag~Cell.Barcode)
  # Give the matrix rownames
  cell.tag.rnm <- new.matrix$CellTag
  new.matrix <- new.matrix[,-1]
  rownames(new.matrix) <- cell.tag.rnm
  # Save the new matrix
  saveRDS(new.matrix, output.file)
  return(new.matrix)
}
