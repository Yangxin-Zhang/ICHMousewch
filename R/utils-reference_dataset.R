
# R/utils-reference_dataset.R

#' integrate mouse RNA seq dataset based on cell type
#'

.integrate_mouse_RNA_seq_dataset <- function() {

  on.exit(gc())

  ref_da <- celldex::MouseRNAseqData()

  anno_data <- as.data.table(colData(ref_da))
  anno_data[,assay_label := colData(ref_da)@rownames]

  cell_type <- anno_data[,label.fine] %>%
    unique()

  integrated_expr_matrix <- data.table()
  integrated_expr_matrix[,gene_name := rownames(assay(ref_da))]
  for (i in 1:length(cell_type)) {

    ass_lb <- anno_data[label.fine == cell_type[i],assay_label]

    expr_matrix <- rowMeans(assay(ref_da)[,ass_lb])

    integrated_expr_matrix[,cell_type[i] := expr_matrix]

  }

  return(integrated_expr_matrix)

}
