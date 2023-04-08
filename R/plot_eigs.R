#' plot_eigs
#'
#' Perform PCA and produce eigenvalue plots for a matrix.
#'
#' @import ggplot2
#' @importFrom Seurat Read10X
#' @importFrom stats prcomp
#' @importFrom factoextra get_eigenvalue
#'
#' @param matrix_path Path to matrix PCA is to be performed on.
#'
#' @return NULL
#' @export
#'
#' @examples
plot_eigs <- function(matrix_path){

  #within-species normalized, combined and integrated cross-species. NOT subsetted.
  #unfiltered raw data
  mat <- Read10X(matrix_path)


  #transpose the matrix for the pca function.
  mat <- as.matrix(mat) %>% t()

  #pca
  pca <- prcomp(mat, center=TRUE, scale=FALSE)

  #eigenvalues, % var, cumulative % var
  eigvals <- get_eigenvalue(pca)

  #line plot of eigenvalues ordered
  pdf(file="1_eigenvalue_vs_PC.pdf", width=8, height=5)
  plot.df <- data.frame(order = 1:length(eigvals$eigenvalue), eigenvalue = eigvals$eigenvalue)
  plot.df.highlight <- plot.df[20,]
  p1 <- ggplot(plot.df, aes(x=order, y=eigenvalue)) + geom_point() + geom_point(data = plot.df.highlight, colour = "red" ) + theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(p1)
  dev.off()

  #cumulative sum
  pdf(file="2_eigenvalue_cumulative_sum.pdf", width=8, height=5)
  plot.df <- data.frame(order = 1:length(cumsum(eigvals$eigenvalue)), eigenvalue_cum_sum = cumsum(eigvals$eigenvalue))
  p2 <- ggplot(plot.df, aes(x=order, y=eigenvalue_cum_sum)) + geom_point() +theme_bw() +xlab("order") +ylab("eigenvalue cumulative sum")
  print(p2)
  dev.off()

  #variance of each component vs eigenvalue
  pdf(file="3_eigenvalue_vs_variance.pdf", width=5.5, height=5)
  plot.df <- data.frame(eigenvalue = eigvals$eigenvalue, percent_variance = eigvals$variance.percent)
  p3 <- ggplot(plot.df, aes(x=eigenvalue, y=percent_variance)) + geom_point() +theme_bw() + xlab("eigenvalue") + ylab("percent variance")
  print(p3)
  dev.off()

  #how much variance is described by the first few PCs?
  #aqhumor_cross_species_integrated_mtx/: 26.65411
  sum(plot.df$percent_variance[1:20])

  #component variance
  pdf(file="4_eigenvalue_component_variance_red_dot_7feb2023.pdf", width=3.4, height=2.2)
  plot.df <- data.frame(component = 1:length(eigvals$variance.percent), percent_variance = eigvals$variance.percent)
  plot.df.highlight <- plot.df[20,]
  p4 <- ggplot(plot.df, aes(x = component, y = percent_variance)) + geom_point(size=0.5)  + xlab("Principal Component") + ylab("Percent Variance") + geom_point(data = plot.df.highlight, colour = "red", size=0.5) + theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(size=6.5), axis.text.y = element_text(size=6.5), axis.title = element_text(size = 7))
  print(p4)
  dev.off()

}
