#' plot_eigs
#'
#' Perform PCA and produce eigenvalue plots for a matrix. These plots are 1) eigenvalue vs PC, 2) eigenvalue cumulative sum, 3) eigenvalue vs percent variance, and 4) percent variance vs PC.
#'
#' @import ggplot2
#' @importFrom Seurat Read10X
#' @importFrom stats prcomp sd
#' @importFrom factoextra get_eigenvalue
#' @importFrom grDevices pdf dev.off
#'
#' @param matrix_path Path to matrix PCA is to be performed on.
#' @param highlight The number of the PC you would like to highlihgt in the plot. The dot representing this PC will be highlighted in red. Default is PC20, used in Mah & Dunn (2023)
#' @param print Boolean. Print out plots.
#' @return plots A list containing the four eigenvalue plots as ggplot objects
#' @export
#'
#' @examples
plot_eigs <- function(matrix_path, highlight=20, print=TRUE){

  #within-species normalized, combined and integrated cross-species. NOT subsetted.
  #unfiltered raw data
  mat <- Read10X(matrix_path)


  #transpose the matrix for the pca function.
  mat <- as.matrix(mat) %>% t()

  #pca
  pca <- prcomp(mat, center=TRUE, scale=FALSE)

  #eigenvalues, % var, cumulative % var
  eigvals <- get_eigenvalue(pca)

  #build plots
  #1. eigenvalue vs PC
  plot.df <- data.frame(order = 1:length(eigvals$eigenvalue), eigenvalue = eigvals$eigenvalue)
  plot.df.highlight <- plot.df[highlight,]
  p1 <- ggplot(plot.df, aes(x=order, y=eigenvalue)) + geom_point() + geom_point(data = plot.df.highlight, colour = "red" ) + theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  #2. eigenvalue cumulative sum
  plot.df <- data.frame(order = 1:length(cumsum(eigvals$eigenvalue)), eigenvalue_cum_sum = cumsum(eigvals$eigenvalue))
  p2 <- ggplot(plot.df, aes(x=order, y=eigenvalue_cum_sum)) + geom_point() +theme_bw() +xlab("order") +ylab("eigenvalue cumulative sum")
  #3 eigenvalue vs percent variance
  plot.df <- data.frame(eigenvalue = eigvals$eigenvalue, percent_variance = eigvals$variance.percent)
  p3 <- ggplot(plot.df, aes(x=eigenvalue, y=percent_variance)) + geom_point() +theme_bw() + xlab("eigenvalue") + ylab("percent variance")

  #how much variance is described by the first few PCs?
  #aqhumor_cross_species_integrated_mtx/: 26.65411
  print(paste0("The amount of variance described by the first ",highlight," PCs is: ",sum(plot.df$percent_variance[1:highlight])))

  # percent variance vs PC
  plot.df <- data.frame(component = 1:length(eigvals$variance.percent), percent_variance = eigvals$variance.percent)
  plot.df.highlight <- plot.df[highlight,]
  p4 <- ggplot(plot.df, aes(x = component, y = percent_variance)) + geom_point(size=0.5)  + xlab("Principal Component") + ylab("Percent Variance") + geom_point(data = plot.df.highlight, colour = "red", size=0.5) + theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(size=6.5), axis.text.y = element_text(size=6.5), axis.title = element_text(size = 7))

  if (print==TRUE){

    if(!dir.exists("eigenvalue_plots")){
      dir.create("eigenvalue_plots")
    }



  #4 variance vs
  #line plot of eigenvalues ordered
  pdf(file="eigenvalue_plots/1_eigenvalue_vs_PC.pdf", width=8, height=5)
  print(p1)
  dev.off()

  #cumulative sum
  pdf(file="eigenvalue_plots/2_eigenvalue_cumulative_sum.pdf", width=8, height=5)
  print(p2)
  dev.off()

  #variance of each component vs eigenvalue
  pdf(file="eigenvalue_plots/3_eigenvalue_vs_variance.pdf", width=5.5, height=5)
  print(p3)
  dev.off()

  #component variance
  pdf(file="eigenvalue_plots/4_eigenvalue_component_variance.pdf",  width=8, height=5)
  print(p4)
  dev.off()
  }

  plots <- list("1_eigenvalue_vs_PC" = p1, "2_eigenvalue_cumulative_sum" = p2, "3_eigenvalue_vs_variance" = p3, "4_eigenvalue_component_variance" = p4)

  return(plots)

}
