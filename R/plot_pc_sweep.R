#' plot_pc_sweep
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @param rank The rank index of the trees to plot. This must match the trees submitted as input to calc_pc_sweep_var. Eg. If you want to plot trees inferred from PC sweep matrices 3-100, rank is c(3:100).
#' @param var_df The dataframe of variable produced by cellphylo::calc_pc_sweep_var
#'
#' @return List containing plots for tree length, average edge length and star-ness score
#' @export
#'
#' @examples
plot_pc_sweep <- function(rank, var_df){
  #assemble plot dfs
  plot.df <- data.frame(rank = rank, int.branch.len  = unlist(var_df$ave_int_edge), tip.branch.len = unlist(var_df$ave_tip_edge), star.measure = unlist(var_df$star_measure), TL = unlist(var_df$TL),sum.int.edge = unlist(var_df$sum_int_edge), sum.tip.edge = unlist(var_df$sum_tip_edge))

  #plot both on same plot - average branch length
  plot.df.2 <- data.frame(branch_lengths = c(plot.df$int.branch.len,plot.df$tip.branch.len), Edge_Type = c(rep("interior", length(plot.df$int.branch.len)), rep("tip", length(plot.df$tip.branch.len))), x = rep(rank, 2), check.names=FALSE)

  #plot sum of branch lengths instead of average
  plot.df.3 <- data.frame(branch_lengths = c(plot.df$sum.int.edge,plot.df$sum.tip.edge), Edge_Type = c(rep("interior", length(plot.df$int.branch.len)), rep("tip", length(plot.df$tip.branch.len))), x = rep(rank, 2), check.names=FALSE)


  #plots

  #TL
  tree_len.plot <- ggplot(data = plot.df, aes(x=rank, y=TL)) + geom_path() + geom_point() + xlab("PC") + ylab("Tree Length") + theme_bw() + scale_colour_grey(start = 0, end = .9)

  pdf(file="1_tree_length.pdf", width=8, height=5)
  print(tree_len.plot)
  dev.off()

  #ave edgelength
  ave_edge_len.plot <- ggplot(data = plot.df.2, aes(x=x, y=branch_lengths, group=Edge_Type, colour=Edge_Type)) + geom_path() + geom_point(aes(shape=Edge_Type)) + xlab("PC") + ylab("Average Edge Length")  + theme_bw()  + theme(legend.position=c(0.80,0.5),  legend.box.background = element_rect(colour = "black"),legend.key.size = unit(0.25, "cm")) + scale_colour_grey(start = 0, end = .7)

  #print
  pdf(file="2_average_edge_length.pdf", width=8, height=5)
  print(ave_edge_len.plot)
  dev.off()

  #sum of edge lengths
  ave_edge_len.plot <- ggplot(data = plot.df.3, aes(x=x, y=branch_lengths, group=Edge_Type, colour=Edge_Type)) + geom_path() + geom_point(aes(shape=Edge_Type)) + xlab("PC") + ylab("Sum of Edge Lengths")  + theme_bw()  + theme(legend.position=c(0.80,0.5),  legend.box.background = element_rect(colour = "black"),legend.key.size = unit(0.25, "cm")) + scale_colour_grey(start = 0, end = .7)

  pdf(file="2_sum_edge_length.pdf", width=8, height=5)
  print(ave_edge_len.plot)
  dev.off()

  #star measure
  star.plot <- ggplot(data = plot.df, aes(x=rank, y=star.measure)) +  geom_path( ) + geom_point() + xlab("PC") + ylab("Star-ness") + theme_bw() + scale_colour_grey(start = 0, end = .9)

  pdf(file="3_star_measure.pdf", width=8, height=5)
  print(star.plot)
  dev.off()

  return(list(tree_len.plot, ave_edge_len.plot, star.plot))

}
