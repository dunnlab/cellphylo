#' plot_pc_sweep
#'
#' Plots the variables calculated by calc_pc_sweep_var. Produces Figs 2A, 2B, 2C. Default values are what were used in Mah & Dunn (2023).
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @param rank The range of the ranks you want to plot. This sets the x-axis for hte plots.The range must match the number of rows in the data frame of variables calculted by `calc_pc_sweep_var`.
#' @param var_df The dataframe of variable produced by cellphylo::calc_pc_sweep_var
#'
#' @return List containing plots for tree length, average edge length and star-ness score
#' @export
#'
#' @examples
plot_pc_sweep <- function(rank=c(3:100), var_df, print=TRUE){
  #assemble plot dfs
  plot.df <- data.frame(rank = rank, int.branch.len  = unlist(var_df$ave_int_edge), tip.branch.len = unlist(var_df$ave_tip_edge), star.measure = unlist(var_df$star_measure), TL = unlist(var_df$TL),sum.int.edge = unlist(var_df$sum_int_edge), sum.tip.edge = unlist(var_df$sum_tip_edge))

  #plot both on same plot - average branch length
  plot.df.2 <- data.frame(branch_lengths = c(plot.df$int.branch.len,plot.df$tip.branch.len), Edge_Type = c(rep("interior", length(plot.df$int.branch.len)), rep("tip", length(plot.df$tip.branch.len))), x = rep(rank, 2), check.names=FALSE)

  #plot sum of branch lengths instead of average
  plot.df.3 <- data.frame(branch_lengths = c(plot.df$sum.int.edge,plot.df$sum.tip.edge), Edge_Type = c(rep("interior", length(plot.df$int.branch.len)), rep("tip", length(plot.df$tip.branch.len))), x = rep(rank, 2), check.names=FALSE)


  #plots

  #TL
  tree_len.plot <- ggplot(data = plot.df, aes(x=rank, y=TL)) + geom_path() + geom_point() + xlab("PC") + ylab("Tree Length") + theme_bw() + scale_colour_grey(start = 0, end = .9)


  #ave edgelength
  #ave_edge_len.plot <- ggplot(data = plot.df.2, aes(x=x, y=branch_lengths, group=Edge_Type, colour=Edge_Type)) + geom_path() + geom_point(aes(shape=Edge_Type)) + xlab("PC") + ylab("Average Edge Length")  + theme_bw()  + theme(legend.position=c(0.80,0.5),  legend.box.background = element_rect(colour = "black"),legend.key.size = unit(0.25, "cm")) + scale_colour_grey(start = 0, end = .7)

  #sum of edge lengths
  edge_len.plot <- ggplot(data = plot.df.3, aes(x=x, y=branch_lengths, group=Edge_Type, colour=Edge_Type)) + geom_path() + geom_point(aes(shape=Edge_Type)) + xlab("PC") + ylab("Sum of Edge Lengths")  + theme_bw()  + theme(legend.position=c(0.80,0.5),  legend.box.background = element_rect(colour = "black"),legend.key.size = unit(0.25, "cm")) + scale_colour_grey(start = 0, end = .7)



  #star measure
  star.plot <- ggplot(data = plot.df, aes(x=rank, y=star.measure)) +  geom_path( ) + geom_point() + xlab("PC") + ylab("Star-ness") + theme_bw() + scale_colour_grey(start = 0, end = .9)

  if (print==TRUE){

    #create a PC sweep directory if doesn't already exist
    if(!dir.exists("PC_sweep")){
      dir.create("PC_sweep")
    }

    #create a directory for the combined matrix
    if(!dir.exists("PC_sweep/plots")){
      dir.create("PC_sweep/plots")
    }

    #print
    #TL
    pdf(file="1_tree_length.pdf", width=8, height=5)
    print(tree_len.plot)
    dev.off()

    #ave edge length
    #pdf(file="2_average_edge_length.pdf", width=8, height=5)
    #print(ave_edge_len.plot)
    #dev.off()

    #sum of edge lengths
    pdf(file="2_sum_edge_length.pdf", width=8, height=5)
    print(edge_len.plot)
    dev.off()

    #starness score
    pdf(file="3_star_measure.pdf", width=8, height=5)
    print(star.plot)
    dev.off()

    file.rename("1_tree_length.pdf", "PC_sweep/plots/1_tree_length.pdf")
    file.rename("2_sum_edge_length.pdf", "PC_sweep/plots/2_sum_edge_length.pdf")
    file.rename("3_star_measure.pdf", "PC_sweep/plots/3_star_measure.pdf")


  }

  return(list(tree_len.plot, edge_len.plot, star.plot))

}
