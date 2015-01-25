#' @title readExcelSheet
#' @description Function reads data from an Excel file from a specified sheet
#' @author Lech Madeyski
#' @export readExcelSheet
#' @param path Path to an Excel file, e.g. /User/lma/datasets/MyDataSet.xls
#' @param sheet Name of a sheet within an Excel file we want to read
#' @param colNames If TRUE, first row of data will be used as column names.
#' @examples
#' myPath=system.file("extdata", "DataSet.xlsx", package = "reproducer")
#' Madeyski15SQJ.NDC<-readExcelSheet(path=myPath, sheet="Madeyski15SQJ.NDC", colNames=FALSE)
readExcelSheet <- function(path, sheet, colNames){
  dataset = openxlsx::read.xlsx(xlsxFile=path, sheet=sheet, colNames=colNames)
  return(dataset)
}

#' @title densityCurveOnHistogram
#' @description Density curve overlaid on histogram
#' @author Lech Madeyski
#' @export densityCurveOnHistogram
#' @param df Data frame with data to be displayed
#' @param colName Name of the selected column in a given data frame
#' @param limLow the limit on the lower side of the displayed range
#' @param limHigh the limit on the higher side of the displayed range
#' @return A figure being a density curve overlaid on histogram
#' @examples
#' densityCurveOnHistogram(Madeyski15EISEJ.PropProjects, "STUD", 0, 100)
#' densityCurveOnHistogram(data.frame(x<-rnorm(50, mean=50, sd=5)), "x", 0, 100)
densityCurveOnHistogram <- function(df, colName, limLow, limHigh) {
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=colName), environment = environment()) +
    ggplot2::geom_histogram(ggplot2::aes_string(y="..density.."), environment = environment(), fill="cornsilk", colour="grey60", size=.2) +
    ggplot2::geom_density(fill = 'green', alpha = 0.4) +
    ggplot2::xlim(limLow, limHigh)
  return(p1)
}

#' @title boxplotHV
#' @description Box plot
#' @author Lech Madeyski
#' @export boxplotHV
#' @param df Data frame with data to be displayed
#' @param colName Name of the selected column in a given data frame
#' @param limLow the limit on the lower side of the displayed range
#' @param limHigh the limit on the higher side of the displayed range
#' @param isHorizontal Boolean value to control whether the box plot should be horizontal or not (i.e., vertical)
#' @return A box plot
#' @examples
#' boxplotHV(Madeyski15EISEJ.PropProjects, "STUD", 0, 100, TRUE)
boxplotHV <- function(df, colName, limLow, limHigh, isHorizontal) {
  p2 <- ggplot2::ggplot(df, ggplot2::aes_string(x=1, y = colName), environment = environment()) +
    ggplot2::geom_boxplot(fill = 'orange') + ggplot2::theme_bw() + ggplot2::ylim(limLow, limHigh)  +
    ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, size=4, fill="white") +
    ggplot2::scale_x_continuous(breaks=NULL)

  if(isHorizontal) {
    p2 <- p2 + ggplot2::coord_flip()
  }
  return(p2)
}

#' @title boxplotAndDensityCurveOnHistogram
#' @description Boxplot and density curve overlaid on histogram
#' @author Lech Madeyski
#' @export boxplotAndDensityCurveOnHistogram
#' @param df Data frame with data to be displayed
#' @param colName Name of the selected column in a given data frame
#' @param limLow the limit on the lower side of the displayed range
#' @param limHigh the limit on the higher side of the displayed range
#' @return A figure being a density curve overlaid on histogram
#' @examples
#' library(grid)
#' library(gridExtra)
#' library(ggplot2)
#' boxplotAndDensityCurveOnHistogram(Madeyski15EISEJ.PropProjects, "STUD", 0, 100)
#' boxplotAndDensityCurveOnHistogram(Madeyski15SQJ.NDC, "simple", 0, 100)
boxplotAndDensityCurveOnHistogram <- function(df, colName, limLow, limHigh) {
  p1 <- densityCurveOnHistogram(df, colName, limLow, limHigh)
  p2 <- boxplotHV(df, colName, limLow, limHigh, TRUE)

  #arrange the plots together, with appropriate height and width for each row and column
  p1 <- p1 + ggplot2::ylab("") + ggplot2::labs(
    axis.title.x = ggplot2::element_blank(),
    text = ggplot2::element_text(),
    x = "")
  p1 <- p1  + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  p2 <- p2 + ggplot2::xlab("")
  p2 <- p2 + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

  plot_empty <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1,1), colour="white") +
  ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )
  plot_y_labels <- p1 +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(),
      axis.text.x = ggplot2::element_text(colour="white")
      ) +
    ggplot2::xlab("")
  gridExtra::grid.arrange(plot_y_labels, p1, plot_empty, p2, ncol=2, nrow=2, widths=c(1, 19), heights=c(5, 1))
}

