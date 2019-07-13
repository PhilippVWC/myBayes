#' @title myGrid draws grids
#' @description Extension of R's core function grid(). This function plots a grid only within the box constraints specified
#' by xlim and ylim. The grid will be drawn to the currently open graphics device.
#' @param xlim Vector of type double - The bounds along the x axis for the grid.
#' @param ylim Vector of type double - The bounds along the y axis for the grid.
#' @param nx Integer - Number of vertical lines to be drawn.
#' @param ny Integer - Number of horizontal lines to be drawn. Defaults to nx.
#' @param lty Integer - Specifies the linetype similarly to R's core function grid().
#' @param lwd Integer - Specifies the linewidth similarly to R's core function grid().
#' @param col Character string - The color to be used for the grid lines.
#' @return No return value.
#' @author P.v.W. Crommelin
#' @examples
#' x = seq(from = 0,
#'              to = 3,
#'              by = 0.01)
#' plot(x = x,
#'      y = x,
#'      type = "l")
#' myGrid(xlim = c(1,2),
#'        ylim = c(0.5,2.5),
#'        nx = 10,
#'        ny = 15,
#'        lty = 2,
#'        lwd = 2,
#'        col = "blue")
#' @export
myGrid = function(xlim,ylim,nx,ny=nx,lty = DASHED,lwd = 1,col = "black"){
  dx = (xlim[2]-xlim[1])/(nx-1)
  xticks = seq(from = xlim[1],
               to = xlim[2],
               by = dx)
  dy = (ylim[2]-ylim[1])/(ny-1)
  yticks = seq(from = ylim[1],
               to = ylim[2],
               by = dy)
  sapply(X = xticks,
         FUN = function(tick){
           lines(x = rep(tick,2),
                 y = ylim,
                 lwd = lwd,
                 lty = lty,
                 col = col)
         })
  sapply(X = yticks,
         FUN = function(tick){
           lines(x = xlim,
                 y = rep(tick,2),
                 lwd = lwd,
                 lty = lty,
                 col = col)
         })
}
