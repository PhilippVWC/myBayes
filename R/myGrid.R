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
