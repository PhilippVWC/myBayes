#' @title Optimize in multiple dimensions
#' @description This optimization routine performs minimum search in more than two dimensions. Variuous optimization algorithms
#' such as (generalized) simmulated annealing, BFGS (Broyden-Fletcher-Goldfarb-Shanno), Rvmmin, (hjkb) Hook-Jeeves and Nelder-Mead are subsequently performed.
#' @param candidates data frame - With one set of parameters (par) for each row.
#' @param  fn function of par - The function that is to be minimized.
#' @param lower Vector of type double - Lower bounds. The variable lower has the same dimension as par (Number of columns of candidates).
#' @param upper Vector of type double - Upper bounds. The variable upper has the same dimension as par (Number of columns of candidates).
#' @param gr A function to return the gradient for BFGS. If set to NULL, finite difference approximations
#' of the gradient are used instead (see documentation of package "optim")
#' @param maxit Maximum amount of iterations for BFGS, RVmmin and general simmulated Annealing.
#' @param statusMessage Boolean - If set to true (default) the output is more verbose.
#' @details This function performs adequate error handling. In problematic cases, where optimization fails, the program is not interrupted but may
#' return an empty data frame. Many optimization algorithms perform their tasks, but only those, which do successfully contribute to a column entry in
#' the resulting data frame.
#' @return data frame - One column for every successfully performed optimization algorithm. Every row corresponds to an optimized parameter. Last row corresponds
#' to the function value of the resulting optimized parameters.
#' @author Philipp van Wickevoort Crommelin
#' @examples
#' fn = function(x) x[1]^2+x[2]^2
#' candidates = data.frame("x" = c(1,-30),
#'                         "y" = c(-5,10))
#' lower = c(-100,-100)
#' upper = c(100,100)
#'
#' res = myOpt_N(candidates = candidates,
#'               fn = fn,
#'               lower = lower,
#'               upper = upper,
#'               maxit = 200,
#'               statusMessages = TRUE)
#' print(res)# Columns are sorted in increasing order
#' @export
myOpt_N = function(candidates,fn,lower,upper,gr=NULL,maxit = 200,statusMessages = TRUE){
  MX = apply(X = as.matrix(candidates),
             MARGIN = 1,
             FUN = function(par){
               dim = length(par)

               newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.

               # Generalized simmulated Annealing
               tryCatch(expr = {
                 if (statusMessages) print("GenSA")
                 res = GenSA::GenSA(par = par,
                                    fn = fn,
                                    lower = lower,
                                    upper = upper,
                                    control = list("maxit" = maxit,
                                                   "max.time" = 10))
                 newVals = c(newVals,res$par,res$value,3)
               },
               warning = function(w){
                 if (statusMessages) print(paste0("GenSA::GenSA WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("GenSA::GenSA ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Simulated Annealing - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("stats::optim SANN")
                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "SANN",
                                    control = list(
                                      temp = 1, #Starting temperature
                                      tmax = 5 #number of function evaluations per temperature step
                                    ))
                 if (!FALSE %in% (res$par <= lower & res$par >= upper)){
                   newVals = c(newVals,res$par,res$value,4)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim Simulated annealing WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim Simulated annealing ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Quasi Newton algorithm with approximated Hessian matrix - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("stats::optim BFGS")
                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "BFGS",
                                    gr = NULL, # if no analytical gradient function is provided, a finite difference method is used
                                    control = list(
                                      maxit = maxit,
                                      reltol = 1e-8
                                    ))
                 if (!FALSE %in% (res$par <= upper & res$par >= lower)){
                   newVals = c(newVals,res$par,res$value,5)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim BFGS WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim BFGS ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Quasi Newton algorithm with approximated Hessian matrix - BOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("Rvmmin")
                 res = Rvmmin::Rvmmin(par = par,
                                      fn = fn,
                                      gr = "grfwd", # "grfwd" = finite difference forward gradient (numerical gradient)
                                      lower = lower,
                                      upper = upper,
                                      maxit = maxit) #Stop computation, if estimators are out of bounds
                 newVals = c(newVals,res$par,res$value,6)
               },
               warning = function(w){
                 if (statusMessages) print(paste0("RVMMIN WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){

                 if (statusMessages) print(paste0("RVMMIN ERROR: ", e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})


               # Hooke-Jeeves algorithm (Derivative free) - BOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("dfoptim::hjkb")
                 res = dfoptim::hjkb(par = par,
                                     fn = fn,
                                     lower = lower,
                                     upper = upper,
                                     control = list())
                 newVals = c(newVals,res$par,res$value,7)
               },
               warning = function(w){
                 if (statusMessages) print(paste0("DFOPTIM::HJKB WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("DFOPTIM::HJKB ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Nelder Mead (Derivative free) - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("stats::optim Nelder-Mead")
                 res = stats::optim(par = par,
                                    method = "Nelder-Mead",
                                    fn = fn)
                 newVals = c(newVals,res$par,res$value,8)
                 if (!FALSE %in% (res$par <= upper & res$par >= lower)){
                   newVals = c(newVals,res$value,res$par,8)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim Nelder-Mead WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim Nelder-Mead ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})


               # Lightweight BFGS - BFGS with box constraints
               tryCatch(expr= {
                 if (statusMessages) print("stats::optim L-BFGS-B")
                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "L-BFGS-B",
                                    lower = lower,
                                    upper = upper,
                                    control = list(maxit = maxit)
                 )
                 newVals = unname(c(newVals,res$par,res$value,9))
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim L-BFGS-B WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim L-BFGS-B ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})
               return(newVals)
             })
  names = c(names(candidates),"value","methodCode")
  methods = c("invalid","optimize","GenSA","SA","BFGS","Rvmmin","Hooke-Jeeves","Nelder-Mead","L-BFGS-B")
  Nrows = length(names)
  MX_df = data.frame(matrix(data = MX,
                            nrow = Nrows))
  rownames(x = MX_df) = names
  colnames(x = MX_df) = methods[unlist(MX_df[Nrows,])] #Take last row with method codes to name the columns.
  MX_df = MX_df[-Nrows,] #Remove method code row.
  MX_df = MX_df[,colnames(MX_df) != "invalid"] #Remove invalid columns.
  MX_df = MX_df[,order(MX_df[nrow(MX_df),])]
  return(MX_df)
}

