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
#' @param statusMessage Boolean - Set to TRUE for more verbose output (default to FALSE).
#' @details This function performs adequate error handling. In problematic cases, where optimization fails, the program is not interrupted but may
#' return an empty data frame. Many optimization algorithms perform their tasks, but only those, which do successfully contribute to a column entry in
#' the resulting data frame.
#' @return data frame - One column for every successfully performed optimization algorithm. Every row corresponds to an optimized parameter. Last row corresponds
#' to the function value of the resulting optimized parameters.
#' @author Philipp van Wickevoort Crommelin
#' @examples
#' install.packages(c("GenSA","dfoptim")
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
myOpt_N = function(candidates,fn,lower,upper,gr=NULL,maxit = 200,statusMessages = FALSE){

  #Built matrix. Every entry of candidates correspond to one single column in matrix.
  MX = apply(X = as.matrix(candidates),
             MARGIN = 1,
             FUN = function(par){
               dim = length(par)
               locEnv = environment()
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
                 newVals = c(newVals,res$par,res$value,2)
               },
               warning = function(w){
                 if (statusMessages) print(paste0("GenSA::GenSA WARNING: ",w))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("GenSA::GenSA ERROR: ",e))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
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
                 if (!(FALSE %in% (res$par <= upper & res$par >= lower))){
                   newVals = c(newVals,res$par,res$value,3)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim Simulated annealing WARNING: ",w))
                 return()
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim Simulated annealing ERROR: ",e))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               finally = {})

               # Quasi Newton algorithm with approximated Hessian matrix - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("stats::optim BFGS")
                 print("####################")
                 print("BFGS expr block")
                 print(environment())
                 print(parent.env(environment()))
                 print("####################")

                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "BFGS",
                                    gr = NULL, # if no analytical gradient function is provided, a finite difference method is used
                                    control = list(
                                      maxit = maxit,
                                      reltol = 1e-8
                                    ))
                 if (!(FALSE %in% (res$par <= upper & res$par >= lower))){
                   newVals = c(newVals,res$par,res$value,4)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim BFGS WARNING: ",w))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim BFGS ERROR: ",e))
                 print("####################")
                 print("BFGS error block")
                 print(environment())
                 print(parent.env(environment()))
                 print("####################")
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
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
                 newVals = c(newVals,res$par,res$value,5)
               },
               warning = function(w){
                 if (statusMessages) print(paste0("DFOPTIM::HJKB WARNING: ",w))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("DFOPTIM::HJKB ERROR: ",e))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               finally = {})

               # Nelder Mead (Derivative free) - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessages) print("stats::optim Nelder-Mead")
                 res = stats::optim(par = par,
                                    method = "Nelder-Mead",
                                    fn = fn)
                 if (!(FALSE %in% (res$par <= upper & res$par >= lower))){
                   newVals = c(newVals,res$value,res$par,6)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim Nelder-Mead WARNING: ",w))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim Nelder-Mead ERROR: ",e))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
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
                                    control = list()
                 )
                 newVals = unname(c(newVals,res$par,res$value,7))
               },
               warning = function(w){
                 if (statusMessages) print(paste0("stats::optim L-BFGS-B WARNING: ",w))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               error = function(e){
                 if (statusMessages) print(paste0("stats::optim L-BFGS-B ERROR: ",e))
                 #proper environment
                 assign(x = "newVals",
                        value = c(newVals,0,rep(0,dim),1),#Method code 1 is a (non informative) gap filler - Will be later erased.
                        envir = locEnv)
                 return()
               },
               finally = {})
               return(newVals)
             })
  print(MX)
  nms = c(names(candidates),"value","methodCode")#names for rows of data frame (see below)
  methods = c("invalid","GenSA","SA","BFGS","Hooke-Jeeves","Nelder-Mead","L-BFGS-B")
  Nrows = length(nms)
  #Feed each set of optimizer output (set of 2+#par values) into dataframe
  MX_df = data.frame(matrix(data = MX,
                            nrow = Nrows))
  rownames(MX_df) = nms
  colnames(MX_df) = methods[unlist(MX_df[Nrows,])] #Take last row with method codes to name the columns.
  MX_df = MX_df[,MX_df[Nrows,]!=1] #remove invalid column entries
  MX_df = MX_df[-Nrows,] #Remove method code row.
  MX_df = MX_df[,order(MX_df[nrow(MX_df),])] #sort according to optimization results. Lowest values is left outermost entry.
  return(MX_df)
}

