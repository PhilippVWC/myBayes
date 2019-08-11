#' @export
myOpt_N = function(candidates,fn,lower,upper,control=list(),gr=NULL,maxit = 200,statusMessages = TRUE){
  MX = apply(X = as.matrix(candidates),
             MARGIN = 1,
             FUN = function(par){
               dim = length(par)

               newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.

               # Generalized simmulated Annealing
               tryCatch(expr = {
                 if (statusMessage) print("GenSA")
                 res = GenSA::GenSA(par = par,
                                    fn = fn,
                                    lower = lower,
                                    upper = upper,
                                    control = list("maxit" = maxit,
                                                   "max.time" = 10))
                 newVals = c(newVals,res$par,-res$value,3)
               },
               warning = function(w){
                 if (statusMessage) print(paste0("GenSA::GenSA WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessage) print(paste0("GenSA::GenSA ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Simulated Annealing - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessage) print("stats::optim SANN")
                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "SANN",
                                    control = list(
                                      temp = 1, #Starting temperature
                                      tmax = 5 #number of function evaluations per temperature step
                                    ))
                 if (!FALSE %in% (res$par <= lower & res$par >= upper)){
                   newVals = c(newVals,res$par,-res$value,4)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessage) print(paste0("stats::optim Simulated annealing WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessage) print(paste0("stats::optim Simulated annealing ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Quasi Newton algorithm with approximated Hessian matrix - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessage) print("stats::optim BFGS")
                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "BFGS",
                                    gr = NULL, # if no analytical gradient function is provided, a finite difference method is used
                                    control = list(
                                      maxit = maxit,
                                      reltol = 1e-8
                                    ))
                 if (!FALSE %in% (res$par <= upper & res$par >= lower)){
                   newVals = c(newVals,res$par,-res$value,5)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessage) print(paste0("stats::optim BFGS WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessage) print(paste0("stats::optim BFGS ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Quasi Newton algorithm with approximated Hessian matrix - BOUNDED
               tryCatch(expr = {
                 if (statusMessage) print("Rvmmin")
                 res = Rvmmin::Rvmmin(par = par,
                                      fn = fn,
                                      gr = "grcentral", # "grfwd" = finite difference forward gradient (numerical gradient)
                                      lower = lower,
                                      upper = upper,
                                      maxit = maxit) #Stop computation, if estimators are out of bounds
                 newVals = c(newVals,res$par,-res$value,6)
               },
               warning = function(w){
                 if (statusMessage) print(paste0("RVMMIN WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){

                 if (statusMessage) print(paste0("RVMMIN ERROR: ", e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})


               # Hooke-Jeeves algorithm (Derivative free) - BOUNDED
               tryCatch(expr = {
                 if (statusMessage) print("dfoptim::hjkb")
                 res = dfoptim::hjkb(par = par,
                                     fn = fn,
                                     lower = lower,
                                     upper = upper,
                                     control = list())
                 newVals = c(newVals,res$par,-res$value,7)
               },
               warning = function(w){
                 if (statusMessage) print(paste0("DFOPTIM::HJKB WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessage) print(paste0("DFOPTIM::HJKB ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})

               # Nelder Mead (Derivative free) - UNBOUNDED
               tryCatch(expr = {
                 if (statusMessage) print("stats::optim Nelder-Mead")
                 res = stats::optim(par = par,
                                    method = "Nelder-Mead",
                                    fn = fn,
                                    control = control)
                 newVals = c(newVals,res$par,-res$value,8)
                 if (!FALSE %in% (res$par <= upper & res$par >= lower)){
                   newVals = c(newVals,-res$value,res$par,8)
                 }else{
                   newVals = c(newVals,0,rep(0,dim),1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 if (statusMessage) print(paste0("stats::optim Nelder-Mead WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessage) print(paste0("stats::optim Nelder-Mead ERROR: ",e))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               finally = {})


               # Lightweight BFGS - BFGS with box constraints
               tryCatch(expr= {
                 if (statusMessage) print("stats::optim L-BFGS-B")
                 res = stats::optim(par = par,
                                    fn = fn,
                                    method = "L-BFGS-B",
                                    lower = lower,
                                    upper = upper,
                                    control = list(maxit = maxit)
                 )
                 newVals = unname(c(newVals,res$par,-res$value,9))
               },
               warning = function(w){
                 if (statusMessage) print(paste0("stats::optim L-BFGS-B WARNING: ",w))
                 newVals = c(0,rep(0,dim),1) #Method code 1 is a (non informative) gap filler - Will be later erased.
                 return()
               },
               error = function(e){
                 if (statusMessage) print(paste0("stats::optim L-BFGS-B ERROR: ",e))
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
  MX_df = MX_df[,order(MX_df[Nrows,],
                       decreasing = TRUE)]
  return(MX_df)
}
