#' @export
myOpt = function(candidates,fn,maximum = TRUE,maxit=200,bounds){
  MX = apply(X = as.matrix(candidates),
             MARGIN = 1,
             FUN = function(row){
               if(!maximum) fn = function(x) -fn(x)
               fn_minus = function(x) -fn(x) # for minimum search
               dim = length(row)
               par = row[1:length(dim)]
               newVals = c(0,0,1) #Method code 1 is a (non informative) gap filler - Will be later erased.

               if(dim==1){
                 #Combination of golden section search and parabolic interpolation.
                 del = abs(0.5*par)
                 tryCatch(expr = {
                   res = stats::optimize(f = fn,
                                         lower = par - del,
                                         upper = par + del,
                                         maximum = TRUE)
                   print(paste0("result = ",res))
                   newVals = c(newVals,res$objective,res$maximum,2) #Last value is method code.
                 },
                 warning = function(w){
                   print(paste0("stats::optimize WARNING: ",w))
                   newVals = c(newVals,0,0,1)
                   return()
                 },
                 error = function(e){
                   print(paste0("stats::optimize ERROR: ",e))
                   return()
                 })
               }

               # Generalized simmulated Annealing
               tryCatch(expr = {
                 res = GenSA::GenSA(par = par,
                                    fn = fn_minus,
                                    lower = bounds[1],
                                    upper = bounds[2],
                                    control = list(maxit = maxit))
                 newVals = c(newVals,-res$value,res$par,3)
               },
               warning = function(w){
                 print(paste0("GenSA::GenSA WARNING: ",w))
                 newVals = c(newVals,0,0,1)
                 return()
               },error = function(e){
                 print(paste0("GenSA::GenSA ERROR: ",e))
                 newVals = c(newVals,0,0,1)
                 return()
               })
               # Simulated Annealing - UNBOUNDED
               tryCatch(expr = {
                 res = stats::optim(par = par,
                                    fn = fn_minus,
                                    method = "SANN",
                                    control = list(
                                      temp = 1, #Starting temperature
                                      tmax = 5 #number of function evaluations per temperature step
                                    ))
                 if(res$par <= domain_last & res$par >= domain_first){
                   newVals = c(newVals,-res$value,res$par,4)
                 }else{
                   newVals = c(newVals,0,0,1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 print(paste0("stats::optim SANN WARNING: ",w))
                 newVals = c(newVals,0,0,1)
                 return()
               },
               error = function(e){
                 print(paste0("stats::optim SANN ERROR: ",e))
                 newVals = c(newVals,0,0,1)
                 return()
               })

               # Quasi Newton algorithm with approximated Hessian matrix - UNBOUNDED
               tryCatch(expr = {
                 res = stats::optim(par = par,
                                    fn = fn_minus,
                                    method = "BFGS",
                                    gr = NULL, # if no analytical gradient function is provided, a finite difference method is used
                                    control = list(
                                      maxit = maxit,
                                      reltol = 1e-8
                                    ))
                 if(res$par <= domain_last & res$par >= domain_first){
                   newVals = c(newVals,-res$value,res$par,5)
                 }else{
                   newVals = c(newVals,0,0,1)#method code 1 is a (non informative) gap filler - Will be erased later
                 }
               },
               warning = function(w){
                 print(paste0("stats::optim BFGS Warning: ",w))
                 newVals = c(newVals,0,0,1)
                 return()
               },
               error = function(e){
                 print(paste0("stats::optim BFGS ERROR: ",e))
                 newVals = c(newVals,0,0,1)
                 return()
               })

               # Quasi Newton algorithm with approximated Hessian matrix - BOUNDED
               tryCatch(expr = {
                 res = Rvmmin::Rvmmin(par = par,
                                      fn = fn_minus,
                                      gr = "grcentral", # "grfwd" = finite difference forward gradient (numerical gradient)
                                      lower = bounds[1],
                                      upper = bounds[2],
                                      maxit = maxit) #Stop computation, if estimators are out of bounds
                 newVals = c(newVals,-res$value,res$par,6)
               },
               warning = function(w){
                 print(paste0("Rvmmin::Rvmmin Warning: ",w))
                 newVals = c(newVals,0,0,1)
                 return()
               },
               error = function(e){
                 print(paste0("Rvmmin::Rvmmin ERROR: ",e))
                 newVals = c(newVals,0,0,1)
                 return()
               })

               # Lightweight BFGS - BFGS with box constraints
               tryCatch(expr = {
                 res = stats::optim(par = par,
                                    fn = fn_minus,
                                    method = "L-BFGS-B",
                                    lower = bounds[1],
                                    upper = bounds[2],
                                    control = list(maxit = maxit)
                 )
                 newVals = unname(c(newVals,-res$value,res$par,7))
               },
               warning = function(w){
                 print(paste0("stats::optim L-BFGS-B Warning: ",w))
                 newVals = c(newVals,0,0,1)
                 return()
               },
               error = function(e){
                 print(paste0("stats::optim L-BFGS-B ERROR: ",e))
                 newVals = c(newVals,0,0,1)
                 return()
               })
               return(newVals)
             })
  names = c("value",names(candidates),"methodCode")
  methods = c("invalid","optimize","GenSA","SA","BFGS","Rvmmin","L-BFGS-B")
  cols = length(names)
  MX_df = data.frame(t(matrix(data = MX,
                              nrow = cols)))
  colnames(x = MX_df) = names
  MX_df = cbind(MX_df[,-cols],data.frame("method" = methods[MX_df$methodCode]))
  MX_df = MX_df[!(MX_df$method == "invalid"),] # remove non informative rows
  return(MX_df)
}
