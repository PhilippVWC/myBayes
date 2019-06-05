#### function library ####
#### Dr. J.C. Lemm and P.v.W.Crommelin ###

#### Rcpp test ####
Rcpp::cppFunction('int add2(int x, int y, int z) {
  int sum = x + y + z;
            return sum;
            }')
###
#' @export
collapse = function(vec){
  deltaVec = vec - c(vec[1]+1,vec[1:(length(vec)-1)])
  indices = which(deltaVec != 0)
  df = data.frame(values = vec[indices], indices = indices)
  return(df[order(df$values,
                  decreasing = TRUE),])
}
#' @export
myOpt = function(candidates,fn,maximum = TRUE,maxit=200,bounds){
  MX = apply(X = as.matrix(candidates),
             MARGIN = 1,
             FUN = function(row){
               if(!maximum) fn = function(x) -fn(x)
               fn_minus = function(x) -fn(x) # for minimum search
               dim = length(row)
               par = row[1:length(dim)]
               newVals = c(0,0,1) #method code 1 is a (non informative) gap filler - Will be erased later

               if(dim==1){
                 # combination of golden section search and parabolic interpolation
                 del = abs(0.5*par)
                 res = optimize(f = fn,
                                lower = par - del,
                                upper = par + del,
                                maximum = TRUE)
                 newVals = c(newVals,res$objective,res$maximum,2) # last value is method code
               }

               # Generalized simmulated Annealing
               res = GenSA(par = par,
                           fn = fn_minus,
                           lower = bounds[1],
                           upper = bounds[2],
                           control = list(maxit = maxit))
               newVals = c(newVals,-res$value,res$par,3)

               # Simulated Annealing - UNBOUNDED
               res = optim(par = par,
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

               # Quasi Newton algorithm with approximated Hessian matrix - UNBOUNDED
               res = optim(par = par,
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

               # Quasi Newton algorithm with approximated Hessian matrix - BOUNDED
               res = Rvmmin(par = par,
                            fn = fn_minus,
                            gr = "grcentral", # "grfwd" = finite difference forward gradient (numerical gradient)
                            lower = bounds[1],
                            upper = bounds[2],
                            maxit = maxit) #Stop computation, if estimators are out of bounds
               newVals = c(newVals,-res$value,res$par,6)

               # Lightweight BFGS - BFGS with box constraints
               res = optim(par = par,
                           fn = fn_minus,
                           method = "L-BFGS-B",
                           lower = bounds[1],
                           upper = bounds[2],
                           control = list(maxit = maxit)
               )
               newVals = unname(c(newVals,-res$value,res$par,7))
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





#' @export
kroneckerDelta = function(x,x0){
  if (x == x0) 1
  else 0
}
#' @export
betaBin = function(x,alpha,beta,n){  # Beta binomial distribution
  choose(n,x) * beta(alpha+x,beta+n-x) / beta(alpha,beta)
}
#' @export
blurrVector = function(vec,lambda,sigmaRel = 0.1){
  n = length(vec) - 1
  mu = which(vec == max(vec)) - 1
  N_samp = 1/sigmaRel^2 - 1
  y_s = min(1/(N_samp*mu/n),1/(N_samp*(1-mu/n))) # ensures that (alpha,beta) == 1 <==> lambda == 0
  lambda_t = 1-(1-y_s)*lambda #transform lambda
  alpha = max(1,lambda_t * mu/n * N_samp)
  beta = max(1,lambda_t * (1-mu/n) * N_samp)
  x = 0:n
  result = sapply(X = x,
                  FUN = function(x) (1-lambda)*kroneckerDelta(x = x,
                                                              x0 = mu) + lambda*betaBin(x = x,
                                                                                        alpha = alpha,
                                                                                        beta = beta,
                                                                                        n = n))
  return(result)
}




#' @export
Heaviside = function(x,x0){
  if (x<x0) return(0)
  else return(1)
}
#' @export
rectFun = function(x,a,b,z){
  result = z * (Heaviside(x = x,
                          x0 = a) - Heaviside(x = x,
                                              x0 = b))
  return(result)
}




# Metropolis Hastings Algoritm. Create a random variable, which is distributed according to a custom probability distribution.
#' @export
MH = function(T = 1,x0,trKern=dnorm,sampler,N=1000,bounds){
  a = bounds[1] # left bound of domain to sample from
  b = bounds[2] # right bound
  chain = rep(0,N) # Markov chain
  chain[1] = x0
  for (i in 2:N){
    x = chain[i-1]
    y = x + T*runif(n = 1,
                    min = -(x-a),
                    max = b-x) # sample next candidate within bounds
    A = min( 1 , (sampler(y)*trKern(y,x) / (sampler(x)*trKern(x,y)) ) )
    if (A == 1) { chain[i] = y }
    else {
      if ( A > runif(n = 1,
                     min = 0,
                     max = 1)){
        chain[i] = y
      }
      else {
        chain[i] = x
      }
    }
  }
  return(chain)
}

#' @export
CeBePa = function(mean,sigmaRel){
  N = 1/sigmaRel^2 - 1.0
  alpha = max(1, mean*N)
  beta = max(1, (1-mean)*N)
  return(c(alpha,beta))
}

#' @export
betaDist = function(x,alpha,beta){
  if (x<=1 & x>=0 & alpha > 0 & beta > 0) {
    B = gamma(alpha+beta)/(gamma(alpha)*gamma(beta))
    return(B * x^(alpha-1)*(1-x)^(beta-1))
  }
}

#' @export
vecToVal = function(vec,x_lower=0,x_upper=1){
  N = length(vec)
  dx = (x_upper-x_lower)/(N-1)
  ref = seq(from = x_lower,
            to = x_upper,
            by = dx)
  #s = sum(vec)
  #if (s!=1) print(paste0("Warning, vector is not normed to one. Sum = ",s))
  return(ref%*%vec)
}

#' @export
valToVec = function(val,N=10,x_lower=0,x_upper=1){
  if (val > x_upper | val < x_lower ) print("warning: val is not within bounds of x_lower and x_upper")
  if (N>1){
    dx = (x_upper-x_lower)/(N-1)
    ref = seq(from = x_lower,
              to = x_upper,
              by = dx)
    diff = abs(ref - rep(val,N))
    minInd = which(diff == min(diff))
    result = rep(0,N)
    result[minInd[1]] = 1 #if value lies exactly between to grid points -> condense to left one
    return(result)
  }else{
    print(paste0("N needs to be greater than 1. N is currently set to ",N))
    return(NA)
  }
}

#' @export
packageManager = function(necessaryPackages){
  installedPackages = installed.packages()[,"Package"]
  missingPackages = necessaryPackages[!(necessaryPackages %in% installedPackages)]
  if (length(missingPackages) > 0){
    cat("Installation of the following packages:", missingPackages ,"\n")
    cat("This may take a while\n")
    install.packages(missingPackages)
    installedPackages = installed.packages()[,"Package"]
  }
  successfullyInstalled = missingPackages[missingPackages %in% installedPackages]
  notInstalled = missingPackages[!(missingPackages %in% installedPackages)]
  if (length(successfullyInstalled) > 0){
    cat("The following packages were installed:",successfullyInstalled,"\n")
  }
  if (length(notInstalled) > 0) {
    cat("The following packages were not installed:",notInstalled,"\n")
  }
  lapply(necessaryPackages,require,character.only=TRUE) # load library after installation
  return("packageManager finished")
}
#' @export
discretizeMap = function(N,map,xlim){
  #--- parameters ---#
  #N      number of discrete points (equal to range of matrix)
  #map    map to be discretized. Should be reduced to function of one variable
  #xlim   vector with two entries containing the bounds of the problem
  x_lower = xlim[1]
  x_upper = xlim[2]
  dx = (x_upper-x_lower)/(N-1)
  domain = seq(from = x_lower,
               to = x_upper,
               by = dx)
  codomain = sapply(X = domain,
                    FUN = function(x) map(x))

  helperMatrix = matrix(codomain,
                        ncol = 1) %*% rep(1,N)

  # "condense" function values to grid points
  indices = apply(X = helperMatrix,
                  MARGIN = 1,
                  FUN = function(x) {
                    dist = abs(x - domain) # distance from grid points
                    index = which(dist == min(dist)) # index of
                    return(index)
                  })
  A = matrix(rep(0,N*N),
             ncol = N,
             nrow = N)
  for (i in 1:N){
    A[indices[i],i] = 1
  }
  return(A)
}


#' @export
max1d = function(vec,epsilon=1.0,maximum=TRUE,maxRows=NA){
  if ( is.vector(vec) & epsilon>=0 & epsilon<=1 ){
    l = length(vec)
    if(maximum){
      vec = data.frame(sort(x = vec,
                            decreasing = TRUE,
                            index.return = TRUE))
    }else{
      vec = data.frame(sort(x = vec,
                            decreasing = FALSE,
                            index.return = TRUE))
    }
    names(x = vec) = c("value","indices")
    if (is.na(maxRows)){
      upper_index = round(epsilon*l,0)
    }else{
      upper_index = min(round(epsilon*l,0),maxRows)
    }
    candidates = vec[1:upper_index,]
    return(candidates)
  }else{
    print("Wrong input data")
    return(NULL)
  }
}


# Function for 2D matrices returning extremal values and its corresponding indices
#' @export
max2d = function(mat,epsilon=1.0,maximum=TRUE,maxRows=NA){
  # PARAMETERS:
  # epsilon       relative tolerance level in [0,1].
  #               Zero: consider only maximum values (possibly degenerated)
  #               One: consider every value
  # mat           Matrix at hand
  # maximum       Boolean value indicating, wether to search for a minimum or a maximum
  # maxRows       Upper bound for number of results
  if ( is.matrix(mat) & length(dim(mat))==2 & epsilon>=0 & epsilon<=1 ){
    dim = dim(mat)
    dimX = dim[1]
    mx = max(mat)
    mn = min(mat)
    tol = abs((mx - mn))*epsilon # convert to absolute tolerance value
    if (maximum == TRUE){
      values = mat[mat >= (mx-tol) & mat <= mx]
      indices = which(mat >= (mx-tol) & mat <= mx)
    } else {
      values = mat[mat >= mn & mat <= (mn+tol)]
      indices = which(mat >= mn & mat <= (mn+tol))
    }
    x_indices = indices %% dimX
    x_indices[x_indices == 0] = dimX        #= dimX, because otherwise values within the last row are interpreted to be in the first row
    y_indices = (indices - 0.1) %/% dimX    #-0.1, because otherwise values within the last row are interpreted to be in the next column
    y_indices = y_indices + 1               #= 1, because R indexes start with number 1
    result = data.frame("value"=values, "RowIndex"=x_indices, "ColIndex"=y_indices)
    nRow = nrow(result)
    if ( !is.na(maxRows) & nRow>maxRows ) {
      print("###########################################################################")
      print("### Function <<max2d>>:")
      print(paste0("### Upper bound hit. Available results: ",nRow))
      print("### Increase <<maxRows>> to increase number of results.")
      print("###########################################################################")
      if (maximum == TRUE){
        result = result[order(result$MaxVal,
                              decreasing = TRUE),]
        return(result[1:maxRows,])
      } else {
        result = result[order(result$MaxVal,
                              decreasing = FALSE),]
        return(result[1:maxRows,])
      }
    }else{
      return(result)
    }
  } else {
    print("Wrong input data")
    return(data.frame("value"=NA, "RowIndex"=NA, "ColIndex"=NA))

  }
}

# Takes X and Y coordinates and creates a grid in the form of a vector
# X and Y are both 1D, surface is 2D
# output is a matrix with 3 columns containing the coordinates X,Y,Z as column vectors
#' @export
coordinates = function(X,Y,surface){
  N_X = length(X)
  N_Y = length(Y)
  Z = c(surface)
  X = rep(X,N_Y)
  Y = c(t(matrix(Y) %*% rep(1,N_X)))
  return(matrix(c(X,Y,Z),ncol=3))
}
#' @export
coordinates2 = function(X,Y,surface){
  N_X = length(X)
  N_Y = length(Y)
  X = matrix(rep(X,N_Y),nrow =  N_X)
  Y = c(t(matrix(Y) %*% rep(1,N_X)))
  return(matrix(c(X,Y,Z),ncol=3))
}

#packageManager("inline")
#### general symmetric map (C extension) ####
sig = c(N="integer",x0="numeric",r="numeric",alpha="numeric",skipFirst="integer")
body = "
#include <math.h>
#include <Rinternals.h>


// convert every SEXP object into a c datatype
double r_ = asReal(r);
double a_ = asReal(alpha);
double x0_ = asReal(x0);
int N_ = asInteger(N);
int skipFirst_ = asInteger(skipFirst);


// Allocate an output vector
SEXP X = PROTECT(allocVector(REALSXP,N_));

// create an Array pointing to that output vector
double *X_;
X_ = REAL(X);


//Manipulate the array
if (skipFirst_)  X_[0] = r_*(1.0 - pow(2.0 * fabs(x0_-0.5),a_));
if (!skipFirst_) X_[0] = x0_;

int i;
for ( i=0 ; i<(N_-1) ; i++ ){
    X_[i+1] = r_*(1.0 - pow(2.0 * fabs(X_[i]-0.5),a_));
}

UNPROTECT(1);
//return that array
return X;
"
#' @export
gensymMap_iter_c =  inline::cfunction(sig = sig,
                              body = body,
                              verbose = FALSE,
                              convention = ".Call",
                              language = "C")
#setCMethod(f = "gensymMap_iter_c",
#           sig = sig,
#           body = body,
#           verbose = FALSE)

#test = inline::cfunction(sig = c(x="integer"),
#                         body = 'return(x);')
# dynamisches Rauschen (intrinsisches Rauschen)
# logistische Abbildung mit dynamischem (additiven, betaverteiltem) Rauschen
# (0<x0<1, 0<a<4, 0 < sig)
#' @export
logMap_dyn = function(N,a,x0,sigmaRel=0.1){
  X = c(x0)
  if(N>=1){
    for (i in 1:N){
      if( (X[i] > 0) & (X[i] < 1) & (sigmaRel > 0) ) {
        #Var = max(0,min(sigma^2,x*(1-x))) # Varianz aus absolutem Sigma, muss bei maximaler Varianz gekappt werden
        #Var = x*(1-x)*(sigmaRel^2) # Varianz aus relativem Sigma
        #nu = x*(1-x)/Var - 1.0 # "sample size" aus Mittelwert (mu=x) und Varianz = sigma^2
        nu = 1/sigmaRel^2 - 1.0 # "sample size"
        alpha = max(1e-10, X[i] * nu) # alpha aus sample size (nu) und Mittelwert (mu=x)
        beta = max(1e-10, (1-X[i]) * nu) # beta aus sample size (nu) und Mittelwert (mu=x)
        X[i] = rbeta(1,alpha,beta) # neues x mit Rauschen um mu=x
        #cat("i=",i,"sigma=",sigma,"nu=",nu,"alpha=",alpha, "beta=",beta," x=",x,"\n")
      }
      X[i+1] = a*X[i]*(1-X[i])
    }
  }
  return(X[-1])
}

#### Likelihoods ####
#' @export
LLik = function(X,Y,sigma){
  # strict result = -sum( 0.5*((Y-X)/sigma)^2 ) + length(X)*log(1/(sqrt(2*pi)*sigma))
  return(-sum( 0.5*((Y-X)/sigma)^2 ))
}

#not a likelihood function
#' @export
Lik = function(X,Y,sigma){
  return(-LLik(X = X,
               Y = Y,
               sigma = sigma))
}

#' @export
Lik_gensymMap = function(alpha,r,x0,Y,sigma,N_discr = 0){
  # N_discr only for compatibility reasons
  n = length(Y)
  X = gensymMap_iter(N = n,
                     x0 = x0,
                     r = r,
                     alpha = alpha,
                     N_discr = N_discr,
                     skipFirst = TRUE)
  L = exp(LLik(X = X,
               Y = Y,
               sigma = sigma)/(sqrt(2*pi)*sigma))
  return(L)
}

#' @export
Lik_dGensymMap = function(alpha,r,x0,Y,sigma,N_discr){
  n = length(Y)
  X = dGensymMap_iter(N = n,
                      x0 = x0,
                      r = r,
                      alpha = alpha,
                      N_discr = N_discr,
                      skipFirst = TRUE)
  L = exp(LLik(X = X,
               Y = Y,
               sigma = sigma)/(sqrt(2*pi)*sigma))
  return(L)
}
#Likelihood function for dynamic noise (beta distributed)
#returns a modell Vector X of dimension Dim(Y)-1
#with Y being the data vector, respectively
#' @export
Lik_dyn = function(a,x0,Y,sigmaRel){
  ly = length(Y)
  n = ly-1
  X = a*Y[2:ly]*(1-Y[2:ly])
  Lik = exp(-sum( ( 0.5*(X-Y[2:ly])/(sigmaRel) )^2  ))/(sigmaRel*sqrt(2*pi))^n
  return(Lik)
}

#' @export
discr = function(val,N_discr,lower = 0.0,upper = 1.0){
  dx = (upper-lower)/(N_discr-1)
  ref = seq(from = lower,
            to = upper,
            by = dx)
  diff = abs(ref - rep(val,N_discr))
  factor = which( diff == min(diff) )[1] # if value is exactly between two grid points --> condense to left one
  return((factor-1)*dx)
}

#' @export
discreteMap_iter = function(N,x0,mat,skipFirst = TRUE){
  x0 = valToVec(val = x0,
                N = dim(mat)[2])
  if (skipFirst) x0 = mat %*% x0
  Series = matrix(data = rep(x0,N),
                  nrow = length(x0),
                  ncol = N,
                  byrow = FALSE)
  for (i in 2:N){
    Series[,i] = mat %*% Series[,i-1]
  }
  result = apply(X = Series,
                 MARGIN = 2,
                 FUN = function(colVec){
                   vecToVal(vec = colVec)
                 })
  return(result)
}

#dGensymMap_iter = function(N,x0,r,alpha,N_discr,skipFirst = TRUE,lower = 0.0, upper = 1.0){
#  Series = rep(x0,N)
#  if(skipFirst){
#    Series[1] = discr(gensymMap(x = x0,
#                                r = r,
#                                alpha = alpha),
#                      N_discr = N_discr,
#                      lower = lower,
#                      upper = upper)
#  }
#  if(N>1){
#    for(i in 2:N){
#      Series[i] = discr(val = gensymMap(x = Series[i-1],
#                                        r = r,
#                                        alpha = alpha),
#                        N_discr = N_discr,
#                        lower = lower,
#                        upper = upper)
#    }
#  }
#  return(Series)
#}

# general symmetric map
#' @export
gensymMap_iter = function(N,x0,r,alpha,N_discr = 0,skipFirst=TRUE){
  X = rep(x0,N)
  if(skipFirst){
    X[1] = if(N_discr != 0){
      round(gensymMap(x = x0,
                      r = r,
                      alpha = alpha)*N_discr)/N_discr
    }else{
      gensymMap(x = x0,
                r = r,
                alpha = alpha)
    }
  }
  if(N>1){
    if(N_discr != 0){
      for (i in 2:N) X[i] = round(gensymMap(x = X[i-1],
                                            r = r,
                                            alpha = alpha)*N_discr)/N_discr
    }else{
      for (i in 2:N) X[i] = gensymMap(x = X[i-1],
                                      r = r,
                                      alpha = alpha)
    }
  }
  return(X)
}



# logistic map
#' @export
logMap_iter = function(N,r,x0,skipFirst=TRUE){
  if (N >= 1) {
    X=rep(0.0,N)
    if(skipFirst) X[1] = 4*r*x0*(1-x0)
    else X[1] = x0
    for (n in 2:N){
      X[n] = 4*r*X[n-1]*(1-X[n-1])
    }
    return(X)
  } else {
    cat("invalid number N =",N,"\n")
    return(NULL)
  }
}



#shifts every element of a vector >>offset<< times to the right
#' @export
shift = function(vec,offset){
  return(c(vec[(1+offset):length(vec)],vec[1:offset]))
}

#' @export
simplePlot = function(main="simplePlot",
                      x = NULL,
                      y = NULL,
                      xlim = NULL,
                      ylim = NULL,
                      xlab = "x",
                      ylab = "y",
                      col = "black",
                      bg = NULL,
                      type = "l",
                      sub = NULL,
                      col.axis = COLOR_AXIS,
                      col.lab = NULL,
                      col.main = NULL,
                      col.sub = NULL,
                      fg = COLOR_AXIS,
                      pch = 16,
                      cex = 1,
                      log = "",
                      leg = NULL,
                      position_leg = "topleft",
                      inset_leg = 0.05,
                      col_leg = NULL,
                      lty_leg = 1,
                      lwd_leg = 2,
                      title_leg = NULL,
                      bg_leg = COLOR_BG_LEG,
                      nx_grid = NX_GRID,
                      ny_grid = NY_GRID,
                      col_grid = COLOR_GRID,
                      lty_grid = POINTS,
                      lwd_grid = LWD_GRID,
                      las = NULL){
  plot(main = main,
       sub = sub,
       x = x,
       y = y,
       xlab = xlab,
       ylab = ylab,
       col = col,
       type = type,
       bg = bg,
       xlim = xlim,
       ylim = ylim,
       col.axis = col.axis,
       col.lab = col.lab,
       col.main = col.main,
       col.sub = col.sub,
       fg = fg,
       pch = pch,
       cex = cex,
       log = log,
       las = las)
  #grid(nx = NX_GRID,
  #     ny = NY_GRID,
  #     col = COLOR_GRID,
  #     lty = POINTS,
  #     lwd = LWD_GRID,
  #     equilogs = TRUE)
  grid()

  if (!is.null(leg)) {
    legend(position_leg,
           inset = inset_leg,
           legend = leg,
           col = col_leg,
           lty = lty_leg,
           lwd = lwd_leg,
           title = title_leg,
           bg = bg_leg)
  }
}

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
           abline(v = tick,
                  lwd = lwd,
                  lty = lty,
                  col = col)
         })
  sapply(X = yticks,
         FUN = function(tick){
           abline(h = tick,
                  lwd = lwd,
                  lty = lty,
                  col = col)
         })
}


#### Abbildungen/Maps ####
#' @export
logMap = function(x,r){
  return(4*r*x*(1-x))
}

#### Variante für Testzwecke
#' @export
logisticMap2 = function(x,r){
  return(4*r*(x-x*x))
}

#### nach Iteration gekappte Variante
#' @export
logisticMapPostcut = function(x,r){
  return(min(1,max(0,4*r*x*(1-x))))
}

#### vor Iteration gekappte Variante
#' @export
logisticMapPrecut = function(x,r){
  x=min(1,max(0,x))
  return(min(1,max(0,4*r*x*(1-x))))
}

#### kubische Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
cubeMap = function(x,r){
  return(r*(1-8*abs(x-0.5)^3)) # = r*(1-(2*abs(x-1/2))^3)
}

#### quartische Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
quartMap = function(x,r){
  return(r*(1-16*(x-0.5)^4))
  # r*(1-(2*(x-1/2))^4) =  r*(1-16*(x-0.5)^4)
  # = r*(1-16*(x^2 -x +1/4)^2) = r*(1 - 16*(x^4 +x^2 +1/16 -2*x^3 +x^2/2 -x/2))
  # = r*(1 -1 -2^3*x -3*2^3*x^2 -2^5*x^3 -2^4*x^4)
  # = 2^3*r*(x -3*x^2 -2^2*x^3 -2*x^4) = 8*r*(x -3*x^2 +4*x^3 -2*x^4)
  # = 8*r*(x -x^2 -2*x^2 +2*x^3 +2*x^3 -2*x^4)
  # = 8*r*x*(1-x)*(1-2*x+2*x^2) (hier ist die Nullstelle bei x=0 und x=1 direkt abzulesen)
  # = 8*r*(x-x^2)*(1-2*x+2*x^2)
  # = logisticMap(x,r)*2*(1-2*x+2*x^2)

}

#### Lorenz-Abbildung (Wurzel)
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
lorenzMap = function(x,r){
  return(r*(1-sqrt(abs(x-0.5))*sqrt(2))) # = r*(1-(2*abs(x-1/2))^(1/2))
}
####

#### Kubische Wurzel-Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
cuberootMap = function(x,r){
  return(r*(1-(2*abs(x-0.5))^(1/3)))
}

#### Sinus-Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
sineMap = function(x,r){
  return(r*sin(pi*x))
}

#### tent map (Zelt-Abbildung)
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
tentMap = function(x,r){
  return(r*(1-2*abs(x-0.5)))
}

#### general symmetric map (S. Sprott: Chaos and Time-Series Analysis, S. 39)
# (0 <= x0 <= 1, 0 <= r <= 1, logisticMap: alpha = 2, lorenzMap: alpha = 1/2 etc. )
#' @title General symmetric map
#' @description a generalization of a one dimensional map, which incorporates (among others) the lorenz-, logistic and cubeMap
#' @param x input value for one iteration
#' @author J.C. Lemm, P.v.W. Crommelin
#' @references S. Sprott, Chaos and Time-series analysis
#' @param r control parameter
#' @param alpha exponent
#' @return The result of one single iteration
#' @export
gensymMap = function(x,r,alpha){
  return(r*(1-(2*abs(x-0.5))^alpha))
}

#### Error-Functions
## 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## and the inverses
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)


#### Gaussian white chaotic map (S. Sprott: Chaos and Time-Series Analysis, S. 420)
#' @export
gaussianWhiteChaoticMap = function(x,r){
  return(r*erfinv(1-2*erf(x/r)))
}

#### Binary shift map
#' @export
binaryShiftMap = function(x,r){
  return((2*r*x)%%1)
}

#plotBifurcation (funct=binaryShiftMap, xStart=c(0.2,0.4), nrLevels=100, nPlot = nLevels, nStart=nLevels)
#plotIteratedMaps(mapNames=list("binaryShiftMap"),r=r0,nList=1:5)
#plotTimeSeries(mapName=list("binaryShiftMap"),x=0.2,r=1,n=200)
#plotTimeSeries(mapName=list("binaryShiftMap"),x=0.2,r=1,n=20,cobweb=TRUE)


# Anzahl Level f?r diskrete Funktionen
nLevels = 101  # n > 1

#### discrete logistic map
# Tabelle diskrete logistische Abbildung (bei festem r so schneller)
# Dlogvals = round(logisticMap((0:(nLevels-1))/(nLevels-1),r=r0)*(nLevels-1))
# --- Diskrete Funktion
#' @export
discreteLogisticMap =function(x,r,m=nLevels){
  #  0 < x <= 1, 0 <= r <= 1
  #return(Dlogvals[round(x*nLevels)]/nLevels)
  return(round(logisticMap(round(x*(m-1))/(m-1),r)*(m-1))/(m-1))
}
#discreteLogisticMap(x=seq(0,1,length.out = nLevels),r0)



#### discrete random map
# Tabelle diskrete Zufallsfunktion
#Drvals = pmax(0,ceiling(runif(nLevels)*nLevels)-1)/(nLevels-1) # Levels von 0 bis nLevels-1 jew. zwischen 0 und 1
Drvals = pmax(0,ceiling(runif(nLevels)*nLevels)-1) # Levels von 0 bis nLevels-1
Drvals = Drvals - min(Drvals) # Starte immer bei 0
# --- Zufallsfunktion
#' @export
discreteRandomMap = function(x,r){
  #  0 < x <= 1, 0 <= r <= 1
  return(ceiling(r*Drvals[pmax(1,ceiling(x*nLevels))])/(nLevels-1))
}




#### Logit-Normal-Verteilung
# Wird parametrisiert durch die Parameter mu und sd der zugrundeliegenden Normalverteilung
#  F?r die Momente der Logit-Normal-Verteilung gibt es allerdings keine analytische L?sung.
#  Gesucht w?re ein Funktion, die aus vorgegebenem Mittelwert und SD,
#  die entsprechenden mu, sd der Normalverteilung berechnet.
#' @export
rlogitnormal = function(...){
  x = exp(rnorm(...))
  x / (1+x)
}
#samp = rlogitnormal(1000, 0, 1)
#summary(samp);1/(1+exp(0))
#samp = rlogitnormal(1000, 1, 1)
#summary(samp);1/(1+exp(-1))
#samp = rlogitnormal(1000, 1, 10)
#summary(samp);1/(1+exp(-1))



#### probabilistische Einschritt-Abbildungen

#### allgemeines dynamisches Zustandsrauschen
#' @export
noisyStateMap = function(x,r
                         ,funct  = logisticMap
                         ,rfunct = function(z) rnorm(n=1,mean=z,sd=SD0) ){
  return(rfunct(funct(x,r)))
}
# Examples:
#noisyStateMap(x0,r0)
#noisyStateMap(x0,r0,rfunct = function(z) rnorm(n=1,mean=z,sd=1)) # Bsp. Normalverteilung
#noisyStateMap(x0,r0,rfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=1)))) # Bsp. Normalverteilung mit Postcut
#noisyStateMap(x0,r0,rfunct = function(z) rbeta_MuSDrel(n=1, x=x0, sdRel=SDrel0)) # Bsp. Betaverteilung
#noisyStateMap(x0,r0,rfunct = rbeta_MuSDrel0) # Bsp. Betaverteilung mit festem sdRel=SDrel0

#### Spezialfall: dynamisches Zustandsrauschen mit Betaverteilung
#' @export
noisyStateMap_BetaMuSDrel0 = function(x,r,funct=logisticMap) noisyStateMap(x, r, funct, rfunct=rbeta_MuSDrel0)
# Examples:
#noisyStateMap_BetaMuSDrel0(x0,r0)
#noisyStateMap_BetaMuSDrel0(x0,r=1)

#### allgemeines dynamisches Parameterrauschen
#' @export
noisyParmMap = function(x,r
                        ,funct  = logisticMap
                        ,prfunct = function(z) rnorm(n=1,mean=z,sd=SD0) ){
  return(funct(x,prfunct(r)))
}
# Examples:
#noisyParmMap(x0,r0)
#noisyParmMap(x0,r0,prfunct = function(z) rnorm(n=1,mean=z,sd=1))
#noisyParmMap(x0,r0,prfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=1))))


#### Spezialfall: dynamisches Parameterrauschen mit Betaverteilung
#' @export
noisyParmMap_BetaMuSDrel0 = function(x,r,funct=logisticMap) noisyParmMap(x, r, funct, prfunct=rbeta_MuSDrel0)
# Examples:
#noisyParmMap_BetaMuSDrel0(x0,r0)
#noisyParmMap_BetaMuSDrel0(x0,r=1)


#### allgemeines dynamisches (Zustands- und Parameter-)Rauschen
# (ersetzt noisyStateMap und noisyParmMap)
#' @export
noisyMap = function(x,r
                    ,funct  = logisticMap  # definiert deterministische univariate Abbildung
                    ,rfunct = function(z) rnorm(n=1,mean=z,sd=SD0)  # definiert dynamisches Zustandsrauschen nach Iterationsschritt
                    ,prfunct = function(z) rnorm(n=1,mean=z,sd=SD0) # definiert dynamisches Parameterrauschen
                    ,xrfunct = identity # definiert dynamisches Zustandsrauschen vor Iterationsschritt
){
  #  return(rfunct(funct(x,prfunct(r))))
  return(rfunct(funct(xrfunct(x),prfunct(r)))) # Variante, bei auch auch x vor der Iteration durch xrfunct() verrauscht wird
}
# Examples:
#noisyMap(x0,r0)
#noisyMap(x0,r=r0,funct=logisticMap)
#noisyMap(x0,r=r0,funct=tentMap)


#### Spezialfall mit gekappter Normalverteilung
# (mit vorgegebener SD f?r dynamisches Zustands- und Parameterrauschen)
#' @export
noisyMap_NormalCut0 = function(x,r,funct=logisticMap)
  noisyMap(x,r,funct
           , rfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=SD0)))
           ,prfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=SD0)))
  )
# Examples:
#noisyMap_NormalCut0(x0,r0)
#noisyMap_NormalCut0(x0,r0,funct=tentMap)

#### Spezialfall mit Betaverteilung
# (mit vorgegebener SD f?r dynamisches Zustands- und Parameterrauschen)
#' @export
noisyMap_BetaMuSDrel0 = function(x,r,funct=logisticMap)
  noisyMap(x,r,funct
           , rfunct = rbeta_MuSDrel0
           ,prfunct = rbeta_MuSDrel0
  )
# Examples:
#noisyMap_BetaMuSDrel0(x0,r0)
#noisyMap_BetaMuSDrel0(x0,r0,funct=tentMap)


#### Iteration

#### einfache Iteration
# (auch geeignet f?r dynamische Rauschvarianten)
#' @export
iteratedMap = function(x,r,n
                       ,funct=logisticMap ){
  if(n>=1){
    for (i in 1:n){
      x = funct(x,r)
    }
  }
  return(x)
}
# Examples:
#iteratedMap(x0, r0, 2)
#iteratedMap(x0, r0, 2, funct = logisticMap)
#iteratedMap(x0, r0, 2, funct = function(x,r) noisyMap(x, r, funct=logisticMap))
#iteratedMap(x0, r0, 2, funct = noisyMap_BetaMuSDrel0)



#### Iteration mit Inputrauschen und abschlie?endem Beobachtungsrauschen (Ausgangsrauschen)
#' @export
iteratedMap_ioNoise = function( x
                                ,r
                                ,n
                                ,funct  = logisticMap
                                ,irfunct = identity # input noise function(u) rnorm(n=1,mean=u,sd=SD0)
                                ,orfunct = identity # output noise function(u) rnorm(n=1,mean=u,sd=SD0)
){
  #  return(rfunct(funct(x,r)))
  return(orfunct(iteratedMap(irfunct(x),r,n,funct)))
}
# Examples:
#iteratedMap_ioNoise(x0,r0,1)
#iteratedMap_ioNoise(x0,r0,1,irfunct=function(u) min(1,max(0,rnorm(n=1,mean=u,sd=0.1))))
#iteratedMap_ioNoise(x0,r0,1,irfunct=function(u) rbeta_MuSDrel(n=1,mean=u,sdRel=0.1))
#iteratedMap_ioNoise(x0,r0,1,orfunct=function(u) rnorm(n=1,mean=u,sd=0.1))
#iteratedMap_ioNoise(x0,r0,1,funct=function(y,r) iteratedMap(y,r,n=10,funct=logisticMap))
#iteratedMap_ioNoise(x0,r0,10,funct=logisticMap)
#iteratedMap_ioNoise(x0,r0,10,funct=tentMap)



#### Zeitreihe
# (eindimensionaler, homogener, diskreter hidden Markovprozess)
#' @export
oneDimhHMM = function(x,r,n
                      ,funct=logisticMap # function(x,r)
                      ,irfunct = identity # function(x)
                      ,orfunct = identity# function(x)
){
  y=rep(NA,n)
  x = irfunct(x)
  if(n>=1){
    for (i in 1:n){
      x = funct(x,r)
      y[i] = orfunct(x)
    }
  }
  return(y)
}
# Examples:
#oneDimhHMM(x0,r0,10)
#plot(oneDimhHMM(x0,0.999,100),type="l",col="blue",lwd=1)
#lines(oneDimhHMM(x0,0.999,100,funct=discreteLogisticMap),type="l",col="orange",lty=3)
#lines(oneDimhHMM(x0,0.999,100,funct=tentMap),type="l",col="red",lwd=1,lty=2)
#lines(oneDimhHMM(x0,0.999,100,funct=lorenzMap),type="l",col="green",lwd=2)













