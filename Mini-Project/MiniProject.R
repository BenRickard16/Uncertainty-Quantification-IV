## 1. Find an Interesting Model

# Function we want to emulate
fried <- function(xx)
{
  ##########################################################################
  #
  # FRIEDMAN FUNCTION
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, x3, x4, x5)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  x5 <- xx[5]
  
  term1 <- 10 * sin(pi*x1*x2)
  term2 <- 20 * (x3-0.5)^2
  term3 <- 10*x4
  term4 <- 5*x5
  
  y <- term1 + term2 + term3 + term4
  return(y)
}

## 2. Explore the Model
library(ggplot2)
# Vary x1 whilst keeping all other inputs constant
x1_values <- seq(0,1,0.05)
length(x1_values)
x1_vary <- data.frame(x1_values, rep(0.5,21), rep(0.5,21), rep(0.5,21), rep(0.5,21))
plot(x1_values, t(fried(x1_vary)), xlab=expression('x'[1], 'values'), ylab='Output', type='l')


# Vary x2 whilst keeping all other inputs constant
x2_values <- seq(0,1,0.05)
length(x2_values)
x2_vary <- data.frame(rep(0.5,21), x2_values, rep(0.5,21), rep(0.5,21), rep(0.5,21))
plot(x2_values, t(fried(x2_vary)), xlab=expression('x'[2], 'values'), ylab='Output', type='l')

# Vary x3 whilst keeping all other inputs constant
x3_values <- seq(0,1,0.05)
length(x3_values)
x3_vary <- data.frame(rep(0.5,21), rep(0.5,21), x3_values, rep(0.5,21), rep(0.5,21))
plot(x3_values, t(fried(x3_vary)), xlab=expression('x'[3], 'values'), ylab='Output', type='l')

# Vary x4 whilst keeping all other inputs constant
x4_values <- seq(0,1,0.05)
length(x4_values)
x4_vary <- data.frame(rep(0.5,21), rep(0.5,21), rep(0.5,21), x4_values, rep(0.5,21))
plot(x4_values, t(fried(x4_vary)), xlab=expression('x'[4], 'values'), ylab='Output', type='l')

# Vary x5 whilst keeping all other inputs constant
x5_values <- seq(0,1,0.05)
length(x5_values)
x5_vary <- data.frame(rep(0.5,21), rep(0.5,21), rep(0.5,21), rep(0.5,21), x5_values)
plot(x5_values, t(fried(x5_vary)), xlab=expression('x'[5], 'values'), ylab='Output', type='l')


## 4. Construct a 1D Emulation

# Emulate with x3 varying and everything else fixed
fried_x3 <- function(x3)
{
  x1 <- 0.5
  x2 <- 0.5
  x3 <- x3
  x4 <- 0.5
  x5 <- 0.5
  
  term1 <- 10 * sin(pi*x1*x2)
  term2 <- 20 * (x3-0.5)^2
  term3 <- 10*x4
  term4 <- 5*x5
  
  y <- term1 + term2 + term3 + term4
  return(y)
}

# Define simple Bayes Linear emulator for single input ###
simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta,      # the correlation lengths
                                  sigma,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f         # prior expectation of f: E(f(x)) = 0 
){
  # store length of runs D  
  n <- length(D)
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2) 
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))}

# Evaluate the function over 6 evenly spaced runs
xD <- c(0,0.2,0.4,0.6,0.8,1)
D <- fried_x3(xD)

# Evaluate emulator over 201 prediction points xP 
xP <- seq(0.001,1,len=201)
nP <- length(xP)

em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=2,E_f=15))

plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ty="l",col="blue",lwd=2.5,
       xlab=expression("x"[3]),ylab="f(x)",main=maintitle, ylim=c(13,20)) # main argument: plot title
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,fried_x3(xP),lwd=2,lty=1)
  
  ### Plot the runs ###
  points(xD,D,pch=21,col=1,bg="green",cex=1.5)
  legend('bottomleft',legend=c("Emulator Expectation","Emulator Prediction Interval",
                                "True function f(x)","Model Evaluations"),
         lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
}


plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output")


## 5. Demonstrate 2D Emulation

# Define simple Bayes Linear emulator for single input in 2D ###
simple_BL_emulator_v2 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  # store length of runs D  
  n <- length(D)
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  # Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)    # XXX Old 1D version
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX New 2D version
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])  # XXX Old 1D version
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX New 2D version
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])    # XXX Old 1D version
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX New 2D version
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))}

# True function we are emulating by varying only x2 and x3
f <- function(x){
  x1 <- 0.5
  x2 <- x[,1]
  x3 <- x[,2]
  x4 <- 0.5
  x5 <- 0.5
  
  term1 <- 10 * sin(pi*x1*x2)
  term2 <- 20 * (x3-0.5)^2
  term3 <- 10*x4
  term4 <- 5*x5
  
  y <- term1 + term2 + term3 + term4
  return(y)}

### Define run locations ###
D_grid <- c(0.05,0.35,0.65,0.95)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid))
D <- f(xD)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))

simple_BL_emulator_v2(x=c(0.2,0.7),xD=xD,D=D,theta=0.5,sigma=3,E_f=15)

em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.5,sigma=3,E_f=15)) 
head(em_out)

E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

filled.contour(x_grid,x_grid,E_D_fx_mat,xlab="x1",ylab="x2")

cont_levs <- seq(6,22,1)      # define contour levels for both functions
filled.contour(x_grid,x_grid,E_D_fx_mat,xlab="x1",ylab="x2",levels=cont_levs,
               plot.axes={
                 axis(1);axis(2) 
                 contour(x_grid,x_grid,E_D_fx_mat,levels=cont_levs,add=TRUE,lwd=0.8)
                 points(xD,pch=21,col=1,bg="green",cex=1.5)
               }
)

# define filled contour plot function for emulator output ###
emul_fill_cont <- function(
    cont_mat,            # matrix of values we want contour plot of 
    cont_levs=NULL,      # contour levels (NULL: automatic selection)
    nlev=20,             # approx no. of contour levels for auto select  
    plot_xD=TRUE,        # plot the design runs TRUE or FALSE
    xD=NULL,             # the design points if needed
    xD_col="green",      # colour of design runs
    x_grid,              # grid edge locations that define xP
    ...                  # extra arguments passed to filled.contour
){
  
  ### Define contour levels if necessary ###
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)     
  
  ### create the filled contour plot ###
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab=expression("x"[2]),ylab=expression("x"[3]),...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8)   # plot contour lines
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(6,22,1),xD=xD,x_grid=x_grid,
               main="Emulator Adjusted Expectation E_D[f(x)]")

library(viridisLite)

emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(6,22,1),xD=xD,x_grid=x_grid,
               color.palette=magma,        # this sets the colour scheme
               main="Emulator Adjusted Expectation")

var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,
               main="Emulator Adjusted Variance")

fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(6,22,1),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,
               main="Friedman Function")

S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",
               color.palette=diag_cols,
               main="Emulator Diagnostics")


####### Latin Hypercube Design
# Create Type 1 LHD ###
set.seed(10)
nl <- 16
x_lhd <- cbind("x1"=sample(0:(nl-1)),"x2"=sample(0:(nl-1))) / nl  +  0.5/nl

### Maximin loop: performs swaps on 1st of two closest points with another random point
for(i in 1:10000){
  mat <- as.matrix(dist(x_lhd)) + diag(10,nl)   # creates matrix of distances between points 
  # note the inflated diagonal 
  closest_runs <- which(mat==min(mat),arr.ind=TRUE) # finds pairs of closest runs
  ind <- closest_runs[sample(nrow(closest_runs),1),1] # chooses one of close runs at random
  swap_ind <- sample(setdiff(1:nl,ind),1)       # randomly selects another run to swap with
  x_lhd2 <- x_lhd                               # creates second version of LHD
  x_lhd2[ind[1],1]   <- x_lhd[swap_ind,1] # swaps x_1 values between 1st close run and other run
  x_lhd2[swap_ind,1] <- x_lhd[ind[1],1]   # swaps x_1 values between 1st close run and other run
  if(min(dist(x_lhd2)) >= min(dist(x_lhd))-0.00001) {  # if min distance between points is same or better
    x_lhd <- x_lhd2                                      # we replace LHD with new LHD with the swap
    cat("min dist =",min(dist(x_lhd)),"Iteration = ",i,"\n") # write out min dist between points
  }
}

### plot maximin LHD ###
plot(x_lhd,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="blue",xlab=expression("x"[2]),ylab=expression("x"[3]),cex=1.4)
abline(h=(0:nl)/nl,col="grey60")
abline(v=(0:nl)/nl,col="grey60")


### Define run locations as the Maximin LHD found above ###
xD <- x_lhd 

### Perform 16 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.5,sigma=4,E_f=15))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Evaluate true function and store in matrix for diagnostic comparisons ###
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

### Evaluate diagnostics S_D(x) and store in matrix ###
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)

emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(6,22,1),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation")

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(6,22,1),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Friedman Function")

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance")

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",color.palette=diag_cols,main="Emulator Diagnostics")


################ History Matching ###################
lhd_maximin <- function(nl=16){                    # nl = number of points in LHD 
  
  x_lhd <- cbind("x1"=sample(0:(nl-1)),"x2"=sample(0:(nl-1))) / nl  +  0.5/nl  # create LHD
  
  ### Maximin loop: performs swaps on 1st of two closest points with another random point
  for(i in 1:1000){
    mat <- as.matrix(dist(x_lhd)) + diag(10,nl) # creates matrix of distances between points
    # note the inflated diagonal 
    closest_runs <- which(mat==min(mat),arr.ind=TRUE)   # finds pairs of closest runs
    ind <- closest_runs[sample(nrow(closest_runs),1),1] # chooses one of close runs at random
    swap_ind <- sample(setdiff(1:nl,ind),1)       # randomly selects another run to swap with
    x_lhd2 <- x_lhd                               # creates second version of LHD
    x_lhd2[ind[1],1]   <- x_lhd[swap_ind,1] # swaps x_1 vals between 1st close run & other run
    x_lhd2[swap_ind,1] <- x_lhd[ind[1],1]   # swaps x_1 vals between 1st close run & other run
    if(min(dist(x_lhd2)) >= min(dist(x_lhd))-0.00001) {  # if min distance between points is same or better
      x_lhd <- x_lhd2                                    # we replace LHD with new LHD with the swap
      # cat("min dist =",min(dist(x_lhd)),"Iteration = ",i,"\n") # write out min dist 
    }
  }
  
  ### plot maximin LHD ###
  plot(x_lhd,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="blue",
       xlab="x1",ylab="x2",cex=1.4)
  abline(h=(0:nl)/nl,col="grey60")
  abline(v=(0:nl)/nl,col="grey60")
  return(x_lhd)
}

set.seed(1)
xD_w1 <- lhd_maximin(nl=10)

x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))

### Defined Extra Objects for HM and Implausibility ###
z <- 11.2
sigma_e <- 0.1
sigma_epsilon <- 0.75

### Define current run locations ###
xD <- xD_w1

### Perform 14 runs of model and store as D (this would take days for realistic example!) ###
D <- f(xD)

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.5,sigma=4,E_f=12))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 

### Define colours and levels for implausibility plots ###
imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),25)
### define filled contour plot function for emulator output ###
emul_fill_cont_V2 <- function(
    cont_mat,            # matrix of values we want contour plot of 
    cont_levs=NULL,      # contour levels (NULL: automatic selection)
    cont_levs_lines=NULL,   # contour levels for lines (NULL: automatic selection)
    nlev=20,             # approx no. of contour levels for auto select  
    plot_xD=TRUE,        # plot the design runs TRUE or FALSE
    xD=NULL,             # the design points if needed
    xD_col="green",      # colour of design runs
    x_grid,              # grid edge locations that define xP
    ...                  # extra arguments passed to filled.contour
){
  
  ### Define contour levels if necessary ###
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)
  
  ### create the filled contour plot ###
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab=expression("x"[2]),ylab=expression("x"[3]),...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   if(is.null(cont_levs_lines)) contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8) # plot usual contour lines 
                   if(!is.null(cont_levs_lines)) {
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.4,labels="")   # plot thin contour lines 
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs_lines,lwd=2)   # plot thick contour lines
                   }
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}
### plot wave 1 implausibility and wave 1 runs only ###
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3*sqrt(sigma_e^2 + sigma_epsilon^2),xD=xD,x_grid=x_grid,
                  xD_col="purple",color.palette=imp_cols,main="Implausibility: Wave 1")
### plot wave 1 emulator expectation ###
emul_fill_cont_V2(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                  color.palette=exp_cols,main="Emulator Adjusted Expectation")
### plot wave 1 emulator variance ###
emul_fill_cont_V2(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                  color.palette=var_cols,main="Emulator Adjusted Variance")


xD_w2 <- matrix(     # the 6 point wave 2 design (chosen by hand)
  c(0.98,0.02,
    0.03,0.95,
    0.26,0.4,
    0.1,0.8,
    0.13,0.65,
    0.07,0.1,
    0.2,0.24),ncol=2,byrow=TRUE
)
### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3*sqrt(sigma_e^2 + sigma_epsilon^2),xD=rbind(xD,xD_w2),
                  x_grid=x_grid,xD_col=rep(c("purple","pink"),c(nrow(xD),nrow(xD_w2))),
                  color.palette=imp_cols,main="Implausibility: Wave 1")

### loop over adding the wave 2 runs: add k runs ###
for(k in 0:8){                          # k=0: wave 1, k>0 add k wave 2 runs sequentially
  
  xD <- xD_w1                           # the 14 point wave 1 design
  if(k>0) xD <- rbind(xD,xD_w2[1:k,])   # k=0: wave 1, k>0 add wave 2 runs sequentially
  
  ### Perform 14 + k runs of model and store as D (would take days for realistic example!) ###
  D <- f(xD)
  
  ### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.5,sigma=4,E_f=12))   
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Evaluate true function and store in matrix for diagnostic comparisons ###
  fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
  Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  ### Calculate Imp Measure for True f(x) Over All 50x50 = 2500 input points in xP ###
  Imp_true_mat <- sqrt( (fxP_mat - z)^2 / (sigma_e^2 + sigma_epsilon^2) ) 
  
  ### Define colours and levels for implausibility plots ###
  imp_cols <- function(n) turbo(n,begin=0.15,end=1)
  imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),25)
  
  ### if k=0 plot wave 1 runs only ###
  if(k==0) emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3*sqrt(sigma_e^2 + sigma_epsilon^2),xD=xD,
                             x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                             main="Implausibility: Wave 1")
  ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
  emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3*sqrt(sigma_e^2 + sigma_epsilon^2),xD=rbind(xD,xD_w2),
                    x_grid=x_grid,xD_col=rep(c("purple","pink"),c(nrow(xD)+k,nrow(xD_w2)-k)),  # cover unevaluated w2 points in pink
                    color.palette=imp_cols,main="Implausibility: Wave 2")
  ### once last run done so k=8, plot implausibility for true function f(x) to compare ###
  if(k==nrow(xD_w2)) emul_fill_cont_V2(cont_mat=Imp_true_mat,cont_levs=imp_levs,
                                       cont_levs_lines=3*sqrt(sigma_e^2 + sigma_epsilon^2),xD=xD,x_grid=x_grid,xD_col="purple",plot_xD=FALSE,
                                       color.palette=imp_cols,main="Implausibility using Friedman Function")
}
