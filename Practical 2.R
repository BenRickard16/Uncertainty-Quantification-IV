## 2. Creating an Emulator for 2-Dimensional Inputs

# Define simple Bayes Linear emulator for single input in 2D 

simple_BL_emulator_v2 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # store length of runs D  
  n <- length(D)
  
  # Define Covariance structure of f(x): Cov[f(x),f(xdash)]
  # Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)    # XXX Old 1D version
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX New 2D version
  
  
  # Define 5 objects needed for BL adjustment
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
  
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x)
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  # Return emulator expectation and variance 
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}



## 3. Applying the 2-Dimensional Emulator

# Define actual 2D computer model/simulator
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))

# Define run locations
D_grid <- c(0.05,0.35,0.65,0.95)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid))
xD

# Perform 16 runs of model and store as D
D <- f(xD)
D

# Define 50x50 grid of prediction points xP for emulator evaluation 
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
dim(xP)     # Check the dimensions of xP

# Predict at single point x = (0.2, 0.7)
simple_BL_emulator_v2(x=c(0.2,0.7),xD=xD,D=D,theta=0.45,sigma=1,E_f=0)

# Evaluate emulator over 50x50=2500 prediction points xP
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.45,sigma=1,E_f=0))   
head(em_out)

# Store emulator output as matrices to aid plotting
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 



## 4. Plotting 2D Emulator Output

# Want to use 2D coloured contour plots
filled.contour(x_grid,x_grid,E_D_fx_mat,xlab="x1",ylab="x2")

# Adding the initial runs
filled.contour(x_grid,x_grid,E_D_fx_mat,xlab="x1",ylab="x2",
               plot.axes={
                 axis(1);axis(2) 
                 points(xD,pch=21,col=1,bg="green",cex=1.5)})

# Adding solid contour lines with numbers
cont_levs <- seq(-2,2,0.2)      # define contour levels for both functions
filled.contour(x_grid,x_grid,E_D_fx_mat,xlab="x1",ylab="x2",levels=cont_levs,
               plot.axes={
                 axis(1);axis(2) 
                 contour(x_grid,x_grid,E_D_fx_mat,levels=cont_levs,add=TRUE,lwd=0.8)
                 points(xD,pch=21,col=1,bg="green",cex=1.5)})

# Create a plotting function to do all of this for us
# Define filled contour plot function for emulator output
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
  
  # Define contour levels if necessary 
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)     
  
  # Create the filled contour plot 
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab="x1",ylab="x2",...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8)   # plot contour lines
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}

# The ... means can add extra arguments
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               main="Emulator Adjusted Expectation E_D[f(x)]")



## 5. Colour Schemes for 2D Contour Plots

# Will use library viridisLite and function hcl.colors()
library(viridisLite)

exp_cols <- magma
var_cols <- function(n) hcl.colors(n, 'YlOrRd', rev=TRUE)
diag_cols <- turbo

# Trying some colour schemes out
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,        # this sets the colour scheme
               main="Emulator Adjusted Expectation E_D[f(x)]")

# Now for emulator variance
emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,
               main="Emulator Adjusted Variance Var_D[f(x)]")

# Evaluate true function and store in matrix for diagnostic comparisons to compare emulator against
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,
               main="True Computer Model Function f(x)")

# Evaluate diagnostics S_D(x) and store in matrix
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",
               color.palette=diag_cols,
               main="Emulator Diagnostics S_D[f(x)]")



## 6. Changing Correlation Lengths

# Investigating impact of correlation length on emulator variance and expectation
# Define a vector of theta values to use 
theta_seq <- c(0.05,0.1,0.15,0.2,0.35,0.4,0.45)

# Loop over vector of theta values 
for(i in 1:length(theta_seq)){
  
  # Evaluate emulator over 201 prediction points xP and store in matrices 
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  # Plot filled contour plot of emulator expectation 
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
                 color.palette=magma,main=paste("Emul. Adjusted Expectation E_D[f(x)], theta =",theta_seq[i]))
}

# Now checking emulator variance
# Define a vector of theta values to use 
theta_seq <- c(0.05,0.1,0.15,0.2,0.35,0.4,0.45)

# Loop over vector of theta values 
for(i in 1:length(theta_seq)){
  
  # Evaluate emulator over 201 prediction points xP and store in matrices 
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  # Plot filled contour plot of emulator expectation 
  emul_fill_cont(cont_mat=Var_D_fx_mat,nlev=12,xD=xD,x_grid=x_grid,color.palette=var_cols,
                 main=paste("Emul. Adjusted Variance Var_D[f(x)], theta =",theta_seq[i]))
}



## 7. Investigating Grid Designs

# Define run locations
D_grid <- c(0.2,0.5,0.8)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid))

# Perform 9 runs of model and store as D (this would takes days for realistic example!)
D <- f(xD)

# Evaluate emulator over 50x50=2500 prediction points xP and store as matrices
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.45,sigma=1,E_f=0))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate true function and store in matrix for diagnostic comparisons 
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate diagnostics S_D(x) and store in matrix 
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)


# Plot emulator expectation, variance, true function, and diagnostics
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="True Computer Model Function f(x)")

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",color.palette=diag_cols,main="Emulator Diagnostics S_D[f(x)]")


## Exercise 7.1 - redo code with slightly larger grid based on
D_grid <- c(0.1, 0.5, 0.9)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid))

# Perform 9 runs of model and store as D (this would takes days for realistic example!)
D <- f(xD)

# Evaluate emulator over 50x50=2500 prediction points xP and store as matrices
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.45,sigma=1,E_f=0))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate true function and store in matrix for diagnostic comparisons 
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate diagnostics S_D(x) and store in matrix 
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)


# Plot emulator expectation, variance, true function, and diagnostics
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="True Computer Model Function f(x)")

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",color.palette=diag_cols,main="Emulator Diagnostics S_D[f(x)]")

## Exercise 7.1 - create a 5x5 grid design
D_grid <- seq(0, 1, 0.25)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid))

# Perform 9 runs of model and store as D (this would takes days for realistic example!)
D <- f(xD)

# Evaluate emulator over 50x50=2500 prediction points xP and store as matrices
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.45,sigma=1,E_f=0))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate true function and store in matrix for diagnostic comparisons 
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate diagnostics S_D(x) and store in matrix 
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)


# Plot emulator expectation, variance, true function, and diagnostics
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="True Computer Model Function f(x)")

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",color.palette=diag_cols,main="Emulator Diagnostics S_D[f(x)]")



## 8. Benefit of Latin Hypercube Designs

# Define actual 2D computer model/simulator 
f <- function(x) sin(3*pi*x[,1]) + (1/20)*cos(2*pi*x[,2]) 

# Define run locations 
D_grid <- c(0,1/3,2/3,1) 
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid))

# Perform 16 runs of model and store as D (this would takes days for realistic example!) 
D <- f(xD)

# Evaluate emulator over 50x50=2500 prediction points xP and store as matrices 
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.35,sigma=1,E_f=0))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate true function and store in matrix for diagnostic comparisons 
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate diagnostics S_D(x) and store in matrix 
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)


# Now plot the emulator expectation, true function, emulator variance, and diagnostics
# Emulator expectation
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

# True function
emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="True Computer Model Function f(x)")

# Emulator variance
emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")

# Diagnostics
emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",color.palette=diag_cols,main="Emulator Diagnostics S_D[f(x)]")

# No we try a LHD using the algorithm from lectures
# Set number of points in LHD 
nl <- 16        # size of LHC

# For each input, sample a permutation from 0:(nl-1)/nl, and cbind together 
x_temp <- cbind("x1"=sample(0:(nl-1)),"x2"=sample(0:(nl-1))) / nl 
x_temp

# Plot with grid to see points are in bottom left of each box 
plot(x_temp,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="red",xlab="x1",ylab="x2",cex=1.4)
abline(h=(0:nl)/nl,col="grey60")
abline(v=(0:nl)/nl,col="grey60")

# For Type 1 we add 1/(2*nl) to recenter in boxes 
x_lhd <- x_temp + 0.5/nl
x_lhd

# Plot with grid to see points are now in centre of each box 
plot(x_lhd,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="red",xlab="x1",ylab="x2",cex=1.4)
abline(h=(0:nl)/nl,col="grey60")
abline(v=(0:nl)/nl,col="grey60")

# We now attempt to improve the LHD by moving points close together further apart to get us closer
# to LHD with maximum minimum distance between points

# Create Type 1 LHD
x_lhd <- cbind("x1"=sample(0:(nl-1)),"x2"=sample(0:(nl-1))) / nl  +  0.5/nl

# Maximin loop: performs swaps on 1st of two closest points with another random point
for(i in 1:1000){
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

# Plot maximin LHD 
plot(x_lhd,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="blue",xlab="x1",ylab="x2",cex=1.4)
abline(h=(0:nl)/nl,col="grey60")
abline(v=(0:nl)/nl,col="grey60")

# Now we use new Maximin LHD
# Define run locations as the Maximin LHD found above 
xD <- x_lhd 

# Perform 16 runs of model and store as D (this would takes days for realistic example!) 
D <- f(xD)

# Evaluate emulator over 50x50=2500 prediction points xP and store as matrices 
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.35,sigma=1,E_f=0))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate true function and store in matrix for diagnostic comparisons 
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

# Evaluate diagnostics S_D(x) and store in matrix 
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)

# Plot emulator expectation, variance, true function and diagnostics 
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="True Computer Model Function f(x)")

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",color.palette=diag_cols,main="Emulator Diagnostics S_D[f(x)]")

