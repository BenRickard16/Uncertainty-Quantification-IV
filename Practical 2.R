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





