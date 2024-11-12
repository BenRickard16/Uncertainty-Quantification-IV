## 2. Creating an Emulator

# Defining the actual computer model / simulator
f <- function(x) sin(2*pi*x)

# Define the run locations
xD <- seq(0, 1, 0.2)
xD

# Perform 6 runs of model f and store as D
D <- f(xD)
D

n <- length(D)
n

# Define the Gaussian covariance structure of the emulator
theta <- 0.25     # correlation length
sigma <- 0.5      # prior standard deviation

Cov_fx_fxdash <- function(x, xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)     # covariance structure

E_f = 0         # prior expectation

# Now we define the 5 objects required to perform the Bayes Linear adjustment
# Expectation of the runs D
E_D <- rep(E_f, n)
E_D

# Variance matrix for runs D
Var_D <- matrix(0, nrow=n, ncol=n)
for(i in 1:n) {
  for(j in 1:n) {
    Var_D[i, j] <- Cov_fx_fxdash(xD[i], xD[j])
  }
}
Var_D

# Expectation of f(x)
E_fx <- E_f
E_fx

# Variance of f(x)
Var_fx <- sigma^2
Var_fx

# Decide on the x value to evaluate the emulator at, then construct Cov[f(x), D] row vector
x <- 0.25

Cov_fx_D <- matrix(0, nrow=1, ncol=n)
for(j in 1:n){
  Cov_fx_D[1,j] <- Cov_fx_fxdash(x, xD[j])
}

# Now we can perform the Bayes Linear adjustment
# The BL adjusted expectation
ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)

# The BL adjusted variance
VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)



## 3. Emulating Many Points: Creating an Emulator Function

# Simple Bayes Linear emulator for a single input x
simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
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
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}

# Check that we get the same resutls as before
# Same 6 runs
xD <- seq(0, 1, 0.2)
D <- f(xD)

# Using the function at x=0.25
em_out_1pt <- simple_BL_emulator_v1(x=0.25, xD=xD, D=D, theta=0.25, sigma=0.5, E_f=0)
em_out_1pt

# Constructing a 3-sigma prediction interval
low1 <- em_out_1pt["ExpD_f(x)"] - 3*sqrt(em_out_1pt["VarD_f(x)"])
up1 <- em_out_1pt["ExpD_f(x)"] + 3*sqrt(em_out_1pt["VarD_f(x)"])
c(low1, up1, use.names=FALSE)     # Note true answer is in this interval

# Evaluating emulator over 201 prediction points
# New prediction points
xP <- seq(0.001, 0.999, len=201)

# Number of prediction points
nP <- length(xP)

# Set up matrix for results
em_out <- matrix(0, nrow=nP, ncol=2, dimnames=list(NULL, c("ExpD_f(x)","VarD_f(x)")))

# Filling in with adjusted expectations and variances
for(i in 1:nP) em_out[i,] <- simple_BL_emulator_v1(x=xP[i],xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0)
head(em_out)  
  
 # Can do the same thing with sapply() function
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
head(em_out)  
  
  

## 4. Plotting the Emulator Output

# Plot the emulator output, uncommenting top and bottom lines will generate a pdf
# pdf(file="fig1A_sin_emulator.pdf",height=7,width=8.5)    # uncomment if you want a pdf file
plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.2,1.2),ty="l",col="blue",lwd=2.5,xlab="Input parameter x",ylab="Output f(x)")
lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)

# Plot true function: we would not normally be able to do this! ###
lines(xP,f(xP),lwd=2,lty=1)

# Plot expectation and prediction interval points at x=0.25 ###
abline(v=0.25,lty=2)
points(0.25,em_out_1pt["ExpD_f(x)"],col=1,bg="blue",pch=21)
points(c(0.25,0.25),c(low1,up1),col=1,bg="red",pch=21)

# Plot the runs 
points(xD,D,pch=21,col=1,bg="green",cex=1.5)

# Add a legend to label everything 
legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                           "True function f(x)","Model Evaluations","x = 0.25 line"),
       lty=c(1,1,1,NA,2),pch=c(NA,NA,NA,16,NA),col=c("blue","red",1,"green",1),
       lwd=c(2.5,2.5,2.5,2.5,1),pt.cex=1.3, cex=0.7)
# dev.off()       # uncomment if you want to make a pdf file

# Function to plot simple emulator output
plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))
    maintitle=NULL, # title of the plot, blank if NULL
    cex=0.7
){

  # Plot emulator output 
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.2,1.2),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) 
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  # Plot true function: we would not normally be able to do this!
  lines(xP,f(xP),lwd=2,lty=1)
  
  ### Plot the runs 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5)
  legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                             "True function f(x)","Model Evaluations"),
         lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3, cex=0.7)
}

par(mfrow=c(1,1))
# Checking it works as expected
plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output")



## 5. Exploring Emulator Behaviour

# Add 2 more runs
# Define run locations 
xD <- c(seq(0,1,0.2),0.7,0.9)    # original 6 runs and two extra runs at x=0.7 and x=0.9

# Perform 8 runs of model and store as D 
D <- f(xD)

# Evaluate emulator over 201 prediction points xP 
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  

# Plot emulator output 
plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output: 8 runs")


# Back to 6 runs
# Define run locations 
xD <- seq(0,1,0.2)    # original 6 runs

# Perform 6 runs of model and store as D 
D <- f(xD)

# Testing the effect of varying sigma
# Make sequence of sigma values to use for emulation 
sigma_seq <- c(0.25,0.5,0.75,1,1.5)

# For loop over different sigma values in sigma_seq 
for(i in 1:length(sigma_seq)){
  # Evaluate emulator over 201 prediction points xP with sigma=sigma_seq[i] 
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=sigma_seq[i],E_f=0))
  
  # Plot emulator output in each case, note use of "paste" for plot title 
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Sigma =",sigma_seq[i]))
}


## 5.1 Investigate the affect of varying correlation length parameter theta
# Sequence of theta values to use for emulation
theta_seq <- c(0.0005, 0.05, 0.2, 0.5, 1)

# For loop over different theta values in sigma_seq 
for(i in 1:length(theta_seq)){
  # Evaluate emulator over 201 prediction points xP with sigma=sigma_seq[i] 
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_seq[i],sigma=0.5,E_f=0))
  
  # Plot emulator output in each case, note use of "paste" for plot title 
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Theta =",theta_seq[i]))
}
 

## 5.2 Explore different locations for the 6 runs
# Define a function based off the sequence of initial run locations
location_variation <- function(xD){
  # Perform 6 runs of model and store as D 
  D <- f(xD)
  
  # Evaluate emulator over 201 prediction points xP 
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
  
  # Plot emulator output 
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output: 8 runs")
  
}

# Try the emulator for different initials runs
location_variation(c(0.4,0.1,0.2,0.8,0.9,0.6))
