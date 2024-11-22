## 2. 1D History Matching: Wave 1


# Define simple Bayes Linear emulator for single input
simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # Store length of runs D  
  n <- length(D)
  
  # Define Covariance structure of f(x): Cov[f(x),f(xdash)] 
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2) 
  
  
  # Define 5 objects needed for BL adjustment 
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
  
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) 
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  # Return emulator expectation and variance 
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
}


# Function to plot simple emulator output
plot_BL_emulator_V2 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))
    maintitle=NULL,  # title of the plot, blank if NULL
    z=NULL,         # the observed data z
    sigma_e=NULL,   # the observation errors SD s.t. Var[e] = sigma_e^2
    sigma_epsilon=NULL, # the model discrepancy SD s.t. Var[epsilon] = sigma_epsilon^2
    plot_true=FALSE   # don't plot true function unless this is TRUE
){
  
  # Plot emulator output 
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.3,1.2),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  # Plot true function: we would not normally be able to do this! 
  if(plot_true) lines(xP,f(xP),lwd=2,lty=1)
  
  # Plot the runs 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5)
  
  # Plot the Observed data plus errors due to obs and MD 
  if(!is.null(z)){
    abline(h=z,lwd=1.4)
    abline(h=z+3*sqrt(sigma_e^2+ sigma_epsilon^2),lty=2,lwd=1.2)
    abline(h=z-3*sqrt(sigma_e^2+ sigma_epsilon^2),lty=2,lwd=1.2)
  }
  
  if(plot_true) legend('topright',legend=c("Emulator Expectation",
                                           "Emulator Prediction Interval",
                                           "True function f(x)","Model Evaluations"), 
                       lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
  if(!is.null(z)) legend('topright',legend=c("Emulator Expectation",
                                             "Emulator Prediction Interval",
                                             "Model Evaluations","Observation z",
                                             "3 sigma interval"), 
                         lty=c(1,1,NA,1,2),pch=c(NA,NA,16,NA,NA),col=c("blue","red","green",1,1),
                         lwd=c(2.5,2.5,NA,1.4,1.2),pt.cex=1.3)
}

# Loading colour schemes
library(viridisLite)

# 1D HM example wave 1
# Define actual computer model/simulator
f <- function(x) sin(2*pi*x)

# Evaluate emulator over 201 prediction points xP
xP <- seq(0.001, 0.999, len=201)

# Define run locations
xD <- c(0, 0.2, 0.4, 0.6, 0.86, 1)    # Shifted x^(5) location

# Perform 6 runs of model and store as D (would realistically take days)
D <- f(xD)
