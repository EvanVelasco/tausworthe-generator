library(ggplot2)
library(plotly)
library(dplyr)


# Functions

# helper function for tausworthe generator
bitToBase10 <- function(x) {
  bits <- as.numeric(unlist(strsplit(as.character(x),"")))
  expo = rev(seq(1,length(bits),1)-1)
  result = sum(bits*2^expo)
  return (result)
}


# primary tausworthe PRN function
tausworthe <- function(seed,n,r,q,l){
  
  # Unpack seed into initial values for B
  B = as.numeric(unlist(strsplit(as.character(seed), "")))
  
  # Determine number of bits required to genera te n random variables
  n = l*n
  
  # Add additional bits
  for (i in (q+1):n){
    B[i] = xor(B[i - r], B[i - q]) 
  }
  
  # Create sequence of indices to split bits on
  s = seq(from = 1, to = n, by = l)
  
  results = c()
  for (i in s){
    # Grab relevant group of bits
    values = B[i:(i+l-1)]
    # Collapse
    values = (paste0(values, collapse=""))
    # Use bitToBase10 to calculate values
    results = c(results, bitToBase10(values) / 2^l)
  }
  return(results)
}


# function to create dataframe for adjacent PRNs in 2 and 3 dimensions
adjacentPrns <- function(x, dimensions=2) {
  if (dimensions == 2){
    df <- data.frame(u = x[1:(length(x) - 1)],
                     u1 = x[2:(length(x))])
  } else if (dimensions == 3){
    df <- data.frame(u = x[1:(length(x)-2)],
                     u1 = x[2:(length(x)-1)],
                     u2 = x[3:length(x)])
  }
  return(df)
}


# chi square gof test for uniform
uniformChiSq <- function(x, alpha=0.05, k=5){
  
  # count obs in each fold
  folds = seq(0, 1, 1/k)
  obs_results = numeric(k)
  for (i in 1:k){
    obs_results[i] = sum(folds[i] <= x & x < folds[i+1])
  }
  
  # edge case. if any obs == 1, add those obs to the final fold
  if (sum(obs_results)!= length(x)) {
    obs_results[k] = obs_results[k] + (length(x) - sum(obs_results))
  }
  
  # create vector of expected probabilities
  exp_results = rep(1/k, k)
  
  # compute test statistic and p value
  X0 = chisq.test(x = obs_results, p = exp_results)$statistic
  Xa = qchisq(1-alpha, df=4)
  pvalue = 1-pchisq(X0, 4)
  
  # determine outcome
  if (X0 > Xa){ outcome = "reject"}
  else {outcome = "fail to reject"}
  
  # print results
  print("Chi Sq Test for Goodness of Fit for Uniforms with k=5")
  print("  H0: PRNs are uniform(0,1)")
  print("  HA: PRNs are NOT uniform(0,1)")
  print("---------------------------------------------------------------------------")
  print(paste0("We ", outcome, " the null hypothesis at a=", alpha, " level of significance."))
  print(paste0("X0 = ", X0, ", Xa = ", round(Xa,4), ", p-value = ", round(pvalue,4)))
}

# Runs Test
runsUpDown <- function(prn, alpha=0.05){
  n = length(prn)
  
  # define sequence of up and down
  runs = numeric(length(prn)-1)
  for (i in 1:(length(prn)-1)){
    if(prn[i] < prn[i+1]){
      runs[i] = 1
    } else{
      runs[i] = 0
    }
  }
  
  # count the number of ups and downs
  A = 1
  for (i in 1:(length(prn)-2)){
    if (runs[i]!=runs[i+1]){
      A = A + 1
    }
  }
  
  # compute test statistics and p value
  A_mean = (2*n-1)/3
  A_var = (16*n-29)/90
  
  z0 = (A - A_mean)/sqrt(A_var)
  
  za = qnorm(1-alpha/2)
  
  pvalue = 2*(1-pnorm(abs(z0)))
  
  # determine outcome
  if (abs(z0) > za){ outcome = "reject"}
  else {outcome = "fail to reject"}
  
  print("Runs Test Up and Down for Independence of Uniforms")
  print("  H0: PRNs are independent of one another")
  print("  HA: PRNs are NOT independent of one another")
  print("---------------------------------------------------------------------------")
  print(paste0("We ", outcome, " the null hypothesis at a=", alpha, " level of significance."))
  print(paste0("runs = ", A, ", |z0| = ", round(abs(z0),4), ", za = ", round(za,4),
               ", p-value = ", round(pvalue,4)))
}

runsAboveBelow <- function(prn, alpha=0.05){
  n1 = 0
  n2 = 0
  n = length(prn)
  
  # determine number of obs above and below the mean
  for (i in 1:length(prn)){
    if (prn[i] >= 0.5){
      n1 = n1 + 1
    } else if (prn[i] < 0.5){
      n2 = n2 + 1
    }
  }
  
  # calculate runs above and below the mean
  runs = numeric(length(prn))
  for (i in 1:(length(prn)-1)){
    if(prn[i] >= 0.5){
      runs[i] = 1
    } else{
      runs[i] = 0
    }
  }
  
  # count the number of runs
  B = 1
  for (i in 1:(length(prn)-1)){
    if (runs[i]!=runs[i+1]){
      B = B + 1
    }
  }
  
  # calculate test statistics and p value
  B_mean = ((2*n1*n2)/n) + 1/2
  B_var = (2*n1*n2*(2*n1*n2 - n)) / ((n^2)*(n-1))
  
  z0 = (B - B_mean)/sqrt(B_var)
  
  za = qnorm(1-alpha/2)
  
  pvalue = 2*(1-pnorm(abs(z0)))
  
  # determine outcome
  if (abs(z0) > za){ outcome = "reject"}
  else {outcome = "fail to reject"}
  
  print("Runs Test Above and Below the Mean for Independence of Uniforms")
  print("  H0: PRNs are independent of one another")
  print("  HA: PRNs are NOT independent of one another")
  print("---------------------------------------------------------------------------")
  print(paste0("We ", outcome, " the null hypothesis at a=", alpha, " level of significance."))
  print(paste0("runs = ", B, ", |z0| = ", round(abs(z0),4), ", za = ", round(za,4),
               ", p-value = ", round(pvalue,4)))
}

# Correlation Test
correlationTest <- function(prn, alpha=0.05){
  n = length(prn)
  
  # calculate the sum of the product of adjacent prns
  sums = 0
  for (i in 1:(n-1)){
    sums = sums + prn[i]*prn[i+1]
  }
  
  # calculate test statistics and p value
  rho = ((12/(n-1))*sums)-3
  rho_mean = 0
  rho_var = (13*n-19) / (n-1)^2
  
  z0 = rho/sqrt(rho_var)
  za = qnorm(1-alpha/2)
  
  pvalue = 2*(1-pnorm(abs(z0)))
  
  # determine results
  if (abs(z0) > za){ outcome = "reject"}
  else {outcome = "fail to reject"}
  
  print("Correlation Test for Independence of Uniforms")
  print("  H0: PRNs are independent of one another")
  print("  HA: PRNs are NOT independent of one another")
  print("---------------------------------------------------------------------------")
  print(paste0("We ", outcome, " the null hypothesis at a=", alpha, " level of significance."))
  print(paste0("|z0| = ", round(abs(z0),4), ", za = ", round(za,4),
               ", p-value = ", round(pvalue,4)))
}

# Normal Variate Generation
boxMullerNorm <- function(u1, u2){
  Z = sqrt(-2*log(u1))*cos(2*pi*u2)
  return(Z)
}


