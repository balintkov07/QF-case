#Read excel imports
install.packages("readxl")
library("readxl")

#For skewness & kurtosis analysis
install.packages("GGally")
library("GGally")

#--------------------------DATA PROCESSING --------------------------


#Extract Data from Excel, confirm correct columns, exclude first observation since NA
df_init <- read_excel("/Users/janekczajnik/Desktop/Erasmus/B3/Introductory Seminar CS/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

#Divide the days into positive, negative or zero returns 
df_pos <- df_init[df_init$`CC Return (%)` > 0, ]
df_neg <- df_init[df_init$`CC Return (%)` < 0, ]
df_none <- df_init[df_init$`CC Return (%)` == 0, ]

#ensure correct filtering
head(df_pos)
head(df_neg)

#Summary in total number of days
print(cat(("| Number of Positive Days: "), nrow(df_pos), ("| Number of Negative Days: "), nrow(df_neg), ("| Number of Zero-return Days: "), nrow(df_none), "\n"))

#Arbitrary - baseline plots -> Returns, Realized Variance, VIX
plot(df_init$`CC Return (%)`)
plot(df_init$`RV5_SS × 10^4`)
plot(df_init$VIX)

#Baseline summary statistics of key variables
summary(df_init$`CC Return (%)`)
summary(df_init$`RV5_SS × 10^4`)
summary(df_init$VIX)
print(cat(("Skewness and Kurtosis (respectively) for Close to Close returns"), skewness(df_init$`CC Return (%)`), kurtosis(df_init$`CC Return (%)`)))
print(cat(("Skewness and Kurtosis (respectively) for Realized variance * 10^4"), skewness(df_init$`RV5_SS × 10^4`), kurtosis(df_init$`RV5_SS × 10^4`)))
print(cat(("Skewness and Kurtosis (respectively) for VIX index"), skewness(df_init$VIX), kurtosis(df_init$VIX)))

#--------------------------GARCH SETUP --------------------------
rt <- as.numeric(df_init$`CC Return (%)`)
mu <- mean(rt)
T <- length(rt)

garch11 <- function(par, rt, mu) {
  
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || alpha + beta >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- var(rt)
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  # Generate conditional variances recursively
  for (t in 1:(T - 1)) {
    h[t + 1] <- omega + alpha * (rt[t] - mu)^2 + beta * h[t]
  }
  
  # Now check h after it has been generated
  if (any(h <= 0) || any(is.na(h)) || any(is.infinite(h))) {
    return(1e10)
  }
  
  # Negative log-likelihood
  garch11ll <- 0.5 * sum(log(2 * pi) + log(h) + ((rt - mu)^2 / h))
  
  if (!is.finite(garch11ll)) {
    return(1e10)
  }
  
  return(garch11ll)
}

# Starting values
start_par <- c(
  omega = 0.05,
  alpha = 0.10,
  beta  = 0.85
)

# Estimate model
garch_fit <- optim(
  par = start_par,
  fn = garch11,
  rt = rt,
  mu = mu,
  method = "L-BFGS-B",
  #Omega > 0 to inf, remaining are bounded by 0 and 1 
  lower = c(1e-8, 0, 0),
  upper = c(Inf, 1, 1)
)

#post - optimisation parameter estimates 
garch_fit$par
garch_fit$value
garch_fit$convergence

h_long = (as.numeric(garch_fit$par[1]))/(as.numeric(1 - garch_fit$par[2] - garch_fit$par[3]))

#model implied unconditional variance 
print(h_long)

#compare to empirical unconditional variance
var(rt)



#GJR-GARCH
garch_gjr <- function(par, rt, mu) {
  
  omega  <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  beta   <- par[4]
  
  # Parameter restrictions
  # omega > 0, alpha1 >= 0, alpha2 >= 0, beta >= 0
  # Stationarity: alpha1 + 0.5 * alpha2 + beta < 1
  if (omega <= 0 || alpha1 < 0 || alpha2 < 0 || beta < 0 ||
      alpha1 + 0.5 * alpha2 + beta >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance
  h[1] <- var(rt)
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  shock <- rt - mu
  
  for (t in 1:(T - 1)) {
    
    indicator <- ifelse(shock[t] < 0, 1, 0)
    
    h[t + 1] <- omega +
      alpha1 * shock[t]^2 +
      alpha2 * indicator * shock[t]^2 +
      beta * h[t]
  }
  
  if (any(h <= 0) || any(is.na(h)) || any(is.infinite(h))) {
    return(1e10)
  }
  
  ll <- 0.5 * sum(log(2 * pi) + log(h) + shock^2 / h)
  
  if (!is.finite(ll)) {
    return(1e10)
  }
  
  return(ll)
}


start_par_gjr <- c(
  omega  = 0.05,
  alpha1 = 0.03,
  alpha2 = 0.07,
  beta   = 0.85
)

gjr_fit <- optim(
  par = start_par_gjr,
  fn = garch_gjr,
  rt = rt,
  mu = mu,
  method = "L-BFGS-B",
  lower = c(1e-8, 0, 0, 0),
  upper = c(10, 1, 1, 1)
)

gjr_fit$par
gjr_fit$value
gjr_fit$convergence

shock <- rt - mu

# Time-varying effective alpha
alpha1_hat <- gjr_fit$par[2]
alpha2_hat <- gjr_fit$par[3]
alpha_eff <- (alpha1_hat + alpha2_hat) / 2

h_long_gjr <- (as.numeric(gjr_fit$par[1]))/(as.numeric(1 - alpha_eff - gjr_fit$par[4]))
print(h_long_gjr)
