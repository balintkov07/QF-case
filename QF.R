#Read excel imports
install.packages("readxl")
library("readxl")

#For skewness & kurtosis analysis
install.packages("GGally")
library("GGally")

install.packages("e1071")
library(e1071)

#--------------------------DATA PROCESSING --------------------------

# Janek:
# df_init <- read_excel("/Users/janekczajnik/Desktop/Erasmus/B3/Introductory Seminar CS/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

# Balint
#df_init <- read_excel("/Users/balintkovacs/Documents/GitHub/QF-case/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

# Luca
setwd("C:/Users/lucam/Dropbox/dad&mum/University/Erasmus/BSC 3/BLK 5/Intro to Seminars/Case2_QF/econometrics & QF case study 2026/data")

#Extract Data from Excel, confirm correct columns, exclude first observation since NA

df_init <- read_excel("data.xlsx")[-1, ]

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
      (alpha1 + alpha2) / 2 + beta >= 1) {
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
    
    h[t+1] <- omega +
      alpha1 * indicator * shock[t]^2 +        # alpha1 when negative
      alpha2 * (1 - indicator) * shock[t]^2 +  # alpha2 when positive
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
  alpha1 = 0.15,  # negative shock reaction — larger
  alpha2 = 0.05,  # positive shock reaction — smaller
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

RTgarch11 <- function(par, rt, mu) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  psi <- par[4]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || psi < 0 || beta + psi >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- omega / (1 - beta - psi)
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  
  
  
  # Generate conditional variances recursively
  for (t in 1:(T - 1)) {
    h[t + 1] <- 0.5*(omega + beta*h[t] + alpha*(rt[t] - mu)^2) + 0.5*sqrt((omega + beta*h[t] + alpha*(rt[t] - mu)^2)^2 + 4*psi*h[t]*(rt[t+1] - mu)^2)
  }
  
  # Now check h after it has been generated
  if (any(h <= 0) || any(is.na(h)) || any(is.infinite(h))) {
    return(1e10)
  }
  
  
  h_tm1 <- h[-T]
  h <- h[-1]
  rt <- rt[-1]
  
  # Negative log-likelihood
  RTgarch11ll <- sum(0.5*log(2*pi)+0.5*(rt-mu)^2/h-log(sqrt(h)/(h+psi*h_tm1*(rt-mu)^2/h)))
  
  
  if (!is.finite(RTgarch11ll)) {
    return(1e10)
  }
  
  return(RTgarch11ll)
}


start_par_RT <- c(
  omega = 0.04590705,
  alpha = 0.19226977 ,
  beta  = 0.77380626 ,
  psi = 0.01
)

RTgarch_fit <- optim(
  par = start_par_RT,
  fn = RTgarch11,
  rt = rt,
  mu = mu,
  method = "L-BFGS-B",
  #Omega > 0 to inf, remaining are bounded by 0 and 1 
  lower = c(1e-8, 0, 0, 0),
  upper = c(Inf, 1, 1, 1)
)

RTgarch_fit$par
RTgarch_fit$value
RTgarch_fit$convergence

h_long = (as.numeric(RTgarch_fit$par[1]))/(as.numeric(1 - RTgarch_fit$par[3] - RTgarch_fit$par[4]))

#model implied unconditional variance 
print(h_long)

RTgjr_garch <- function(par, rt, mu) {
  
  omega  <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  beta   <- par[4]
  phi1   <- par[5]
  phi2   <- par[6]
  
  phi_bar   <- (phi1 + phi2) / 2
  alpha_bar <- (alpha1 + alpha2) / 2
  
  if (omega <= 0 || alpha1 < 0 || alpha2 < 0 || beta < 0 ||
      phi1 < 0 || phi2 < 0 ||
      beta + phi_bar >= 1) return(1e10)
  
  T <- length(rt)
  h <- numeric(T)
  h[1] <- omega / (1 - beta - phi_bar)
  
  if (!is.finite(h[1]) || h[1] <= 0) return(1e10)
  
  for (t in 1:(T-1)) {
    alpha_t <- ifelse(rt[t]   <= mu, alpha1, alpha2)
    phi_t   <- ifelse(rt[t+1] <= mu, phi1,   phi2)
    g_t     <- omega + beta*h[t] + alpha_t*(rt[t]-mu)^2
    h[t+1]  <- 0.5*g_t + 0.5*sqrt(g_t^2 + 4*phi_t*h[t]*(rt[t+1]-mu)^2)
  }
  
  if (any(h <= 0) || any(is.na(h)) || any(is.infinite(h))) return(1e10)
  
  h_tm1   <- h[-T]
  h_cur   <- h[-1]
  rt_cur  <- rt[-1]
  phi_vec <- ifelse(rt_cur <= mu, phi1, phi2)
  
  ll <- sum(
    0.5*log(2*pi) +
      0.5*(rt_cur-mu)^2/h_cur -
      log(sqrt(h_cur)/(h_cur + phi_vec*h_tm1*(rt_cur-mu)^2/h_cur))
  )
  
  if (!is.finite(ll)) return(1e10)
  return(ll)
}

start_par_RTgjr <- c(
  omega  = RTgarch_fit$par[1],
  alpha1 = RTgarch_fit$par[2] * 1.2,
  alpha2 = RTgarch_fit$par[2] * 0.8,
  beta   = RTgarch_fit$par[3],
  phi1   = RTgarch_fit$par[4] * 1.1,
  phi2   = RTgarch_fit$par[4] * 0.9
)

start_par_RTgjr <- c(
  omega  = 0.05,
  alpha1 = 0.08,
  alpha2 = 0.02,
  beta   = 0.75,
  phi1   = 0.10,
  phi2   = 0.08
)

# RTgjr_fit <- optim(
#   par     = start_par_RTgjr,
#   fn      = RTgjr_garch,
#   rt      = rt,
#   mu      = mu,
#   method  = "L-BFGS-B",
#   lower   = c(1e-8, 0, 0, 0, 0, 0),
#   upper   = c(Inf,  1, 1, 1, 1, 1),
#   control = list(maxit = 1000, factr = 1e7)
# )

RTgjr_fit <- optim(
  par     = start_par_RTgjr,
  fn      = RTgjr_garch,
  rt      = rt,
  mu      = mu,
  method  = "L-BFGS-B",
  lower   = c(1e-8, 0, 1e-8, 0, 0, 1e-8),
  upper   = c(Inf,  1, 1,    1, 1, 1),
  control = list(maxit = 1000)
)

RTgjr_fit$par
RTgjr_fit$value
RTgjr_fit$convergence

#Information Criterion:

#Loop for best starting values:

bestValuesGARCH <- c(omega = 0,
                alpha = 0,
                beta = 0,
                llvalue = Inf)

bestValuesGJRGARCH <- c(omega = 0,
                     alpha1 = 0,
                     alpha2 = 0,
                     beta = 0,
                     llvalue = Inf)

ValOptimG <- numeric(6^4)
ValOptimGGJR <- numeric(6^4)

OmeOptimG <- numeric(6^4)
OmeOptimGGJR <- numeric(6^4)

Alp1OptimG <- numeric(6^4)
alp1OptimGGJR <- numeric(6^4)

Alp2OptimG <- numeric(6^4)
Alp2OptimGGJR <- numeric(6^4)

AlpOptimG <- numeric(6^4)
AlpOptimGGJR <- numeric(6^4)

BetOptimG <- numeric(6^4)
BetOptimGGJR <- numeric(6^4)

num <- 0

for (o in seq(from = 0, to=1, by=0.2)) {
  for(a1 in seq(from=0, to=1, by=0.2)) {
   for(a2 in seq(from=0, to=1, by=0.2)) {
     for(b in seq(from=0, to=1, by=0.2)){
       
       start_par <- c(
         omega = o,
         alpha = (a1+a2)/2,
         beta  = b
       )
       
       start_par_gjr <- c(
         omega  = o,
         alpha1 = a1,  # negative shock reaction — larger
         alpha2 = a2,  # positive shock reaction — smaller
         beta   = b
       )
       
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
       
       gjr_fit <- optim(
         par = start_par_gjr,
         fn = garch_gjr,
         rt = rt,
         mu = mu,
         method = "L-BFGS-B",
         lower = c(1e-8, 0, 0, 0),
         upper = c(10, 1, 1, 1)
       )
       
       
       
       ValOptimG[num] <- garch_fit$value
       ValOptimGGJR[num] <- gjr_fit$value
       
       OmeOptimG[num] <- garch_fit$par["omega"]
       OmeOptimGGJR[num] <- gjr_fit$par["omega"]
       
       alp1OptimGGJR[num] <- gjr_fit$par["alpha1"]
       
       Alp2OptimGGJR[num] <- gjr_fit$par["alpha2"]
       
       AlpOptimG[num] <- garch_fit$par["alpha"]
       
       BetOptimG[num] <- garch_fit$par["beta"]
       BetOptimGGJR[num] <- gjr_fit$par["beta"]
       
       
       
       if (garch_fit$value < bestValuesGARCH["llvalue"]) {
         bestValuesGARCH["omega"] <- garch_fit$par["omega"]
         bestValuesGARCH["alpha"] <- garch_fit$par["alpha"]
         bestValuesGARCH["beta"] <- garch_fit$par["beta"]
         bestValuesGARCH["llvalue"] <- garch_fit$value
       }
       
       if (gjr_fit$value < bestValuesGJRGARCH["llvalue"]) {
         bestValuesGJRGARCH["omega"] <- gjr_fit$par["omega"]
         bestValuesGJRGARCH["alpha1"] <- gjr_fit$par["alpha1"]
         bestValuesGJRGARCH["alpha2"] <- gjr_fit$par["alpha2"]
         bestValuesGJRGARCH["beta"] <- gjr_fit$par["beta"]
         bestValuesGJRGARCH["llvalue"] <- gjr_fit$value
       }
       
       
       
       
       num <- num + 1
       
       print(paste("Iterations completed: ", num))
    }  
   }
  }
}


hist(ValOptimG[ValOptimG > 0])
hist(ValOptimGGJR[ValOptimGGJR > 0])

hist(OmeOptimG[OmeOptimG > 0])
hist(OmeOptimGGJR[OmeOptimGGJR > 0])

hist(alp1OptimGGJR[alp1OptimGGJR > 0])
hist(Alp2OptimGGJR[Alp2OptimGGJR > 0])
hist(AlpOptimG[AlpOptimG > 0])

hist(BetOptimG[BetOptimG > 0])
hist(BetOptimGGJR[BetOptimGGJR > 0])


print("For Garch:")

bestValuesGARCH["omega"]
bestValuesGARCH["alpha"]
bestValuesGARCH["beta"]
bestValuesGARCH["llvalue"]

print("For GJR Garch:")

bestValuesGJRGARCH["omega"]
bestValuesGJRGARCH["alpha1"]
bestValuesGJRGARCH["alpha2"]
bestValuesGJRGARCH["beta"]
bestValuesGJRGARCH["llvalue"]



AIC_GARCH    <- 2*3 + 2*garch_fit$value      # 3 params
AIC_GJR      <- 2*4 + 2*gjr_fit$value        # 4 params
AIC_RTGARCH  <- 2*4 + 2*RTgarch_fit$value    # 4 params
AIC_RTgjr    <- 2*6 + 2*RTgjr_fit$value      # 6 params

print(c(GARCH    = AIC_GARCH,
        GJR      = AIC_GJR,
        RT_GARCH = AIC_RTGARCH,
        RT_GJR   = AIC_RTgjr))


