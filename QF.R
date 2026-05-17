#Read excel imports
install.packages("readxl")
library("readxl")

#For skewness & kurtosis analysis
install.packages("GGally")
library("GGally")

install.packages("e1071")
library(e1071)

library(tidyr)

install.packages("numDeriv")
library(numDeriv)

#--------------------------DATA PROCESSING --------------------------

# Janek:
# df_init <- read_excel("/Users/janekczajnik/Desktop/Erasmus/B3/Introductory Seminar CS/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

# Balint
#df_init <- read_excel("/Users/balintkovacs/Documents/GitHub/QF-case/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

# Luca
setwd("C:/Users/lucam/Dropbox/dad&mum/University/Erasmus/BSC 3/BLK 5/Intro to Seminars/Case2_QF/econometrics & QF case study 2026/data")

#Filip
#df_init <- read_excel("C:/Users/filip/Downloads/econometrics & QF case study 2026/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

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

EXPANDEDRTgarch11 <- function(par, rt, mu) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  psi <- par[4]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || psi < 0 || beta + psi + alpha + 2*psi*alpha >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- omega / (1 - beta - psi - alpha - 2*psi*alpha)
  
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

EXPANDEDRTgarch_fit <- optim(
  par = start_par_RT,
  fn = EXPANDEDRTgarch11,
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

EXPANDEDRTgarch_fit$par
EXPANDEDRTgarch_fit$value
EXPANDEDRTgarch_fit$convergence


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

EXPANDEDRTgjr_garch <- function(par, rt, mu) {
  
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
      beta + phi_bar + alpha_bar - alpha_bar*phi_bar + 3/2*(alpha1*phi1 + alpha2*phi2) >= 1) return(1e10)
  
  T <- length(rt)
  h <- numeric(T)
  h[1] <- omega / (1 - beta - phi_bar - alpha_bar + alpha_bar*phi_bar - 3/2*(alpha1*phi1 + alpha2*phi2))
  
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
  omega  = 0.04149673,
  alpha1 = 0.26501891,
  alpha2 = 0.01882367,
  beta   = 0.81790257,
  phi1   = 0.015,
  phi2   = 0.005
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
  lower   = c(1e-8, 0, 0, 0, 0, 0),
  upper   = c(Inf,  1, 1, 1, 1, 1),
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

BIC_GARCH    <- 3*log(T-1) + 2*garch_fit$value      # 3 params
BIC_GJR      <- 4*log(T-1) + 2*gjr_fit$value        # 4 params
BIC_RTGARCH  <- 4*log(T-1) + 2*RTgarch_fit$value    # 4 params
BIC_RTgjr    <- 6*log(T-1) + 2*RTgjr_fit$value      # 6 params

print(c(GARCH    = BIC_GARCH,
        GJR      = BIC_GJR,
        RT_GARCH = BIC_RTGARCH,
        RT_GJR   = BIC_RTgjr))

#--------------------------MARKET BREAK DUMMIES-----------------------------------


CompleteYear <- c()
CompleteMonth <- c()
CompleteDay <- c()

for (month in 10:12){
  if(month %in% c(11)){
    for(day in 1:30){
      CompleteYear <- c(CompleteYear, 2009)
      CompleteMonth <- c(CompleteMonth, month)
      CompleteDay <- c(CompleteDay, day)
    }
  } else {
    for(day in 1:31){
      CompleteYear <- c(CompleteYear, 2009)
      CompleteMonth <- c(CompleteMonth, month)
      CompleteDay <- c(CompleteDay, day)
    }
  }
}


for(year in 2010:2025){
  leap <- 0
  
  if(year %in% c(2010, 2016, 2020, 2024)){
    leap <- 1
  }
  
  for (month in 1:12){
    if(month %in% c(4,6,9,11)){
      for(day in 1:30){
        CompleteYear <- c(CompleteYear, year)
        CompleteMonth <- c(CompleteMonth, month)
        CompleteDay <- c(CompleteDay, day)
      }
    } else if (month %in% c(1,3,5,7,8,10,12)) {
      for(day in 1:31){
        CompleteYear <- c(CompleteYear, year)
        CompleteMonth <- c(CompleteMonth, month)
        CompleteDay <- c(CompleteDay, day)
      }
    }else{
      for(day in 1:(28+leap)){
        CompleteYear <- c(CompleteYear, year)
        CompleteMonth <- c(CompleteMonth, month)
        CompleteDay <- c(CompleteDay, day)
      }
    }
  }
}

for (month in 1:3){
  if (month %in% c(1,3)) {
    for(day in 1:31){
      CompleteYear <- c(CompleteYear, 2026)
      CompleteMonth <- c(CompleteMonth, month)
      CompleteDay <- c(CompleteDay, day)
    }
  }else{
    for(day in 1:28){
      CompleteYear <- c(CompleteYear, 2026)
      CompleteMonth <- c(CompleteMonth, month)
      CompleteDay <- c(CompleteDay, day)
    }
  }
}

CompleteMonth <-sprintf("%02d", CompleteMonth)
CompleteDay <-sprintf("%02d", CompleteDay)


CompleteCalender <- as.data.frame(cbind(CompleteYear, CompleteMonth, CompleteDay))


CompleteCalender <- unite(CompleteCalender,"Date" ,CompleteYear, CompleteMonth, CompleteDay, sep = "-")

DayInData <- integer(nrow(CompleteCalender))
AfterHoliday <- integer(nrow(CompleteCalender))
for(i in 1:nrow(CompleteCalender)){
  DayInData[i] <- ifelse(CompleteCalender$Date[i] %in% substring(df_init$Date,1,10), 1, 0)
  if (i > 1) {
    if (DayInData[i] == 1 & DayInData[i-1] == 0) {
      AfterHoliday[i] <- 1
    }
  }
}

CompleteCalender$DayInData <- DayInData
CompleteCalender$AfterHoliday <- AfterHoliday

DummyAfterHoliday <- integer(T)

for(i in 1:T) {
  if(substring(df_init$Date,1,10)[i] %in% CompleteCalender$Date){
    DummyAfterHoliday[i] <- CompleteCalender[CompleteCalender$Date == substring(df_init$Date,1,10)[i], 3]
  }
}
df_init$DummyAfterHoliday <- DummyAfterHoliday


lengthOfHolidays <- c()
k <- 0


for(i in 3:nrow(CompleteCalender)){
  if(CompleteCalender$DayInData[i] == 0){
    k <- k + 1
  } else if(CompleteCalender$DayInData[i] == 1 & 
            CompleteCalender$DayInData[i-1] == 0) {
    lengthOfHolidays <- c(lengthOfHolidays, k)
    k <- 0
  } else{
    k <- 0
  }
}

AllPrices <- c()

for(i in 1:4178){
  AllPrices <- c(AllPrices, df_init$Open[i])
  AllPrices <- c(AllPrices, df_init$Close[i])
}

OpenMarketReturns <- 0
CloseMarketReturns <- 0

for(i in 2:length(AllPrices)){
  k <- AllPrices[i] - AllPrices[i-1]
  if((i %% 2) == 0){
    OpenMarketReturns <- OpenMarketReturns + k
  } else {
    CloseMarketReturns <- CloseMarketReturns + k
  }
}

OO_rt <- as.numeric(100*log(df_init$Open[-1]/df_init$Open[-T]))

DummyRTgarch11 <- function(par, rt, mu, DUM) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  psi <- par[4]
  delta <- par[5]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || psi < 0 || beta + psi + alpha + delta*mean(DUM) + 2*psi*(alpha + delta*mean(DUM)) >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- omega / (1 - beta - psi - alpha - delta*mean(DUM) - 2*psi*(alpha + delta*mean(DUM)))
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  # Generate conditional variances recursively
  for (t in 1:(T - 1)) {
    h[t + 1] <- 0.5*(omega + beta*h[t] + (alpha+delta*DUM[t+1])*(rt[t] - mu)^2) + 0.5*sqrt((omega + beta*h[t] + (alpha+delta*DUM[t+1])*(rt[t] - mu)^2)^2 + 4*(psi)*h[t]*(rt[t+1] - mu)^2)
  }
  
  # Now check h after it has been generated
  if (any(h <= 0) || any(is.na(h)) || any(is.infinite(h))) {
    return(1e10)
  }
  
  
  h_tm1 <- h[-T]
  h <- h[-1]
  rt <- rt[-1]
  DUM <- DUM[-1]
  
  # Negative log-likelihood
  RTgarch11ll <- sum(0.5*log(2*pi)+0.5*(rt-mu)^2/h-log(sqrt(h)/(h+(psi)*h_tm1*(rt-mu)^2/h)))
  
  
  if (!is.finite(RTgarch11ll)) {
    return(1e10)
  }
  
  return(RTgarch11ll)
}


DUMMYRTgjr_garch <- function(par, rt, mu, DUM) {
  
  omega  <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  beta   <- par[4]
  phi1   <- par[5]
  phi2   <- par[6]
  delta  <- par[7]
  
  phi_bar   <- (phi1 + phi2) / 2
  alpha_bar <- (alpha1 + alpha2) / 2
  
  if (omega <= 0 || alpha1 < 0 || alpha2 < 0 || beta < 0 ||
      phi1 < 0 || phi2 < 0 ||
      beta + phi_bar + alpha_bar + delta*mean(DUM) - (alpha_bar + delta*mean(DUM))*phi_bar + 3/2*(alpha1*phi1 + alpha2*phi2) + 3*delta*mean(DUM)*phi_bar >= 1) return(1e10)

  T <- length(rt)
  h <- numeric(T)
  h[1] <- omega / (1 - beta - phi_bar - alpha_bar - delta*mean(DUM) + (alpha_bar + delta*mean(DUM))*phi_bar - 3/2*(alpha1*phi1 + alpha2*phi2) - 3*delta*mean(DUM)*phi_bar)
  
  if (!is.finite(h[1]) || h[1] <= 0) return(1e10)
  
  for (t in 1:(T-1)) {
    alpha_t <- ifelse(rt[t]   <= mu, alpha1, alpha2)
    phi_t   <- ifelse(rt[t+1] <= mu, phi1,   phi2)
    g_t     <- omega + beta*h[t] + (alpha_t+delta*DUM[t+1])*(rt[t]-mu)^2
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


start_par_RT <- c(
  omega = 0.04590705,
  alpha = 0.19226977 ,
  beta  = 0.77380626 ,
  psi = 0.01,
  delta = 0.01
)

DummyRTgarch_fit <- optim(
  par = start_par_RT,
  fn = DummyRTgarch11,
  rt = rt,
  mu = mu,
  DUM = DummyAfterHoliday,
  method = "L-BFGS-B",
  #Omega > 0 to inf, remaining are bounded by 0 and 1 
  lower = c(1e-8, 0, 0, 0),
  upper = c(Inf, 1, 1, 1)
)

start_par_RTgjr <- c(
  omega  = 0.04149673,
  alpha1 = 0.26501891,
  alpha2 = 0.01882367,
  beta   = 0.81790257,
  phi1   = 0.015,
  phi2   = 0.005,
  delta = 0.01
)


DummyRTgarchGJR_fit <- optim(
  par = start_par_RTgjr,
  fn = DUMMYRTgjr_garch,
  rt = rt,
  mu = mu,
  DUM = DummyAfterHoliday,
  method = "L-BFGS-B",
  #Omega > 0 to inf, remaining are bounded by 0 and 1 
  lower   = c(1e-8, 0, 0, 0, 0, 0),
  upper   = c(Inf,  1, 1, 1, 1, 1),
  control = list(maxit = 1000)
)





DummyRTgarch_fit$par
DummyRTgarch_fit$value
DummyRTgarch_fit$convergence

DummyRTgarchGJR_fit$par
DummyRTgarchGJR_fit$value
DummyRTgarchGJR_fit$convergence


#--------------------------HESSIANS--------------------------------------------

rt <- df_init$`CC Return (%)`
mu <- mean(rt)
T <- length(rt)


FischerGARCH <- matrix(0,3,3)
for (i in 1:T){

  ObsLLGARCH <- function(par) {
    
    omega <- par[1]
    alpha <- par[2]
    beta  <- par[3]
    
    h <- numeric(i)
    
    # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
    h[1] <- omega/(1-alpha-beta)
    

    
    
    if (i > 1){
      for (t in 2:i){
        h[t] <- omega + alpha * (rt[t-1] - mu)^2 + beta * h[t-1]
      }
    }
    
    # log-likelihood
    garch11ll <- 0.5 * (-log(2 * pi) - log(h[i]) - ((rt[i] - mu)^2 / h[i]))
    
    return(garch11ll)
  }
  print(i)
  
  FischerGARCH <- FischerGARCH + hessian(ObsLLGARCH, garch_fit$par, method = "complex")
}

FischerGARCH <- solve(-FischerGARCH)



#--------------------------h_t PLOTTING ---------------------------------------

#GARCH recursion

h_1 <- numeric(T)
h_1[1] <- var(rt)

omega <- garch_fit$par[1]
alpha <- garch_fit$par[2]
beta <- garch_fit$par[3]

for (t in 1:(T - 1)) {
  h_1[t + 1] <- omega + alpha*(rt[t] - mu)^2 + beta*h_1[t]
}


#GARCH-GJR recursion

h_garch <- numeric(T)
h_garch[1] <- var(rt)

omega <- gjr_fit$par[1]
alpha1 <- gjr_fit$par[2]
alpha2 <- gjr_fit$par[3]
beta <- gjr_fit$par[4]

shock <- rt - mu

for (t in 1:(T - 1)) {
  
  indicator <- ifelse(shock[t] < 0, 1, 0)
  
  h_garch[t+1] <- omega +
    alpha1 * indicator * shock[t]^2 +        # alpha1 when negative
    alpha2 * (1 - indicator) * shock[t]^2 +  # alpha2 when positive
    beta * h_garch[t]
}


#RT-GARCH recursion

h_htgarch <- numeric(T)
g__htgarch <- numeric(T)

omega <- RTgarch_fit$par[1]
alpha <- RTgarch_fit$par[2]
beta <- RTgarch_fit$par[3]
psi <- RTgarch_fit$par[4]

h_htgarch[1] <- omega / (1 - beta - psi)

for (t in 1:(T - 1)) {
  h_htgarch[t + 1] <- 0.5*(omega + beta*h_htgarch[t] + alpha*(rt[t] - mu)^2) + 0.5*sqrt((omega + beta*h_htgarch[t] + alpha*(rt[t] - mu)^2)^2 + 4*psi*h_htgarch[t]*(rt[t+1] - mu)^2)
  g__htgarch[t] <- omega + beta*h_htgarch[t] + alpha*(rt[t] - mu)^2
}

#RT-GARCH-GJR recursion

h_htgarchgrj <- numeric(T)
g_htgarchgrj <- numeric(T)

omega  <- RTgjr_fit$par[1]
alpha1 <- RTgjr_fit$par[2]
alpha2 <- RTgjr_fit$par[3]
beta   <- RTgjr_fit$par[4]
phi1   <- RTgjr_fit$par[5]
phi2   <- RTgjr_fit$par[6]

phi_bar   <- (phi1 + phi2) / 2
alpha_bar <- (alpha1 + alpha2) / 2

h_htgarchgrj[1] <- omega / (1 - beta - phi_bar)

for (t in 1:(T-1)) {
  alpha_t <- ifelse(rt[t]   <= mu, alpha1, alpha2)
  phi_t   <- ifelse(rt[t+1] <= mu, phi1,   phi2)
  g_t     <- omega + beta*h_htgarchgrj[t] + alpha_t*(rt[t]-mu)^2
  h_htgarchgrj[t+1]  <- 0.5*g_t + 0.5*sqrt(g_t^2 + 4*phi_t*h_htgarchgrj[t]*(rt[t+1]-mu)^2)
  g_htgarchgrj[t] <- g_t
}

#actual plots

plot(h_1)
plot(h_garch)
plot(h_htgarch)
plot(h_htgarchgrj)

# 1. Global Y-limit for consistency
y_limit <- c(0, max(c(h_1, h_garch, h_htgarch, h_htgarchgrj), na.rm = TRUE))

# 2. Define Period 1
period1 <- 400:600

# 3. Plotting Grid 1 (Period 1)
par(mfrow = c(2, 2))

plot(period1, h_1[period1], type = "l", col = "blue", ylim = y_limit,
     main = "GARCH (400-600)", xlab = "Time", ylab = "h_t")

plot(period1, h_garch[period1], type = "l", col = "red", ylim = y_limit,
     main = "GARCH-GJR (400-600)", xlab = "Time", ylab = "h_t")

plot(period1, h_htgarch[period1], type = "l", col = "darkgreen", ylim = y_limit,
     main = "RT-GARCH (400-600)", xlab = "Time", ylab = "h_t")

plot(period1, h_htgarchgrj[period1], type = "l", col = "purple", ylim = y_limit,
     main = "RT-GARCH-GJR (400-600)", xlab = "Time", ylab = "h_t")

par(mfrow = c(1, 1)) # Reset

# 4. Define Period 2
period2 <- 2000:3000

# 5. Plotting Grid 2 (Period 2)
par(mfrow = c(2, 2))

plot(period2, h_1[period2], type = "l", col = "blue", ylim = y_limit,
     main = "GARCH (2000-3000)", xlab = "Time", ylab = "h_t")

plot(period2, h_garch[period2], type = "l", col = "red", ylim = y_limit,
     main = "GARCH-GJR (2000-3000)", xlab = "Time", ylab = "h_t")

plot(period2, h_htgarch[period2], type = "l", col = "darkgreen", ylim = y_limit,
     main = "RT-GARCH (2000-3000)", xlab = "Time", ylab = "h_t")

plot(period2, h_htgarchgrj[period2], type = "l", col = "purple", ylim = y_limit,
     main = "RT-GARCH-GJR (2000-3000)", xlab = "Time", ylab = "h_t")

par(mfrow = c(1, 1)) # Reset
#------------------------------------------ESTIMATION/EVALUATION SPLIT---------------------------------------

#Define the Estimation Split
df_init$Date <- as.Date(df_init$Date) 
T_est <- which(df_init$Date == as.Date("2019-12-31"))
if(length(T_est) == 0) T_est <- max(which(df_init$Date <= as.Date("2019-12-31")))
rt_est <- rt[1:T_est]
mu_est <- mean(rt_est)

#Reestimate GARCH
garch_est_fit <- optim(
  par = start_par,
  fn = garch11,
  rt = rt_est,
  mu = mu_est,
  method = "L-BFGS-B",
  lower = c(1e-8, 0, 0),
  upper = c(Inf, 1, 1)
)

#Reestimate GJR-GARCH
gjr_est_fit <- optim(
  par = start_par_gjr,
  fn = garch_gjr,
  rt = rt_est,
  mu = mu_est,
  method = "L-BFGS-B",
  lower = c(1e-8, 0, 0, 0),
  upper = c(10, 1, 1, 1)
)

#Reestimate RT-GARCH
RTgarch_est_fit <- optim(
  par = start_par_RT,
  fn = RTgarch11,
  rt = rt_est,
  mu = mu_est,
  method = "L-BFGS-B",
  #Omega > 0 to inf, remaining are bounded by 0 and 1 
  lower = c(1e-8, 0, 0, 0),
  upper = c(Inf, 1, 1, 1)
)



#Reestimate RT-GJR-GARCH
RTgjr_est_fit <- optim(
  par  = start_par_RTgjr,
  fn   = RTgjr_garch,
  rt   = rt_est,
  mu   = mu_est,
  method  = "L-BFGS-B",
  lower   = c(1e-8, 0, 1e-8, 0, 0, 1e-8),
  upper   = c(Inf,  1, 1,    1, 1, 1),
  control = list(maxit = 1000)
)

#------------------------RECURSIONS FOR ESTIMATION-----------------------------------------------------------

#GARCH recursion

h_1 <- numeric(T)
h_1[1] <- var(rt_est)

omega <- garch_est_fit$par[1]
alpha <- garch_est_fit$par[2]
beta <- garch_est_fit$par[3]

for (t in 1:(T - 1)) {
  h_1[t + 1] <- omega + alpha*(rt[t] - mu_est)^2 + beta*h_1[t]
}


#GARCH-GJR recursion

h_garch <- numeric(T)
h_garch[1] <- var(rt_est)

omega <- gjr_est_fit$par[1]
alpha1 <- gjr_est_fit$par[2]
alpha2 <- gjr_est_fit$par[3]
beta <- gjr_est_fit$par[4]

shock <- rt - mu_est

for (t in 1:(T - 1)) {
  
  indicator <- ifelse(shock[t] < 0, 1, 0)
  
  h_garch[t+1] <- omega +
    alpha1 * indicator * shock[t]^2 +        # alpha1 when negative
    alpha2 * (1 - indicator) * shock[t]^2 +  # alpha2 when positive
    beta * h_garch[t]
}


#RT-GARCH recursion

h_htgarch <- numeric(T)
g_htgarch <- numeric(T)

omega <- RTgarch_est_fit$par[1]
alpha <- RTgarch_est_fit$par[2]
beta <- RTgarch_est_fit$par[3]
psi <- RTgarch_est_fit$par[4]

h_htgarch[1] <- omega / (1 - beta - psi)

for (t in 1:(T - 1)) {
  h_htgarch[t + 1] <- 0.5*(omega + beta*h_htgarch[t] + alpha*(rt[t] - mu_est)^2) + 0.5*sqrt((omega + beta*h_htgarch[t] + alpha*(rt[t] - mu_est)^2)^2 + 4*psi*h_htgarch[t]*(rt[t+1] - mu_est)^2)
  g_htgarch[t] <- omega + beta*h_htgarch[t] + alpha*(rt[t] - mu_est)^2
}

#RT-GARCH-GJR recursion

h_htgarchgjr <- numeric(T)
g_htgarchgjr <- numeric(T)

omega  <- RTgjr_est_fit$par[1]
alpha1 <- RTgjr_est_fit$par[2]
alpha2 <- RTgjr_est_fit$par[3]
beta   <- RTgjr_est_fit$par[4]
phi1   <- RTgjr_est_fit$par[5]
phi2   <- RTgjr_est_fit$par[6]

phi_bar   <- (phi1 + phi2) / 2
alpha_bar <- (alpha1 + alpha2) / 2

h_htgarchgjr[1] <- omega / (1 - beta - phi_bar)

for (t in 1:(T-1)) {
  alpha_t <- ifelse(rt[t]   <= mu_est, alpha1, alpha2)
  phi_t   <- ifelse(rt[t+1] <= mu_est, phi1,   phi2)
  g_t     <- omega + beta*h_htgarchgjr[t] + alpha_t*(rt[t]-mu_est)^2
  h_htgarchgjr[t+1]  <- 0.5*g_t + 0.5*sqrt(g_t^2 + 4*phi_t*h_htgarchgjr[t]*(rt[t+1]-mu_est)^2)
  g_htgarchgjr[t] <- g_t
}




#-----------------------------------------TARGET VALUES-----------------------------------------------------

num_forecast <- T - T_est

# Demeaned squared return
TV_r <- numeric(num_forecast)
indx <- 1

for (t in T_est:(T - 1)){
  TV_r[indx] <- (rt[t+1] - mu_est)^2
  indx <- indx + 1
}

#Realised variance 
TV_rv <- numeric(num_forecast)
indx <- 1

for (t in T_est:(T - 1)){
  TV_rv[indx] <- df_init$`RV5_SS × 10^4`[t + 1]
  indx <- indx + 1
}



#-----------PREDICTED VALUES-----------------------------------------------

#one-day GARCH forecast
PV_GARCH <- numeric(num_forecast)
indx <- 1

for (t in T_est:(T - 1)) {
  PV_GARCH[indx] <- h_1[t + 1]
  indx <- indx + 1
}

#one-day GARCH-GJR forecast
PV_GARCHGJR <- numeric(num_forecast)
indx <- 1

for (t in T_est:(T - 1)) {
  PV_GARCHGJR[indx] <- h_garch[t + 1]
  indx <- indx + 1
}

#one-day RT-GARCH forecast
PV_RTGARCH <- numeric(num_forecast)
kurtosis <- 3
psi <- RTgarch_est_fit$par[4]
indx <- 1

for (t in T_est:(T - 1)) {
  PV_RTGARCH[indx] <- g_htgarch[t] + h_htgarch[t] * kurtosis * psi
  indx <- indx + 1
  
}

#one-day RT-GARCH_GJR forecast
PV_RTGARCHGJR <- numeric(num_forecast)
kurtosis <- 3
phi1   <- RTgjr_est_fit$par[5]
phi2   <- RTgjr_est_fit$par[6]
phi_bar   <- (phi1 + phi2) / 2
indx <- 1

for (t in T_est:(T - 1)) {
  PV_RTGARCHGJR[indx] <- g_htgarchgjr[t] + h_htgarchgjr[t] * kurtosis * phi_bar
  indx <- indx + 1
  
}


#one-day VIX forecast
PV_VIX <- numeric(num_forecast)
indx <- 1

for (t in T_est:(T - 1)) {
  PV_VIX[indx] <- (df_init$VIX[t])^2 / 250
  indx <- indx + 1
}

#HAR-RV forecast
RV_day <- numeric(T)
RV_week <- numeric(T)
RV_month <- numeric(T)

for (t in 1:T){
  RV_day[t] <- df_init$`RV5_SS × 10^4`[t]
  if (t > 4){
    RV_week[t] <- sum(RV_day[(t-4):t]) / 5
  }
  if (t > 20){
    RV_month[t] <- sum(RV_day[(t-20):t]) / 21
  }
}

PV_HAR_RV <- numeric(num_forecast)
PV_HAR_R  <- numeric(num_forecast)

RV_d_est <- RV_day[21:(T_est - 1)]
RV_w_est <- RV_week[21:(T_est - 1)]
RV_m_est <- RV_month[21:(T_est - 1)]

TV_rv_est <- RV_day[22:T_est]
TV_r_est  <- (rt[22:T_est] - mu_est)^2

har_data <- data.frame(RV_d_est,RV_w_est, RV_m_est, TV_rv_est, TV_r_est)

model_rv <- lm(TV_rv_est ~ RV_d_est + RV_w_est + RV_m_est, data = har_data)
model_r <- lm(TV_r_est ~ RV_d_est + RV_w_est + RV_m_est, data = har_data)

indx <- 1

for (t in T_est:(T - 1)) {
  PV_HAR_RV[indx] <- coef(model_rv)[1] + coef(model_rv)[2]*RV_day[t] + coef(model_rv)[3]*RV_week[t] + coef(model_rv)[4]*RV_month[t]
  PV_HAR_R[indx] <- coef(model_r)[1] + coef(model_r)[2]*RV_day[t] + coef(model_r)[3]*RV_week[t] + coef(model_r)[4]*RV_month[t]
  indx <- indx + 1
}

#-------------------------------LOSS FUNCTIONS----------------------------------------

#For cases when forecasts are nonpositive 
safety_val <- 1e-10
PV_HAR_RV[PV_HAR_RV <= 0] <- safety_val
PV_HAR_R[PV_HAR_R <= 0] <- safety_val
PV_GARCH[PV_GARCH <= 0] <- safety_val
PV_GARCHGJR[PV_GARCHGJR <= 0] <- safety_val
PV_RTGARCH[PV_RTGARCH <= 0] <- safety_val
PV_RTGARCHGJR[PV_RTGARCHGJR <= 0] <- safety_val
PV_VIX[PV_VIX <= 0] <- safety_val

#Loss function calculations
calc_losses <- function(PV, TV) {
  evaluation_size <- length(PV)
  
  # Mean Squared Error
  mse <- (1/evaluation_size) * sum((TV - PV)^2)
  
  # QLIKE Loss
  qlike <- (1/evaluation_size) * sum(log(PV) + (TV / PV))
  
  return(c(MSE = mse, QLIKE = qlike))
}


# List of all predicted values
predictions <- list(
  GARCH = PV_GARCH,
  GJR_GARCH = PV_GARCHGJR,
  RT_GARCH = PV_RTGARCH,
  RT_GJR = PV_RTGARCHGJR,
  VIX = PV_VIX,
  HAR = PV_HAR_RV
)


loss_rv <- sapply(predictions, function(p) calc_losses(p, TV_rv))
predictions$HAR <- PV_HAR_R 
loss_r <- sapply(predictions, function(p) calc_losses(p, TV_r))


final_results <- rbind(
  Target_RV_MSE   = loss_rv["MSE", ],
  Target_RV_QLIKE = loss_rv["QLIKE", ],
  Target_R_MSE    = loss_r["MSE", ],
  Target_R_QLIKE  = loss_r["QLIKE", ]
)

print(round(final_results, 6))
