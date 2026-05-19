#Read excel imports
install.packages("readxl")
library("readxl")

#For skewness & kurtosis analysis
install.packages("GGally")
library(GGally)

install.packages("e1071")
library(e1071)

library(tidyr)

install.packages("numDeriv")
library(numDeriv)

install.packages("knitr")
library(knitr)

install.packages("nnet")
library(nnet)


#The coolest graphs ever! 

library(ggplot2)
install.packages("ggplot2")
install.packages("ggplot")

install.packages("reshape2")
library(reshape2)
install.packages("ggthemes")
library(ggthemes)


install.packages("nnet")
library(nnet)

install.packages("forecast")
library(forecast)



#--------------------------DATA PROCESSING --------------------------

# Janek:
# df_init <- read_excel("/Users/janekczajnik/Desktop/Erasmus/B3/Introductory Seminar CS/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

# Balint
#df_init <- read_excel("/Users/balintkovacs/Documents/GitHub/QF-case/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

# Luca
#df_init <- read_excel("C:/Users/lucam/Dropbox/dad&mum/University/Erasmus/BSC 3/BLK 5/Intro to Seminars/Case2_QF/econometrics & QF case study 2026/data/data.xlsx")[-1, ]
library("readxl")
#Filip
df_init <- read_excel("C:/Users/filip/Downloads/econometrics & QF case study 2026/econometrics & QF case study 2026/data/data.xlsx")[-1, ]

#Extract Data from Excel, confirm correct columns, exclude first observation since NA

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
plot(df_init$`RV5_SS × 10^4`, 
     ylab = "RV5_SS x 10^4",
     main = "Realized Variance, daily 2009 - 2026")


plot(df_init$VIX, 
     ylab = "VIX price",
     main = "VIX")



# Pretty Graphs - > CC return, rv5, vix


plot(
  df_init$Date,
  df_init$`CC Return (%)`,
  col = ifelse(df_init$`CC Return (%)` < 0, "blue", "red"),
  pch = 20,
  xlab = "Date",
  ylab = "CC Return",
  main = "Close to Close Returns (daily 2009 - 2026, in %)"
)

abline(h = 0, lty = 2, col = "black")


threshold_75VIX <- quantile(df_init$VIX, 0.75, na.rm = TRUE)

print(threshold_75VIX)

plot(
  df_init$Date,
  df_init$VIX,
  col = ifelse(df_init$VIX < 21.3175, "black", "black"),
  pch = 20,
  xlab = "Date",
  ylab = "VIX",
  main = "VIX (daily 2009 - 2026, in %)"
)

abline(h = 21.3175, lty = 2, col = "black")

threshold_75rv <- quantile(df_init$`RV5_SS × 10^4`, 0.75, na.rm = TRUE)
print(threshold_75rv)

plot(
  df_init$Date,
  df_init$`RV5_SS × 10^4`,
  col = ifelse(df_init$`RV5_SS × 10^4` < 1.074468, "blue", "red"),
  pch = 20,
  xlab = "Date",
  ylab = "RV",
  main = "RV (daily 2009 - 2026, in %)"
)

abline(h = 1.074468, lty = 2, col = "black")

#Baseline summary statistics of key variables
summary(df_init$`CC Return (%)`)
summary(df_init$`RV5_SS × 10^4`)
summary(df_init$VIX)
print(cat(("Skewness and Kurtosis (respectively) for Close to Close returns"), skewness(df_init$`CC Return (%)`), kurtosis(df_init$`CC Return (%)`)))
print(cat(("Skewness and Kurtosis (respectively) for Realized variance * 10^4"), skewness(df_init$`RV5_SS × 10^4`), kurtosis(df_init$`RV5_SS × 10^4`)))
print(cat(("Skewness and Kurtosis (respectively) for VIX index"), skewness(df_init$VIX), kurtosis(df_init$VIX)))


#-----Discrete state space transition probabilities, complement Garch using changes in regimes. 

rt <- as.numeric(df_init$`CC Return (%)`)
rv <- as.numeric(df_init$`RV5_SS × 10^4`)
mu <- mean(rt)
T <- length(rt)

#lag-safe threshold -> top 25% realised variance = high-vol state
threshold <- quantile(rv, 0.75, na.rm = TRUE)

state <- ifelse(rv > threshold, "H", "L")

markov_df <- data.frame(
  state_t   = state[-length(state)],
  state_tp1 = state[-1],
  r_t       = rt[-length(rt)]
)


markov_df$return_sign <- ifelse(markov_df$r_t < 0, "negative", "positive_or_zero")

# Conditional transition probabilities:
tab <- xtabs(~ state_t + return_sign + state_tp1, data = markov_df)

prob <- prop.table(tab, margin = c(1, 2))

prob[, , "H"]
prob[, , "L"]

head(tab)


#__Markov Printing__
high_tbl <- round(prob[, , "H"], 3)
low_tbl  <- round(prob[, , "L"], 3)

colnames(high_tbl) <- c("Negative Return", "Positive / Zero Return")
rownames(high_tbl) <- c("Current State: High", "Current State: Low")

colnames(low_tbl) <- c("Negative Return", "Positive / Zero Return")
rownames(low_tbl) <- c("Current State: High", "Current State: Low")

kable(high_tbl,
      caption = "Probability of Transitioning to High-Volatility State")

kable(low_tbl,
      caption = "Probability of Transitioning to Low-Volatility State")




# ============================================================
# 5-STATE VOLATILITY REGIME EXTENSION
# ============================================================


# ----------------------------
# 1. Define variables
# ----------------------------

rv <- as.numeric(df_init$`RV5_SS × 10^4`)
r  <- as.numeric(df_init$`CC Return (%)`)

# Create 5 volatility states using quintiles of realised variance
q <- quantile(rv, probs = seq(0, 1, 0.2), na.rm = TRUE)

vol_state <- cut(
  rv,
  breaks = q,
  include.lowest = TRUE,
  labels = c("Very Low", "Low", "Medium", "High", "Extreme")
)

# Build transition dataset
markov5_df <- data.frame(
  state_t   = vol_state[-length(vol_state)],
  state_tp1 = vol_state[-1],
  r_t       = r[-length(r)]
)

markov5_df$return_sign <- ifelse(
  markov5_df$r_t < 0,
  "Negative",
  "Positive_or_Zero"
)

markov5_df$abs_return <- abs(markov5_df$r_t)

# Remove missing values
markov5_df <- na.omit(markov5_df)

# Make sure states are ordered correctly
state_levels <- c("Very Low", "Low", "Medium", "High", "Extreme")

markov5_df$state_t <- factor(markov5_df$state_t, levels = state_levels)
markov5_df$state_tp1 <- factor(markov5_df$state_tp1, levels = state_levels)
markov5_df$return_sign <- factor(markov5_df$return_sign)


# ============================================================
# 2. Transition matrices conditional on return sign
# ============================================================

tab5 <- xtabs(~ state_t + return_sign + state_tp1, data = markov5_df)

prob5 <- prop.table(tab5, margin = c(1, 2))

cat("\n====================================================\n")
cat(" 5-STATE CONDITIONAL TRANSITION PROBABILITIES\n")
cat("====================================================\n")

cat("\nP(S[t+1] = j | S[t] = i, Return is Negative)\n\n")
print(round(prob5[, "Negative", ], 3))

cat("\n----------------------------------------------------\n")

cat("\nP(S[t+1] = j | S[t] = i, Return is Positive or Zero)\n\n")
print(round(prob5[, "Positive_or_Zero", ], 3))

cat("\n====================================================\n")
cat(" RAW TRANSITION COUNTS\n")
cat("====================================================\n")

cat("\nCounts conditional on Negative returns:\n\n")
print(tab5[, "Negative", ])

cat("\nCounts conditional on Positive or Zero returns:\n\n")
print(tab5[, "Positive_or_Zero", ])



neg_matrix <- round(prob5[, "Negative", ], 3)
pos_matrix <- round(prob5[, "Positive_or_Zero", ], 3)

kable(
  neg_matrix,
  caption = "Transition probabilities conditional on negative returns"
)

kable(
  pos_matrix,
  caption = "Transition probabilities conditional on positive or zero returns"
)


# ============================================================
# 3. Multinomial logistic regression
# ============================================================

# Model: next volatility state depends on current state,
# return sign, and absolute return size

multi_logit <- multinom(
  state_tp1 ~ state_t + return_sign + abs_return,
  data = markov5_df,
  trace = FALSE
)

cat("\n====================================================\n")
cat(" MULTINOMIAL LOGISTIC REGRESSION RESULTS\n")
cat("====================================================\n\n")

summary_multi <- summary(multi_logit)

print(summary_multi)

# ============================================================
# 4. Approximate z-statistics and p-values
# ============================================================

coefs <- summary_multi$coefficients
ses   <- summary_multi$standard.errors

print(coefs)


z_vals <- coefs / ses
p_vals <- 2 * (1 - pnorm(abs(z_vals)))

cat("\n====================================================\n")
cat(" APPROXIMATE P-VALUES\n")
cat("====================================================\n\n")

print(round(p_vals, 4))


# ============================================================
# 5. Predicted probabilities example
# ============================================================

example_data <- data.frame(
  state_t = factor(
    c("Low", "Low"),
    levels = state_levels
  ),
  return_sign = factor(
    c("Negative", "Positive_or_Zero"),
    levels = levels(markov5_df$return_sign)
  ),
  abs_return = mean(markov5_df$abs_return, na.rm = TRUE)
)

pred_probs <- predict(
  multi_logit,
  newdata = example_data,
  type = "probs"
)

rownames(pred_probs) <- c(
  "Current Low + Negative Return",
  "Current Low + Positive/Zero Return"
)

cat("\n====================================================\n")
cat(" PREDICTED NEXT-STATE PROBABILITIES\n")
cat("====================================================\n\n")

print(round(pred_probs, 3))





diff_matrix <- prob5[, "Negative", ] - prob5[, "Positive_or_Zero", ]

diff_df <- melt(diff_matrix)
colnames(diff_df) <- c("Current_State", "Next_State", "Difference")


ggplot(diff_df, aes(x = Next_State, y = Current_State, fill = Difference)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Difference, 3)), size = 4) +
  scale_fill_gradient2(
    low = "lightblue",
    mid = "white",
    high ="red",
    midpoint = 0
  ) +
  labs(
    title = "Leverage Asymmetry in Volatility Regime Transitions",
    subtitle = "Difference: P(next state | negative return) - P(next state | non-negative return)",
    x = "Next Volatility State",
    y = "Current Volatility State",
    fill = "Difference"
  ) +
  theme_excel()

neg_df <- melt(prob5[, "Negative", ])
pos_df <- melt(prob5[, "Positive_or_Zero", ])

neg_df$return_sign <- "Negative return"
pos_df$return_sign <- "Positive / zero return"

plot_df <- rbind(neg_df, pos_df)
colnames(plot_df)[1:3] <- c("Current_State", "Next_State", "Probability")

ggplot(plot_df, aes(x = Next_State, y = Current_State, fill = Probability)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Probability, 2)), size = 3.5) +
  facet_wrap(~ return_sign) +
  labs(
    title = "Volatility Regime Transition Matrices by Return Sign",
    x = "Next Volatility State",
    y = "Current Volatility State",
    fill = "Probability"
  ) +
  theme_void()

# Probability of moving to High or Extreme next state
markov5_df$upper_next <- markov5_df$state_tp1 %in% c("High", "Extreme")

upper_prob <- aggregate(
  upper_next ~ state_t + return_sign,
  data = markov5_df,
  FUN = mean
)

ggplot(upper_prob, aes(x = state_t, y = upper_next, fill = return_sign)) +
  geom_col(position = "dodge") +
  labs(
    title = "Probability of Moving to High or Extreme Volatility",
    x = "Current Volatility State",
    y = "Probability",
    fill = "Return Sign"
  ) +
  theme_minimal()




pred_df <- as.data.frame(pred_probs)
pred_df$Scenario <- rownames(pred_probs)

pred_long <- reshape2::melt(pred_df, id.vars = "Scenario")
colnames(pred_long) <- c("Scenario", "Next_State", "Probability")

ggplot(pred_long, aes(x = Next_State, y = Probability, fill = Scenario)) +
  geom_col(position = "dodge") +
  labs(
    title = "Predicted Next-State Probabilities from Multinomial Logit",
    x = "Next Volatility State",
    y = "Predicted Probability"
  ) +
  theme_minimal()




# 5 group - > RV splits for regime-like transition probabilities 

# --------------------------
# RV plot with 5 variance groups
# --------------------------

rv <- as.numeric(df_init$`RV5_SS × 10^4`)

# Quintile breakpoints: 0%, 20%, 40%, 60%, 80%, 100%
rv_breaks <- quantile(rv, probs = seq(0, 1, 0.2), na.rm = TRUE)

state_labels <- c("Very Low", "Low", "Medium", "High", "Extreme")

rv_state <- cut(
  rv,
  breaks = rv_breaks,
  include.lowest = TRUE,
  labels = state_labels
)

# Print split definitions
cat("\n============================================\n")
cat(" REALISED VARIANCE STATE DEFINITIONS\n")
cat("============================================\n\n")

for (i in 1:length(state_labels)) {
  cat(
    state_labels[i], ": ",
    round(rv_breaks[i], 4), " to ",
    round(rv_breaks[i + 1], 4), "\n",
    sep = ""
  )
}

cat("\nObservations per state:\n")
print(table(rv_state))

# Plot
state_cols <- c(
  "Very Low" = "darkblue",
  "Low"      = "skyblue",
  "Medium"   = "gold",
  "High"     = "orange",
  "Extreme"  = "red"
)

plot(
  df_init$Date,
  rv,
  col = state_cols[as.character(rv_state)],
  pch = 20,
  xlab = "Date",
  ylab = "RV5_SS × 10^4",
  main = "Realised Variance Split into Five Volatility States"
)

# Add horizontal split lines
abline(h = rv_breaks[2:5], lty = 2, col = "black")

legend(
  "topright",
  legend = state_labels,
  col = state_cols[state_labels],
  pch = 20,
  title = "Volatility State",
  cex = 0.8
)



plot(
  df_init$Date,
  log(rv + 1e-6),
  col = state_cols[as.character(rv_state)],
  pch = 20,
  xlab = "Date",
  ylab = "log(RV5_SS × 10^4)",
  main = "Realised Variance States on Log Scale"
)

abline(h = log(rv_breaks[2:5] + 1e-6), lty = 2, col = "black")

legend(
  "topright",
  legend = state_labels,
  col = state_cols[state_labels],
  pch = 20,
  title = "Volatility State",
  cex = 0.8
)


#--------Regime Graphs -- end 



#--------------------------GARCH SETUP --------------------------



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
  h[1] <- omega / (1-alpha-beta)
  
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
  h[1] <- omega / (1 - (alpha1 + alpha2)/2 - beta)
  
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
  phi <- par[4]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || phi < 0 || beta + phi + alpha + 2*phi*alpha >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- omega / (1 - beta - phi - alpha - 2*phi*alpha)
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  
  
  
  # Generate conditional variances recursively
  for (t in 1:(T - 1)) {
    h[t + 1] <- 0.5*(omega + beta*h[t] + alpha*(rt[t] - mu)^2) + 0.5*sqrt((omega + beta*h[t] + alpha*(rt[t] - mu)^2)^2 + 4*phi*h[t]*(rt[t+1] - mu)^2)
  }
  
  # Now check h after it has been generated
  if (any(h <= 0) || any(is.na(h)) || any(is.infinite(h))) {
    return(1e10)
  }
  
  
  h_tm1 <- h[-T]
  h <- h[-1]
  rt <- rt[-1]
  
  # Negative log-likelihood
  RTgarch11ll <- sum(0.5*log(2*pi)+0.5*(rt-mu)^2/h-log(sqrt(h)/(h+phi*h_tm1*(rt-mu)^2/h)))
  
  
  if (!is.finite(RTgarch11ll)) {
    return(1e10)
  }
  
  return(RTgarch11ll)
}




start_par_RT <- c(
  omega = 0.04557794,
  alpha = 0.19471400,
  beta  = 0.77334072,
  phi = 0.01
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
  omega  = 0.04150796,
  alpha1 = 0.26965124,
  alpha2 = 0.01997420,
  beta   = 0.81654904,
  phi1   = 0.015,
  phi2   = 0.005
)


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

## Dummy Garch
Dummygarch11 <- function(par, rt, mu, DUM) {
  
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  delta <- par[4]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || alpha + delta*mean(DUM) + beta >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- omega / (1 - alpha - delta*mean(DUM) - beta)
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  # Generate conditional variances recursively
  for (t in 1:(T - 1)) {
    h[t + 1] <- omega + (alpha+delta*mean(DUM)) * (rt[t] - mu)^2 + beta * h[t]
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

start_par <- c(
  omega = 0.05,
  alpha = 0.10,
  beta  = 0.85,
  delta = 0.01
)

# Estimate model
Dummygarch_fit <- optim(
  par = start_par,
  fn = Dummygarch11,
  rt = rt,
  mu = mu,
  DUM = DummyAfterHoliday, 
  method = "L-BFGS-B",
  #Omega > 0 to inf, remaining are bounded by 0 and 1 
  lower = c(1e-8, 0, 0),
  upper = c(Inf, 1, 1)
)

#post - optimisation parameter estimates 
Dummygarch_fit$par
Dummygarch_fit$value
Dummygarch_fit$convergence


## Dummy GJR Garch

DummyGARCH_gjr <- function(par, rt, mu, DUM) {
  
  omega  <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  beta   <- par[4]
  delta <- par[5]
  
  # Parameter restrictions
  # omega > 0, alpha1 >= 0, alpha2 >= 0, beta >= 0
  # Stationarity: alpha1 + 0.5 * alpha2 + beta < 1
  if (omega <= 0 || alpha1 < 0 || alpha2 < 0 || beta < 0 ||
      (alpha1 + alpha2) / 2 + delta*mean(DUM) + beta >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance
  h[1] <- omega / (1 - (alpha1 + alpha2)/2 - delta * mean(DUM) - beta)
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  shock <- rt - mu
  
  for (t in 1:(T - 1)) {
    
    indicator <- ifelse(shock[t] < 0, 1, 0)
    
    h[t+1] <- omega +
      (alpha1 + delta*DUM[t]) * indicator * shock[t]^2 +        # alpha1 when negative
      (alpha2 + delta*DUM[t]) * (1 - indicator) * shock[t]^2 +  # alpha2 when positive
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
  beta   = 0.85,
  delta = 0.01
)

Dummygjr_fit <- optim(
  par = start_par_gjr,
  fn = DummyGARCH_gjr,
  rt = rt,
  mu = mu,
  DUM = DummyAfterHoliday,
  method = "L-BFGS-B",
  lower = c(1e-8, 0, 0, 0),
  upper = c(10, 1, 1, 1)
)

Dummygjr_fit$par
Dummygjr_fit$value
Dummygjr_fit$convergence


## Dummy RT GARCH

DummyRTgarch11 <- function(par, rt, mu, DUM) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  phi <- par[4]
  delta <- par[5]
  
  # Penalize invalid parameter values
  if (omega <= 0 || alpha < 0 || beta < 0 || phi < 0 || beta + phi + alpha + delta*mean(DUM) + 2*phi*(alpha + delta*mean(DUM)) >= 1) {
    return(1e10)
  }
  
  T <- length(rt)
  h <- numeric(T)
  
  # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
  h[1] <- omega / (1 - beta - phi - alpha - delta*mean(DUM) - 2*phi*(alpha + delta*mean(DUM)))
  
  if (!is.finite(h[1]) || h[1] <= 0) {
    return(1e10)
  }
  
  # Generate conditional variances recursively
  for (t in 1:(T - 1)) {
    h[t + 1] <- 0.5*(omega + beta*h[t] + (alpha+delta*DUM[t+1])*(rt[t] - mu)^2) + 0.5*sqrt((omega + beta*h[t] + (alpha+delta*DUM[t+1])*(rt[t] - mu)^2)^2 + 4*(phi)*h[t]*(rt[t+1] - mu)^2)
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
  RTgarch11ll <- sum(0.5*log(2*pi)+0.5*(rt-mu)^2/h-log(sqrt(h)/(h+(phi)*h_tm1*(rt-mu)^2/h)))
  
  
  if (!is.finite(RTgarch11ll)) {
    return(1e10)
  }
  
  return(RTgarch11ll)
}


start_par_RT <- c(
  omega = 0.04590705,
  alpha = 0.19226977 ,
  beta  = 0.77380626 ,
  phi = 0.01,
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

DummyRTgarch_fit$par
DummyRTgarch_fit$value
DummyRTgarch_fit$convergence


## Dummy RT GJR garch

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


DummyRTgarchGJR_fit$par
DummyRTgarchGJR_fit$value
DummyRTgarchGJR_fit$convergence


#--------------------------HESSIANS--------------------------------------------

rt <- df_init$`CC Return (%)`
mu <- mean(rt)
T <- length(rt)



## Garch model

FischerGARCH <- matrix(0,3,3)
for (i in 1:T){
  
  ObsLL <- function(par) {
    
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
    ll <- 0.5 * (-log(2 * pi) - log(h[i]) - ((rt[i] - mu)^2 / h[i]))
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, garch_fit$par, method = "complex"), grad(ObsLL, garch_fit$par, method = "complex"))
}

CovGARCH <- solve(FischerGARCH)
tStats_GARCh <- c(omega = garch_fit$par[1]/sqrt(CovGARCH[1,1]),
                  alpha = garch_fit$par[2]/sqrt(CovGARCH[2,2]),
                  beta = garch_fit$par[3]/sqrt(CovGARCH[3,3])
)
print(tStats_GARCh)

## GJR-Garch model

FischerGARCH <- matrix(0,4,4)
for (i in 1:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    beta  <- par[4]
    
    h <- numeric(i)
    
    # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
    h[1] <- omega/(1-(alpha1+alpha2)/2-beta)
    
    shock <- rt - mu
    
    if (i > 1){
      for (t in 2:i){
        indicator <- ifelse(shock[t-1] < 0, 1, 0)
        
        h[t] <- omega +
          alpha1 * indicator * shock[t-1]^2 +        # alpha1 when negative
          alpha2 * (1 - indicator) * shock[t-1]^2 +  # alpha2 when positive
          beta * h[t-1]
      }
    }
    
    # log-likelihood
    ll <- 0.5 * (-log(2 * pi) - log(h[i]) - ((shock[i])^2 / h[i]))
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, gjr_fit$par, method = "complex"), grad(ObsLL, gjr_fit$par, method = "complex"))
}

CovGJRGARCH <- solve(FischerGARCH)
tStats_GJRGARCh <- c(omega = gjr_fit$par[1]/sqrt(CovGJRGARCH[1,1]),
                     alpha1 = gjr_fit$par[2]/sqrt(CovGJRGARCH[2,2]),
                     alpha2 = gjr_fit$par[3]/sqrt(CovGJRGARCH[3,3]),
                     beta = gjr_fit$par[4]/sqrt(CovGJRGARCH[4,4])
)
print(tStats_GJRGARCh)

## RTGarch model

FischerGARCH <- matrix(0,4,4)
for (i in 2:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha <- par[2]
    beta  <- par[3]
    phi <- par[4]
    
    
    h <- numeric(i)
    
    # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
    h[1] <- omega / (1 - beta - phi)
    
    shock <- rt - mu
    
    for (t in 2:i){
      h[t] <- 0.5*(omega + beta*h[t-1] + alpha*(shock[t-1])^2) + 0.5*sqrt((omega + beta*h[t-1] + alpha*(shock[t-1])^2)^2 + 4*phi*h[t-1]*(shock[t])^2)
    }
    
    # log-likelihood
    ll <- -0.5*log(2*pi)-0.5*(shock[i])^2/h[i]+log(sqrt(h[i])/(h[i]+phi*h[i-1]*(shock[i])^2/h[i]))
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, RTgarch_fit$par, method = "complex"), grad(ObsLL, RTgarch_fit$par, method = "complex"))
}

CovRTGARCH <- solve(FischerGARCH)
tStats_RTGARCH <- c(omega = RTgarch_fit$par[1]/sqrt(CovRTGARCH[1,1]),
                    alpha = RTgarch_fit$par[2]/sqrt(CovRTGARCH[2,2]),
                    beta = RTgarch_fit$par[3]/sqrt(CovRTGARCH[3,3]),
                    phi = RTgarch_fit$par[4]/sqrt(CovRTGARCH[4,4])
)
print(tStats_RTGARCH)

## RTGJRGarch model

FischerGARCH <- matrix(0,6,6)
for (i in 2:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    beta <- par[4]
    phi1 <- par[5]
    phi2 <- par[6]
    
    phi_bar   <- (phi1 + phi2) / 2
    alpha_bar <- (alpha1 + alpha2) / 2
    
    h <- numeric(i)
    h[1] <- omega / (1 - beta - phi_bar)
    
    
    for (t in 2:i) {
      alpha_t_1 <- ifelse(rt[t-1]   <= mu, alpha1, alpha2)
      phi_t   <- ifelse(rt[t] <= mu, phi1,   phi2)
      g_t_1     <- omega + beta*h[t-1] + alpha_t_1*(rt[t-1]-mu)^2
      h[t]  <- 0.5*g_t_1 + 0.5*sqrt(g_t_1^2 + 4*phi_t*h[t-1]*(rt[t]-mu)^2)
    }
    
    phi_i   <- ifelse(rt[i] <= mu, phi1,   phi2)
    
    ll <- -0.5*log(2*pi) - 0.5*(rt[i]-mu)^2/h[i] + log(sqrt(h[i])/(h[i] + phi_i*h[i-1]*(rt[i]-mu)^2/h[i]))
    
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, RTgjr_fit$par, method = "complex"), grad(ObsLL, RTgjr_fit$par, method = "complex"))
}

CovRTGJRGARCH <- solve(FischerGARCH)
tStats_RTGJRGARCH <- c(omega = RTgjr_fit$par[1]/sqrt(CovRTGJRGARCH[1,1]),
                       alpha1 = RTgjr_fit$par[2]/sqrt(CovRTGJRGARCH[2,2]),
                       alpha2 = RTgjr_fit$par[3]/sqrt(CovRTGJRGARCH[3,3]),
                       beta = RTgjr_fit$par[4]/sqrt(CovRTGJRGARCH[4,4]),
                       phi1 = RTgjr_fit$par[5]/sqrt(CovRTGJRGARCH[5,5]),
                       phi2 = RTgjr_fit$par[6]/sqrt(CovRTGJRGARCH[6,6])
)
print(tStats_RTGJRGARCH)


## ----------------- Dummmy Garch model

DUM <- DummyAfterHoliday


FischerGARCH <- matrix(0,4,4)
for (i in 1:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha <- par[2]
    beta  <- par[3]
    delta <- par[4]
    
    h <- numeric(i)
    
    # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
    h[1] <- omega/(1 - alpha - delta*mean(DUM) - beta)
    
    
    
    
    if (i > 1){
      for (t in 2:i){
        h[t] <- omega + (alpha + delta*DUM[t-1]) * (rt[t-1] - mu)^2 + beta * h[t-1]
      }
    }
    
    # log-likelihood
    ll <- 0.5 * (-log(2 * pi) - log(h[i]) - ((rt[i] - mu)^2 / h[i]))
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, Dummygarch_fit$par, method = "complex"), grad(ObsLL, Dummygarch_fit$par, method = "complex"))
}

CovDUMMYGARCH <- solve(FischerGARCH)
tStats_DUMMYGARCh <- c(omega = Dummygarch_fit$par[1]/sqrt(CovDUMMYGARCH[1,1]),
                       alpha = Dummygarch_fit$par[2]/sqrt(CovDUMMYGARCH[2,2]),
                       beta = Dummygarch_fit$par[3]/sqrt(CovDUMMYGARCH[3,3]),
                       delta = Dummygarch_fit$par[4]/sqrt(CovDUMMYGARCH[4,4])
)
print(tStats_DUMMYGARCh)

## ----------------- Dummmy GJR-Garch model

FischerGARCH <- matrix(0,5,5)
for (i in 1:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    beta  <- par[4]
    delta <- par[5]
    
    h <- numeric(i)
    
    # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
    h[1] <- omega/(1-(alpha1+alpha2)/2 - delta*mean(DUM) -beta)
    
    shock <- rt - mu
    
    if (i > 1){
      for (t in 2:i){
        indicator <- ifelse(shock[t-1] < 0, 1, 0)
        
        h[t] <- omega +
          (alpha1 + delta*DUM[t-1]) * indicator * shock[t-1]^2 +        # alpha1 when negative
          (alpha2 + delta*DUM[t-1]) * (1 - indicator) * shock[t-1]^2 +  # alpha2 when positive
          beta * h[t-1]
      }
    }
    
    # log-likelihood
    ll <- 0.5 * (-log(2 * pi) - log(h[i]) - ((shock[i])^2 / h[i]))
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, Dummygjr_fit$par, method = "complex"), grad(ObsLL, Dummygjr_fit$par, method = "complex"))
}

CovDUMMYGJRGARCH <- solve(FischerGARCH)
tStats_DUMMYGJRGARCh <- c(omega = Dummygjr_fit$par[1]/sqrt(CovGJRGARCH[1,1]),
                          alpha1 = Dummygjr_fit$par[2]/sqrt(CovGJRGARCH[2,2]),
                          alpha2 = Dummygjr_fit$par[3]/sqrt(CovGJRGARCH[3,3]),
                          beta = Dummygjr_fit$par[4]/sqrt(CovGJRGARCH[4,4]),
                          delta = Dummygjr_fit$par[5]/sqrt(CovDUMMYGJRGARCH[5,5])
)
print(tStats_DUMMYGJRGARCh)

## ----------------- Dummmy RTGarch model

FischerGARCH <- matrix(0,5,5)
for (i in 2:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha <- par[2]
    beta  <- par[3]
    phi <- par[4]
    delta <- par[5]
    
    
    h <- numeric(i)
    
    # Initial conditional variance, we do not require h1 = hhat since it is unknown pre-computation (and doesnt affect convergence)
    h[1] <- omega / (1 - beta - phi - alpha - delta*mean(DUM) - 2*phi*(alpha + delta*mean(DUM)))
    
    shock <- rt - mu
    
    for (t in 2:i){
      h[t] <- 0.5*(omega + beta*h[t-1] + (alpha+delta*DUM[t-1])*(shock[t-1])^2) + 0.5*sqrt((omega + beta*h[t-1] + (alpha+delta*DUM[t-1])*(shock[t-1])^2)^2 + 4*phi*h[t-1]*(shock[t])^2)
    }
    
    # log-likelihood
    ll <- -0.5*log(2*pi)-0.5*(shock[i])^2/h[i]+log(sqrt(h[i])/(h[i]+phi*h[i-1]*(shock[i])^2/h[i]))
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, DummyRTgarch_fit$par, method = "complex"), grad(ObsLL, DummyRTgarch_fit$par, method = "complex"))
}

CovDUMMYRTGARCH <- solve(FischerGARCH)
tStats_DUMMYRTGARCH <- c(omega = DummyRTgarch_fit$par[1]/sqrt(CovDUMMYRTGARCH[1,1]),
                         alpha = DummyRTgarch_fit$par[2]/sqrt(CovDUMMYRTGARCH[2,2]),
                         beta = DummyRTgarch_fit$par[3]/sqrt(CovDUMMYRTGARCH[3,3]),
                         phi = DummyRTgarch_fit$par[4]/sqrt(CovDUMMYRTGARCH[4,4]),
                         delta = DummyRTgarch_fit$par[5]/sqrt(CovDUMMYRTGARCH[5,5])
)
print(tStats_DUMMYRTGARCH)

## ----------------- Dummmy RTGJRGarch model

FischerGARCH <- matrix(0,7,7)
for (i in 2:T){
  
  ObsLL <- function(par) {
    
    omega <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    beta <- par[4]
    phi1 <- par[5]
    phi2 <- par[6]
    delta <- par[7]
    
    phi_bar   <- (phi1 + phi2) / 2
    alpha_bar <- (alpha1 + alpha2) / 2
    
    h <- numeric(i)
    h[1] <- omega / (1 - beta - phi_bar - alpha_bar - delta*mean(DUM) + (alpha_bar + delta*mean(DUM))*phi_bar - 3/2*(alpha1*phi1 + alpha2*phi2) - 3*delta*mean(DUM)*phi_bar)
    
    
    for (t in 2:i) {
      alpha_t_1 <- ifelse(rt[t-1]   <= mu, alpha1, alpha2)
      phi_t   <- ifelse(rt[t] <= mu, phi1,   phi2)
      g_t_1     <- omega + beta*h[t-1] + (alpha_t_1 + delta*mean(DUM))*(rt[t-1]-mu)^2
      h[t]  <- 0.5*g_t_1 + 0.5*sqrt(g_t_1^2 + 4*phi_t*h[t-1]*(rt[t]-mu)^2)
    }
    
    phi_i   <- ifelse(rt[i] <= mu, phi1,   phi2)
    
    ll <- -0.5*log(2*pi) - 0.5*(rt[i]-mu)^2/h[i] + log(sqrt(h[i])/(h[i] + phi_i*h[i-1]*(rt[i]-mu)^2/h[i]))
    
    
    return(ll)
  }
  if (i %% 500 == 0) {
    print(i)
  }
  
  FischerGARCH <- FischerGARCH + outer(grad(ObsLL, DummyRTgarchGJR_fit$par, method = "complex"), grad(ObsLL, DummyRTgarchGJR_fit$par, method = "complex"))
}

CovDUMMYRTGJRGARCH <- solve(FischerGARCH)
tStats_DUMMYRTGJRGARCH <- c(omega = DummyRTgarchGJR_fit$par[1]/sqrt(CovDUMMYRTGJRGARCH[1,1]),
                            alpha1 = DummyRTgarchGJR_fit$par[2]/sqrt(CovDUMMYRTGJRGARCH[2,2]),
                            alpha2 = DummyRTgarchGJR_fit$par[3]/sqrt(CovDUMMYRTGJRGARCH[3,3]),
                            beta = DummyRTgarchGJR_fit$par[4]/sqrt(CovDUMMYRTGJRGARCH[4,4]),
                            phi1 = DummyRTgarchGJR_fit$par[5]/sqrt(CovDUMMYRTGJRGARCH[5,5]),
                            phi2 = DummyRTgarchGJR_fit$par[6]/sqrt(CovDUMMYRTGJRGARCH[6,6]),
                            delta = DummyRTgarchGJR_fit$par[7]/sqrt(CovDUMMYRTGJRGARCH[7,7])
)
print(tStats_RTGJRGARCH)





#--------------------------h_t PLOTTING ---------------------------------------

#GARCH recursion



omega <- garch_fit$par[1]
alpha <- garch_fit$par[2]
beta <- garch_fit$par[3]

h_1 <- numeric(T)
h_1[1] <- omega / (1-alpha-beta)

for (t in 1:(T - 1)) {
  h_1[t + 1] <- omega + alpha*(rt[t] - mu)^2 + beta*h_1[t]
}


#GARCH-GJR recursion


omega <- gjr_fit$par[1]
alpha1 <- gjr_fit$par[2]
alpha2 <- gjr_fit$par[3]
beta <- gjr_fit$par[4]

h_garch <- numeric(T)
h_garch[1] <- omega / (1 - (alpha1+alpha2)/2 - beta) 

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
phi <- RTgarch_fit$par[4]

h_htgarch[1] <- omega / (1 - beta - phi)

for (t in 1:(T - 1)) {
  h_htgarch[t + 1] <- 0.5*(omega + beta*h_htgarch[t] + alpha*(rt[t] - mu)^2) + 0.5*sqrt((omega + beta*h_htgarch[t] + alpha*(rt[t] - mu)^2)^2 + 4*phi*h_htgarch[t]*(rt[t+1] - mu)^2)
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

#Alternative post-covid split
#df_init$Date <- as.Date(df_init$Date) 
#T_est <- which(df_init$Date == as.Date("2022-12-31"))
#if(length(T_est) == 0) T_est <- max(which(df_init$Date <= as.Date("2022-12-31")))
#rt_est <- rt[1:T_est]
#mu_est <- mean(rt_est)

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


omega <- garch_est_fit$par[1]
alpha <- garch_est_fit$par[2]
beta <- garch_est_fit$par[3]

h_1 <- numeric(T)
h_1[1] <- omega / (1 - alpha - beta)


for (t in 1:(T - 1)) {
  h_1[t + 1] <- omega + alpha*(rt[t] - mu_est)^2 + beta*h_1[t]
}


#GARCH-GJR recursion

omega <- gjr_est_fit$par[1]
alpha1 <- gjr_est_fit$par[2]
alpha2 <- gjr_est_fit$par[3]
beta <- gjr_est_fit$par[4]

h_garch <- numeric(T)
h_garch[1] <- omega / (1 - (alpha1+alpha2)/2 - beta)

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
phi <- RTgarch_est_fit$par[4]

h_htgarch[1] <- omega / (1 - beta - phi)

for (t in 1:(T - 1)) {
  h_htgarch[t + 1] <- 0.5*(omega + beta*h_htgarch[t] + alpha*(rt[t] - mu_est)^2) + 0.5*sqrt((omega + beta*h_htgarch[t] + alpha*(rt[t] - mu_est)^2)^2 + 4*phi*h_htgarch[t]*(rt[t+1] - mu_est)^2)
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
phi <- RTgarch_est_fit$par[4]
indx <- 1

for (t in T_est:(T - 1)) {
  PV_RTGARCH[indx] <- g_htgarch[t] + h_htgarch[t] * kurtosis * phi
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

#-----------PLOT FOR PVs and TVs-------------------------------------



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

#---------------DIEBOLD-MARIANO & MINCER-ZARNOWITZ TESTS-----------------------------------


crisis_1 <- 1:760
calm <- 761:1211
crisis_2 <- 1212:1597

#####################################
# MSE LOSSES
####################################

######RV TARGETS#########

# GJR-GARCH vs RT-GJR-GARCH
dm_RV_c1_gjr_rt  <- dm.test((TV_rv[crisis_1] - PV_GARCHGJR[crisis_1]), 
                            (TV_rv[crisis_1] - PV_RTGARCHGJR[crisis_1]), alternative = "two.sided")

# RT-GJR-GARCH vs HAR-RV Baseline
dm_RV_c1_rt_har  <- dm.test((TV_rv[crisis_1] - PV_RTGARCHGJR[crisis_1]), 
                            (TV_rv[crisis_1] - PV_HAR_RV[crisis_1]), alternative = "two.sided")

# RT-GJR-GARCH vs VIX Baseline
dm_RV_c1_rt_vix  <- dm.test((TV_rv[crisis_1] - PV_RTGARCHGJR[crisis_1]), 
                            (TV_rv[crisis_1] - PV_VIX[crisis_1]), alternative = "two.sided")



#calm (2023 to september 2024)

dm_RV_calm_gjr_rt <- dm.test((TV_rv[calm] - PV_GARCHGJR[calm]), 
                             (TV_rv[calm] - PV_RTGARCHGJR[calm]), alternative = "two.sided")

dm_RV_calm_rt_har <- dm.test((TV_rv[calm] - PV_RTGARCHGJR[calm]), 
                             (TV_rv[calm] - PV_HAR_RV[calm]), alternative = "two.sided")

dm_RV_calm_rt_vix <- dm.test((TV_rv[calm] - PV_RTGARCHGJR[calm]), 
                             (TV_rv[calm] - PV_VIX[calm]), alternative = "two.sided")



#Crisis 2 (starting october 2024)
dm_RV_c2_gjr_rt  <- dm.test((TV_rv[crisis_2] - PV_GARCHGJR[crisis_2]), 
                            (TV_rv[crisis_2] - PV_RTGARCHGJR[crisis_2]), alternative = "two.sided")

dm_RV_c2_rt_har  <- dm.test((TV_rv[crisis_2] - PV_RTGARCHGJR[crisis_2]), 
                            (TV_rv[crisis_2] - PV_HAR_RV[crisis_2]), alternative = "two.sided")

dm_RV_c2_rt_vix  <- dm.test((TV_rv[crisis_2] - PV_RTGARCHGJR[crisis_2]), 
                            (TV_rv[crisis_2] - PV_VIX[crisis_2]), alternative = "two.sided")

#entire evaluation window
dm_RV_gjr_rt  <- dm.test((TV_rv - PV_GARCHGJR), 
                         (TV_rv - PV_RTGARCHGJR), alternative = "two.sided")

dm_RV_rt_har  <- dm.test((TV_rv - PV_RTGARCHGJR), 
                         (TV_rv - PV_HAR_RV), alternative = "two.sided")

dm_RV_rt_vix  <- dm.test((TV_rv - PV_RTGARCHGJR), 
                         (TV_rv - PV_VIX), alternative = "two.sided")

#Results
print("ENTIRE WINDOW")
print(dm_RV_gjr_rt)
print(dm_RV_rt_har)
print(dm_RV_rt_vix)

print("WINDOW 1")
print(dm_RV_c1_gjr_rt)
print(dm_RV_c1_rt_har)
print(dm_RV_c1_rt_vix)

print("WINDOW 2")
print(dm_RV_calm_gjr_rt)
print(dm_RV_calm_rt_har)
print(dm_RV_calm_rt_vix)

print("WINDOW 3")
print(dm_RV_c2_gjr_rt)
print(dm_RV_c2_rt_har)
print(dm_RV_c2_rt_vix)


#R TARGETS##########

# GJR-GARCH vs RT-GJR-GARCH
dm_R_c1_gjr_rt  <- dm.test((TV_r[crisis_1] - PV_GARCHGJR[crisis_1]), 
                           (TV_r[crisis_1] - PV_RTGARCHGJR[crisis_1]), alternative = "two.sided")

# RT-GJR-GARCH vs HAR-RV Baseline
dm_R_c1_rt_har  <- dm.test((TV_r[crisis_1] - PV_RTGARCHGJR[crisis_1]), 
                           (TV_r[crisis_1] - PV_HAR_RV[crisis_1]), alternative = "two.sided")

# RT-GJR-GARCH vs VIX Baseline
dm_R_c1_rt_vix  <- dm.test((TV_r[crisis_1] - PV_RTGARCHGJR[crisis_1]), 
                           (TV_r[crisis_1] - PV_VIX[crisis_1]), alternative = "two.sided")



#calm (2023 to september 2024)

dm_R_calm_gjr_rt <- dm.test((TV_r[calm] - PV_GARCHGJR[calm]), 
                            (TV_r[calm] - PV_RTGARCHGJR[calm]), alternative = "two.sided")

dm_R_calm_rt_har <- dm.test((TV_r[calm] - PV_RTGARCHGJR[calm]), 
                            (TV_r[calm] - PV_HAR_RV[calm]), alternative = "two.sided")

dm_R_calm_rt_vix <- dm.test((TV_r[calm] - PV_RTGARCHGJR[calm]), 
                            (TV_r[calm] - PV_VIX[calm]), alternative = "two.sided")



#Crisis 2 (starting october 2024)
dm_R_c2_gjr_rt  <- dm.test((TV_r[crisis_2] - PV_GARCHGJR[crisis_2]), 
                           (TV_r[crisis_2] - PV_RTGARCHGJR[crisis_2]), alternative = "two.sided")

dm_R_c2_rt_har  <- dm.test((TV_r[crisis_2] - PV_RTGARCHGJR[crisis_2]), 
                           (TV_r[crisis_2] - PV_HAR_RV[crisis_2]), alternative = "two.sided")

dm_R_c2_rt_vix  <- dm.test((TV_r[crisis_2] - PV_RTGARCHGJR[crisis_2]), 
                           (TV_r[crisis_2] - PV_VIX[crisis_2]), alternative = "two.sided")

#entire evaluation window
dm_R_gjr_rt  <- dm.test((TV_r - PV_GARCHGJR), 
                        (TV_r - PV_RTGARCHGJR), alternative = "two.sided")

dm_R_rt_har  <- dm.test((TV_r - PV_RTGARCHGJR), 
                        (TV_r - PV_HAR_RV), alternative = "two.sided")

dm_R_rt_vix  <- dm.test((TV_r - PV_RTGARCHGJR), 
                        (TV_r - PV_VIX), alternative = "two.sided")

#Results
print("ENTIRE WINDOW")
print(dm_R_gjr_rt)
print(dm_R_rt_har)
print(dm_R_rt_vix)

print("WINDOW 1")
print(dm_R_c1_gjr_rt)
print(dm_R_c1_rt_har)
print(dm_R_c1_rt_vix)

print("WINDOW 2")
print(dm_R_calm_gjr_rt)
print(dm_R_calm_rt_har)
print(dm_R_calm_rt_vix)

print("WINDOW 3")
print(dm_R_c2_gjr_rt)
print(dm_R_c2_rt_har)
print(dm_R_c2_rt_vix)


#####################################
# QLIKE LOSSES
####################################

qlike_loss <- function(TV, PV) {
  log(PV) + TV / PV
}
#RV TARGETS###############
#crisis 1
# GJR-GARCH vs RT-GJR-GARCH
qdm_RV_c1_gjr_rt  <- dm.test(qlike_loss(TV_rv[crisis_1], PV_GARCHGJR[crisis_1]), 
                             qlike_loss(TV_rv[crisis_1] , PV_RTGARCHGJR[crisis_1]), alternative = "two.sided", h = 1)

# RT-GJR-GARCH vs HAR-RV Baseline
qdm_RV_c1_rt_har  <- dm.test(qlike_loss(TV_rv[crisis_1] , PV_RTGARCHGJR[crisis_1]), 
                             qlike_loss(TV_rv[crisis_1] , PV_HAR_RV[crisis_1]), alternative = "two.sided", h = 1)

# RT-GJR-GARCH vs VIX Baseline
qdm_RV_c1_rt_vix  <- dm.test(qlike_loss(TV_rv[crisis_1] , PV_RTGARCHGJR[crisis_1]), 
                             qlike_loss(TV_rv[crisis_1] , PV_VIX[crisis_1]), alternative = "two.sided", h = 1)



#calm (2023 to september 2024)

qdm_RV_calm_gjr_rt <- dm.test(qlike_loss(TV_rv[calm] , PV_GARCHGJR[calm]), 
                              qlike_loss(TV_rv[calm] , PV_RTGARCHGJR[calm]), alternative = "two.sided", h = 1)

qdm_RV_calm_rt_har <- dm.test(qlike_loss(TV_rv[calm] , PV_RTGARCHGJR[calm]), 
                              qlike_loss(TV_rv[calm] , PV_HAR_RV[calm]), alternative = "two.sided", h = 1)

qdm_RV_calm_rt_vix <- dm.test(qlike_loss(TV_rv[calm] , PV_RTGARCHGJR[calm]), 
                              qlike_loss(TV_rv[calm] , PV_VIX[calm]), alternative = "two.sided", h = 1)



#Crisis 2 (starting october 2024)
qdm_RV_c2_gjr_rt  <- dm.test(qlike_loss(TV_rv[crisis_2] , PV_GARCHGJR[crisis_2]), 
                             qlike_loss(TV_rv[crisis_2] , PV_RTGARCHGJR[crisis_2]), alternative = "two.sided", h = 1)

qdm_RV_c2_rt_har  <- dm.test(qlike_loss(TV_rv[crisis_2] , PV_RTGARCHGJR[crisis_2]), 
                             qlike_loss(TV_rv[crisis_2] , PV_HAR_RV[crisis_2]), alternative = "two.sided", h = 1)

qdm_RV_c2_rt_vix  <- dm.test(qlike_loss(TV_rv[crisis_2] , PV_RTGARCHGJR[crisis_2]), 
                             qlike_loss(TV_rv[crisis_2] , PV_VIX[crisis_2]), alternative = "two.sided", h = 1)

#entire evaluation window
qdm_RV_gjr_rt  <- dm.test(qlike_loss(TV_rv , PV_GARCHGJR), 
                          qlike_loss(TV_rv , PV_RTGARCHGJR), alternative = "two.sided", h = 1)

qdm_RV_rt_har  <- dm.test(qlike_loss(TV_rv , PV_RTGARCHGJR), 
                          qlike_loss(TV_rv , PV_HAR_RV), alternative = "two.sided", h = 1)

qdm_RV_rt_vix  <- dm.test(qlike_loss(TV_rv , PV_RTGARCHGJR), 
                          qlike_loss(TV_rv , PV_VIX), alternative = "two.sided", h = 1)

#Results
print("ENTIRE WINDOW")
print(qdm_RV_gjr_rt)
print(qdm_RV_rt_har)
print(qdm_RV_rt_vix)

print("WINDOW 1")
print(qdm_RV_c1_gjr_rt)
print(qdm_RV_c1_rt_har)
print(qdm_RV_c1_rt_vix)

print("WINDOW 2")
print(qdm_RV_calm_gjr_rt)
print(qdm_RV_calm_rt_har)
print(qdm_RV_calm_rt_vix)

print("WINDOW 3")
print(qdm_RV_c2_gjr_rt)
print(qdm_RV_c2_rt_har)
print(qdm_RV_c2_rt_vix)


#R TARGETS#################

#crisis 1
# GJR-GARCH vs RT-GJR-GARCH
qdm_R_c1_gjr_rt  <- dm.test(qlike_loss(TV_r[crisis_1], PV_GARCHGJR[crisis_1]), 
                            qlike_loss(TV_r[crisis_1] , PV_RTGARCHGJR[crisis_1]), alternative = "two.sided", h = 1)

# RT-GJR-GARCH vs HAR-RV Baseline
qdm_R_c1_rt_har  <- dm.test(qlike_loss(TV_r[crisis_1] , PV_RTGARCHGJR[crisis_1]), 
                            qlike_loss(TV_r[crisis_1] , PV_HAR_RV[crisis_1]), alternative = "two.sided", h = 1)

# RT-GJR-GARCH vs VIX Baseline
qdm_R_c1_rt_vix  <- dm.test(qlike_loss(TV_r[crisis_1] , PV_RTGARCHGJR[crisis_1]), 
                            qlike_loss(TV_r[crisis_1] , PV_VIX[crisis_1]), alternative = "two.sided", h = 1)



#calm (2023 to september 2024)

qdm_R_calm_gjr_rt <- dm.test(qlike_loss(TV_r[calm] , PV_GARCHGJR[calm]), 
                             qlike_loss(TV_r[calm] , PV_RTGARCHGJR[calm]), alternative = "two.sided", h = 1)

qdm_R_calm_rt_har <- dm.test(qlike_loss(TV_r[calm] , PV_RTGARCHGJR[calm]), 
                             qlike_loss(TV_r[calm] , PV_HAR_RV[calm]), alternative = "two.sided", h = 1)

qdm_R_calm_rt_vix <- dm.test(qlike_loss(TV_r[calm] , PV_RTGARCHGJR[calm]), 
                             qlike_loss(TV_r[calm] , PV_VIX[calm]), alternative = "two.sided", h = 1)



#Crisis 2 (starting october 2024)
qdm_R_c2_gjr_rt  <- dm.test(qlike_loss(TV_r[crisis_2] , PV_GARCHGJR[crisis_2]), 
                            qlike_loss(TV_r[crisis_2] , PV_RTGARCHGJR[crisis_2]), alternative = "two.sided", h = 1)

qdm_R_c2_rt_har  <- dm.test(qlike_loss(TV_r[crisis_2] , PV_RTGARCHGJR[crisis_2]), 
                            qlike_loss(TV_r[crisis_2] , PV_HAR_RV[crisis_2]), alternative = "two.sided", h = 1)

qdm_R_c2_rt_vix  <- dm.test(qlike_loss(TV_r[crisis_2] , PV_RTGARCHGJR[crisis_2]), 
                            qlike_loss(TV_r[crisis_2] , PV_VIX[crisis_2]), alternative = "two.sided", h = 1)

#entire evaluation window
qdm_R_gjr_rt  <- dm.test(qlike_loss(TV_r , PV_GARCHGJR), 
                         qlike_loss(TV_r , PV_RTGARCHGJR), alternative = "two.sided", h = 1)

qdm_R_rt_har  <- dm.test(qlike_loss(TV_r , PV_RTGARCHGJR), 
                         qlike_loss(TV_r , PV_HAR_RV), alternative = "two.sided", h = 1)

qdm_R_rt_vix  <- dm.test(qlike_loss(TV_r , PV_RTGARCHGJR), 
                         qlike_loss(TV_r , PV_VIX), alternative = "two.sided", h = 1)

#Results
print("ENTIRE WINDOW")
print(qdm_R_gjr_rt)
print(qdm_R_rt_har)
print(qdm_R_rt_vix)

print("WINDOW 1")
print(qdm_R_c1_gjr_rt)
print(qdm_R_c1_rt_har)
print(qdm_R_c1_rt_vix)

print("WINDOW 2")
print(qdm_R_calm_gjr_rt)
print(qdm_R_calm_rt_har)
print(qdm_R_calm_rt_vix)

print("WINDOW 3")
print(qdm_R_c2_gjr_rt)
print(qdm_R_c2_rt_har)
print(qdm_R_c2_rt_vix)

#---Mincer-Zarnowitz tests----

install.packages("car")
library(sandwich)
library(car)

##GARCH-GJR##

# rv

mz_RV_GARCHGJR <- lm(TV_rv ~ PV_GARCHGJR)


linearHypothesis(mz_RV_GARCHGJR, c("(Intercept) = 0", "PV_GARCHGJR = 1"), 
                 vcov = NeweyWest(mz_RV_GARCHGJR, lag = 1))
summary(mz_RV_GARCHGJR)

# r

mz_R_GARCHGJR <- lm(TV_r ~ PV_GARCHGJR)


linearHypothesis(mz_R_GARCHGJR, c("(Intercept) = 0", "PV_GARCHGJR = 1"), 
                 vcov = NeweyWest(mz_R_GARCHGJR, lag = 1))
summary(mz_R_GARCHGJR)


##RT-GJR##

# rv
mz_RV_RTGARCHGJR <- lm(TV_rv ~ PV_RTGARCHGJR)


linearHypothesis(mz_RV_RTGARCHGJR, c("(Intercept) = 0", "PV_RTGARCHGJR = 1"), 
                 vcov = NeweyWest(mz_RV_RTGARCHGJR, lag = 1))
summary(mz_RV_RTGARCHGJR)

# r

mz_R_RTGARCHGJR <- lm(TV_r ~ PV_RTGARCHGJR)


linearHypothesis(mz_R_RTGARCHGJR, c("(Intercept) = 0", "PV_RTGARCHGJR = 1"), 
                 vcov = NeweyWest(mz_R_RTGARCHGJR, lag = 1))
summary(mz_R_RTGARCHGJR)


#VIX

# rv

mz_RV_VIX <- lm(TV_rv ~ PV_VIX)


linearHypothesis(mz_RV_VIX, c("(Intercept) = 0", "PV_VIX = 1"), 
                 vcov = NeweyWest(mz_RV_VIX, lag = 1))
summary(mz_RV_VIX)

# r

mz_R_VIX <- lm(TV_r ~ PV_VIX)


linearHypothesis(mz_R_VIX, c("(Intercept) = 0", "PV_VIX = 1"), 
                 vcov = NeweyWest(mz_R_VIX, lag = 1))
summary(mz_R_VIX)

#HAR

# rv

mz_RV_HAR_RV <- lm(TV_rv ~ PV_HAR_RV)


linearHypothesis(mz_RV_HAR_RV, c("(Intercept) = 0", "PV_HAR_RV = 1"), 
                 vcov = NeweyWest(mz_RV_HAR_RV, lag = 1))
summary(mz_RV_HAR_RV)

# r

mz_R_HAR_RV <- lm(TV_r ~ PV_HAR_RV)


linearHypothesis(mz_R_HAR_RV, c("(Intercept) = 0", "PV_HAR_RV = 1"), 
                 vcov = NeweyWest(mz_R_HAR_RV, lag = 1))
summary(mz_R_HAR_RV)



