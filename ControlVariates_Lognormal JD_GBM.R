#CV
iter <- 20000
options(scipen = 999)
set.seed(1)

#Main Data
header <- c("Time", "Open", "High", "Low", "Last", "Change", "%Chg", "Volume", "Open Int")
data <- read.csv("harga.csv", header = TRUE, sep = ",")
data <- na.omit(data)
library(lubridate)
data$Time <- mdy(data$Time)
data$Open <- as.numeric(data$Open)
data$High <- as.numeric(data$High)
data$Low <- as.numeric(data$Low)
data$Last <- as.numeric(data$Last)
data <- data[order(data$Time, decreasing = FALSE), ]
rownames(data) <- 1:nrow(data)

days <- 365
K <- c(50, 60, 70, 80, 90, 100)
strike <- as.data.frame(K)

#Library
library(lubridate)
library(ggplot2)
library(Metrics)
library(openxlsx)
library(palmerpenguins)
library("psych") 

#Parameter tambahan untuk BS
w = 1
delta = 0
j = 0

#Fungsi untuk menghitung solusi analitik opsi asia geometrik
BS_asiageo <- function(w, So, K, r, delta, sigma, n, j, T){
  B = 1 
  h = 1 / 365
  Tmunj = (n-j)/n * (T - (h*(n-j-1))/2)
  Tnj = T * ((n-j)/n)^2 - ((n-j)*(n-j-1)*(4*n-4*j+1))/(6 * n^2) * h
  A = exp(-r*(T-Tmunj)-(sigma^2*(Tmunj - Tnj))/2)*B
  d = (log(So/K)+(r-delta-0.5*sigma^2)*Tmunj + log(B))/(sigma*sqrt(Tnj))
  v = w * So * A * exp(-delta*Tmunj) * pnorm(w*d + w*sigma*sqrt(Tnj)) - w*K*exp(-r*T)*pnorm(w*d)
  return(v)
}
#Bangkitkan Price Paths
jump_paths <- function(So, T, r, sigma, lam, mu_y, sigma_y, days, n) {
  S = NULL
  S[1] = So
  dt <- 1 / days
  
  for(i in 2:n){
    z <- rnorm(1,0,1)
    n_jump <- rpois(1, lam)
    sigma_n <- sqrt(sigma^2 + n_jump*sigma_y^2/dt)
    k <- exp(mu_y+1/2*sigma_y^2)-1
    S[i] <- S[i-1] * exp(-lam*k*dt+n_jump*mu_y+n_jump*sigma_y^2/2+(r-1/2*sigma_n^2)*dt+sigma_n*sqrt(dt)*z)
  }
  S <- S[-1]
  return(S)
}

gbm_paths <- function(So, T, r, sigma, days, n) {
  S1 = NULL
  S1[1] = So
  dt <- 1 / days
  
  for(i in 2:n){
    z1 <- rnorm(1,0,1)
    S1[i] <- S1[i-1] * exp((r-1/2*sigma^2)*dt+sigma*sqrt(dt)*z1)
  }
  S1 <- S1[-1]
  return(S1)
}

#Fungsi untuk menghitung harga opsi asia
cv_price_jump_call <- function(w, So, K, r, delta, sigma, n, j, T, lam, mu_y, sigma_y, days, iter){
  #Solusi Analitik BS
  CG <- BS_asiageo(w, So, K, r, delta, sigma, n, j, T)
  
  #Bentuk n price paths
  n_path = iter
  paths <- list()
  for(i in 1:n_path){
    S <- jump_paths(So, T, r, sigma, lam, mu_y, sigma_y, days, n)
    paths[[i]] <- S
  }
  
  S_data <- data.frame(Days = rep(1:(n-1), times = n_path), Stock_Price = unlist(paths),
                       Path = rep(1:n_path, each = (n-1)))
  
  #Hitung harga opsi Asia aritmatika
  C_AP = NULL
  h_A = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_A[i] = mean(dat$Stock_Price)
    C_AP[i] = exp(-r*T) * max(0, h_A[i] - K)
  }
  
  #Hitung harga opsi Asia geometrik
  C_GP = NULL
  h_G = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_G[i] = geometric.mean(dat$Stock_Price)
    C_GP[i] = exp(-r*T) * max(0, h_G[i] - K)
  }
  
  #hitung theta optimal
  theta_c <- NULL
  if(cov(C_AP, C_GP) == 0){
    theta_c = 0
  }else{
    theta_c = cov(C_AP, C_GP)/var(C_GP)
  }
  
  C_CV = NULL
  
  iter <- n_path
  for(i in 1:iter){
    C_CV <- c(C_CV, C_AP[i] + theta_c*(CG - C_GP[i])) 
  }
  
  #std err
  std_err <- sd(C_CV)/sqrt(length(C_CV))
  
  if(mean(C_CV) < 0){
    C_CV = 0
  } else{
    C_CV = mean(C_CV)
  }
  
  corr <- cor(C_AP, C_GP)
  
  Call_Price <- as.data.frame(list(Call_Jump = C_CV, cor = corr, std_error = std_err))
  return(Call_Price)
}
cv_price_jump_put <- function(w, So, K, r, delta, sigma, n, j, T, lam, mu_y, sigma_y, days, iter){
  #Solusi Analitik BS
  PG <- BS_asiageo(w=-1, So, K, r, delta, sigma, n, j, T)
  
  #Bentuk n price paths
  n_path = iter
  paths <- list()
  for(i in 1:n_path){
    S <- jump_paths(So, T, r, sigma, lam, mu_y, sigma_y, days, n)
    paths[[i]] <- S
  }
  
  S_data <- data.frame(Days = rep(1:(n-1), times = n_path), Stock_Price = unlist(paths),
                       Path = rep(1:n_path, each = (n-1)))
  
  #Hitung harga opsi Asia aritmatika
  P_AP = NULL
  h_A = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_A[i] = mean(dat$Stock_Price)
    P_AP[i] = exp(-r*T) * max(0, K - h_A[i])
  }
  
  #Hitung harga opsi Asia geometrik
  P_GP = NULL
  h_G = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_G[i] = geometric.mean(dat$Stock_Price)
    P_GP[i] = exp(-r*T) * max(0, K - h_G[i])
  }
  
  #hitung theta optimal
  theta_p <- NULL
  if(cov(P_AP, P_GP) == 0){
    theta_p = 0
  }else{
    theta_p = cov(P_AP, P_GP)/var(P_GP)
  }
  
  P_CV = NULL
  
  iter <- n_path
  for(i in 1:iter){
    P_CV <- c(P_CV, P_AP[i] + theta_p*(PG - P_GP[i]))
  }
  
  #std err
  std_err <- sd(P_CV)/sqrt(length(P_CV))
  
  if(mean(P_CV) < 0){
    P_CV = 0
  } else{
    P_CV = mean(P_CV)
  }
  
  corr <- cor(P_AP, P_GP)
  
  Put_Price <- as.data.frame(list(Put_Jump =  P_CV, cor = corr, std_error = std_err))
  return(Put_Price)
}
cv_price_gbm_call <- function(w, So, K, r, delta, sigma, n, j, T, days, iter){
  #Solusi Analitik BS
  CG <- BS_asiageo(w, So, K, r, delta, sigma, n, j, T)
  
  #Bentuk n price paths
  n_path = iter
  paths <- list()
  for(i in 1:n_path){
    S <- gbm_paths(So, T, r, sigma, days, n)
    paths[[i]] <- S
  }
  
  S_data <- data.frame(Days = rep(1:(n-1), times = n_path), Stock_Price = unlist(paths),
                       Path = rep(1:n_path, each = (n-1)))
  
  #Hitung harga opsi Asia aritmatika
  C_AP = NULL
  h_A = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_A[i] = mean(dat$Stock_Price)
    C_AP[i] = exp(-r*T) * max(0, h_A[i] - K)
  }
  
  #Hitung harga opsi Asia geometrik
  C_GP = NULL
  h_G = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_G[i] = geometric.mean(dat$Stock_Price)
    C_GP[i] = exp(-r*T) * max(0, h_G[i] - K)
  }
  
  #hitung theta optimal
  theta_c <- NULL
  if(cov(C_AP, C_GP) == 0){
    theta_c = 0
  }else{
    theta_c = cov(C_AP, C_GP)/var(C_GP)
  }
  
  C_CV = NULL
  
  iter <- n_path
  for(i in 1:iter){
    C_CV <- c(C_CV, C_AP[i] + theta_c*(CG - C_GP[i]))
  }
  
  #std err
  std_err <- sd(C_CV)/sqrt(length(C_CV))
  
  C_CV = mean(C_CV)
  
  corr <- cor(C_AP, C_GP)
  
  Call_Price <- as.data.frame(list(Call_GBM = C_CV, cor = corr, std_error = std_err))
  return(Call_Price)
}

cv_price_gbm_put <- function(w, So, K, r, delta, sigma, n, j, T, days, iter){
  #Solusi Analitik BS
  PG <- BS_asiageo(w=-1, So, K, r, delta, sigma, n, j, T)
  
  #Bentuk n price paths
  n_path = iter
  paths <- list()
  for(i in 1:n_path){
    S <- gbm_paths(So, T, r, sigma, days, n)
    paths[[i]] <- S
  }
  
  S_data <- data.frame(Days = rep(1:(n-1), times = n_path), Stock_Price = unlist(paths),
                       Path = rep(1:n_path, each = (n-1)))
  
  #Hitung harga opsi Asia aritmatika
  P_AP = NULL
  h_A = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_A[i] = mean(dat$Stock_Price)
    P_AP[i] = exp(-r*T) * max(0, K - h_A[i])
  }
  
  #Hitung harga opsi Asia geometrik
  P_GP = NULL
  h_G = NULL
  
  for (i in 1:n_path) {
    dat <- subset(S_data, Path == i)
    h_G[i] = geometric.mean(dat$Stock_Price)
    P_GP[i] = exp(-r*T) * max(0, K - h_G[i])
  }
  
  #hitung theta optimal
  theta_p <- NULL
  if(cov(P_AP, P_GP) == 0){
    theta_p = 0
  }else{
    theta_p = cov(P_AP, P_GP)/var(P_GP)
  }
  
  P_CV = NULL
  
  iter <- n_path
  for(i in 1:iter){
    P_CV <- c(P_CV, P_AP[i] + theta_p*(PG - P_GP[i]))
  }
  
  #std err
  std_err <- sd(P_CV)/sqrt(length(P_CV))
  
  P_CV = mean(P_CV)
  
  corr <- cor(P_AP, P_GP)
  
  Put_Price <- as.data.frame(list(Put_GBM =  P_CV, cor = corr, std_error = std_err))
  return(Put_Price)
}

#16 Days
#data
data_asli_16days <- data[c(1462:1473),]
data_16days_1YR <- data[c(1210:1461),]
data_16days_5YR <- data[c(203:1461),]
At_16days <- mean(data_asli_16days$Last)

time_to_maturity_16days <- 16
T_16days <- time_to_maturity_16days/days  
i_16days <- 0.05396 #US Treasury Yield (1 month) 24 OCT
r_16days <- log(1+i_16days)
n_16days <- nrow(data_asli_16days) + 1

call_realprice_16days <- c(37.74, 27.76, 17.79, 7.81, 0.12, 0.01)
put_realprice_16days <- c(0.01, 0.01, 0.01, 0.01, 2.28, 12.14)
real_price_16days <- as.data.frame(list(Call_RealPrice = call_realprice_16days, Put_RealPrice = put_realprice_16days))

real_payoff_call_16days <- NULL
real_payoff_put_16days <- NULL
for (i in 1:length(K)) {
  real_payoff_call_16days <- c(real_payoff_call_16days, max(At_16days - K[i],0)*exp(-r_16days*T_16days))
  real_payoff_put_16days <- c(real_payoff_put_16days, max(K[i] - At_16days,0)*exp(-r_16days*T_16days))
} 

#1YR
#Menghitung Logreturn
log_return_16days_1YR <- NULL
log_return_16days_1YR[1] <- 0
for(i in 2:nrow(data_16days_1YR)){
  log_return_16days_1YR[i] = log(data_16days_1YR$Last[i]/data_16days_1YR$Last[i-1])
}
data_16days_1YR$`log return` <- log_return_16days_1YR

#Volatilitas
vol_16days_1YR <- sd(data_16days_1YR$`log return`)

#Jump Parameter
is.jump_16days_1YR <- which(abs(data_16days_1YR$`log return`) > 2*sd(data_16days_1YR$`log return`))
num_jump_16days_1YR <- length(is.jump_16days_1YR)

mu_y_16days_1YR <- mean(data_16days_1YR$`log return`[is.jump_16days_1YR])
sigma_y_16days_1YR <- sd(data_16days_1YR$`log return`[is.jump_16days_1YR])
lam_16days_1YR <- num_jump_16days_1YR/nrow(data_16days_1YR)
So_16days_1YR <- data_16days_1YR$Last[nrow(data_16days_1YR)] 
sigma_16days_1YR <- sd(data_16days_1YR$`log return`)*sqrt(time_to_maturity_16days)

p_16days_1YR <- ggplot(data_16days_1YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_16days_1YR[is.jump_16days_1YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 13 September 2022 - 13 September 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_16days_1YR)

cv_price_option_jump_16days_1YR_call <- NULL
cv_price_option_jump_16days_1YR_put <- NULL
cv_price_option_gbm_16days_1YR_call <- NULL
cv_price_option_gbm_16days_1YR_put <- NULL

start_16days_1YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_16days_1YR_call <- rbind(cv_price_option_jump_16days_1YR_call, cv_price_jump_call(w, So = So_16days_1YR, K[i], r = r_16days, delta, sigma = sigma_16days_1YR, n = n_16days, j, 
                                                                                                         T = T_16days, lam = lam_16days_1YR, mu_y=mu_y_16days_1YR, sigma_y=sigma_y_16days_1YR, days, iter))
} 
end_16days_1YR_JUMP_call <- Sys.time()

start_16days_1YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_16days_1YR_put <- rbind(cv_price_option_jump_16days_1YR_put, cv_price_jump_put(w, So = So_16days_1YR, K[i], r = r_16days, delta, sigma = sigma_16days_1YR, n = n_16days, j, 
                                                                                                      T = T_16days, lam = lam_16days_1YR, mu_y=mu_y_16days_1YR, sigma_y=sigma_y_16days_1YR, days, iter))
}
end_16days_1YR_JUMP_put <- Sys.time()

start_16days_1YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_16days_1YR_call <- rbind(cv_price_option_gbm_16days_1YR_call, cv_price_gbm_call(w, So = So_16days_1YR, K[i], r = r_16days, delta, sigma = sigma_16days_1YR, n = n_16days, j, T = T_16days, days, iter))
  
}
end_16days_1YR_GBM_call <- Sys.time()

start_16days_1YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_16days_1YR_put <- rbind(cv_price_option_gbm_16days_1YR_put, cv_price_gbm_put(w, So = So_16days_1YR, K[i], r = r_16days, delta, sigma = sigma_16days_1YR, n = n_16days, j, T = T_16days, days, iter))
  
}
end_16days_1YR_GBM_put <- Sys.time()

time_16days_1YR_JUMP_call <- start_16days_1YR_JUMP_call - end_16days_1YR_JUMP_call
time_16days_1YR_JUMP_put <- start_16days_1YR_JUMP_put - end_16days_1YR_JUMP_put
time_16days_1YR_GBM_call <- start_16days_1YR_GBM_call - end_16days_1YR_GBM_call
time_16days_1YR_GBM_put <- start_16days_1YR_GBM_put - end_16days_1YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_16days_1YR_call <- cbind(cv_price_option_jump_16days_1YR_call$Call_Jump, cv_price_option_gbm_16days_1YR_call$Call_GBM)
cv_price_option_16days_1YR_call <- cbind(strike, cv_price_option_16days_1YR_call, real_price_16days$Call_RealPrice, real_payoff_call_16days)
colnames(cv_price_option_16days_1YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_16days_1YR_put <- cbind(cv_price_option_jump_16days_1YR_put$Put_Jump, cv_price_option_gbm_16days_1YR_put$Put_GBM)
cv_price_option_16days_1YR_put <- cbind(strike, cv_price_option_16days_1YR_put, real_price_16days$Put_RealPrice, real_payoff_put_16days)
colnames(cv_price_option_16days_1YR_put) <- col_names_put

rmse_call_jump_16days_1YR <- rmse(cv_price_option_16days_1YR_call$Barchart, cv_price_option_16days_1YR_call$`Call Jump`)
rmse_call_gbm_16days_1YR <- rmse(cv_price_option_16days_1YR_call$Barchart, cv_price_option_16days_1YR_call$`Call GBM`)
rmse_put_jump_16days_1YR <- rmse(cv_price_option_16days_1YR_put$Barchart, cv_price_option_16days_1YR_put$`Put Jump`)
rmse_put_gbm_16days_1YR <- rmse(cv_price_option_16days_1YR_put$Barchart, cv_price_option_16days_1YR_put$`Put GBM`)

rmse_call_jump_real_16days_1YR <- rmse(cv_price_option_16days_1YR_call$`Payoff Realisasi`, cv_price_option_16days_1YR_call$`Call Jump`)
rmse_call_gbm_real_16days_1YR <- rmse(cv_price_option_16days_1YR_call$`Payoff Realisasi`, cv_price_option_16days_1YR_call$`Call GBM`)
rmse_put_jump_real_16days_1YR <- rmse(cv_price_option_16days_1YR_put$`Payoff Realisasi`, cv_price_option_16days_1YR_put$`Put Jump`)
rmse_put_gbm_real_16days_1YR <- rmse(cv_price_option_16days_1YR_put$`Payoff Realisasi`, cv_price_option_16days_1YR_put$`Put GBM`)

#5YR
#Menghitung Logreturn
log_return_16days_5YR <- NULL
log_return_16days_5YR[1] <- 0
for(i in 2:nrow(data_16days_5YR)){
  log_return_16days_5YR[i] = log(data_16days_5YR$Last[i]/data_16days_5YR$Last[i-1])
}
data_16days_5YR$`log return` <- log_return_16days_5YR

#Volatilitas
vol_16days_5YR <- sd(data_16days_5YR$`log return`)

#Jump Parameter
is.jump_16days_5YR <- which(abs(data_16days_5YR$`log return`) > 2*sd(data_16days_5YR$`log return`))
num_jump_16days_5YR <- length(is.jump_16days_5YR)

mu_y_16days_5YR <- mean(data_16days_5YR$`log return`[is.jump_16days_5YR])
sigma_y_16days_5YR <- sd(data_16days_5YR$`log return`[is.jump_16days_5YR])
lam_16days_5YR <- num_jump_16days_5YR/nrow(data_16days_5YR)
So_16days_5YR <- data_16days_5YR$Last[nrow(data_16days_5YR)] 
sigma_16days_5YR <- sd(data_16days_5YR$`log return`)*sqrt(time_to_maturity_16days)

p_16days_5YR <- ggplot(data_16days_5YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_16days_5YR[is.jump_16days_5YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 13 September 2018 - 13 September 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_16days_5YR)

cv_price_option_jump_16days_5YR_call <- NULL
cv_price_option_jump_16days_5YR_put <- NULL
cv_price_option_gbm_16days_5YR_call <- NULL
cv_price_option_gbm_16days_5YR_put <- NULL

start_16days_5YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_16days_5YR_call <- rbind(cv_price_option_jump_16days_5YR_call, cv_price_jump_call(w, So = So_16days_5YR, K[i], r = r_16days, delta, sigma = sigma_16days_5YR, n = n_16days, j, 
                                                                                                         T = T_16days, lam = lam_16days_5YR, mu_y=mu_y_16days_5YR, sigma_y=sigma_y_16days_5YR, days, iter))
} 
end_16days_5YR_JUMP_call <- Sys.time()

start_16days_5YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_16days_5YR_put <- rbind(cv_price_option_jump_16days_5YR_put, cv_price_jump_put(w, So = So_16days_5YR, K[i], r = r_16days, delta, sigma = sigma_16days_5YR, n = n_16days, j, 
                                                                                                      T = T_16days, lam = lam_16days_5YR, mu_y=mu_y_16days_5YR, sigma_y=sigma_y_16days_5YR, days, iter))
}
end_16days_5YR_JUMP_put <- Sys.time()

start_16days_5YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_16days_5YR_call <- rbind(cv_price_option_gbm_16days_5YR_call, cv_price_gbm_call(w, So = So_16days_5YR, K[i], r = r_16days, delta, sigma = sigma_16days_5YR, n = n_16days, j, T = T_16days, days, iter))
  
}
end_16days_5YR_GBM_call <- Sys.time()

start_16days_5YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_16days_5YR_put <- rbind(cv_price_option_gbm_16days_5YR_put, cv_price_gbm_put(w, So = So_16days_5YR, K[i], r = r_16days, delta, sigma = sigma_16days_5YR, n = n_16days, j, T = T_16days, days, iter))
  
}
end_16days_5YR_GBM_put <- Sys.time()

time_16days_5YR_JUMP_call <- start_16days_5YR_JUMP_call - end_16days_5YR_JUMP_call
time_16days_5YR_JUMP_put <- start_16days_5YR_JUMP_put - end_16days_5YR_JUMP_put
time_16days_5YR_GBM_call <- start_16days_5YR_GBM_call - end_16days_5YR_GBM_call
time_16days_5YR_GBM_put <- start_16days_5YR_GBM_put - end_16days_5YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_16days_5YR_call <- cbind(cv_price_option_jump_16days_5YR_call$Call_Jump, cv_price_option_gbm_16days_5YR_call$Call_GBM)
cv_price_option_16days_5YR_call <- cbind(strike, cv_price_option_16days_5YR_call, real_price_16days$Call_RealPrice, real_payoff_call_16days)
colnames(cv_price_option_16days_5YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_16days_5YR_put <- cbind(cv_price_option_jump_16days_5YR_put$Put_Jump, cv_price_option_gbm_16days_5YR_put$Put_GBM)
cv_price_option_16days_5YR_put <- cbind(strike, cv_price_option_16days_5YR_put, real_price_16days$Put_RealPrice, real_payoff_put_16days)
colnames(cv_price_option_16days_5YR_put) <- col_names_put

rmse_call_jump_16days_5YR <- rmse(cv_price_option_16days_5YR_call$Barchart, cv_price_option_16days_5YR_call$`Call Jump`)
rmse_call_gbm_16days_5YR <- rmse(cv_price_option_16days_5YR_call$Barchart, cv_price_option_16days_5YR_call$`Call GBM`)
rmse_put_jump_16days_5YR <- rmse(cv_price_option_16days_5YR_put$Barchart, cv_price_option_16days_5YR_put$`Put Jump`)
rmse_put_gbm_16days_5YR <- rmse(cv_price_option_16days_5YR_put$Barchart, cv_price_option_16days_5YR_put$`Put GBM`)

rmse_call_jump_real_16days_5YR <- rmse(cv_price_option_16days_5YR_call$`Payoff Realisasi`, cv_price_option_16days_5YR_call$`Call Jump`)
rmse_call_gbm_real_16days_5YR <- rmse(cv_price_option_16days_5YR_call$`Payoff Realisasi`, cv_price_option_16days_5YR_call$`Call GBM`)
rmse_put_jump_real_16days_5YR <- rmse(cv_price_option_16days_5YR_put$`Payoff Realisasi`, cv_price_option_16days_5YR_put$`Put Jump`)
rmse_put_gbm_real_16days_5YR <- rmse(cv_price_option_16days_5YR_put$`Payoff Realisasi`, cv_price_option_16days_5YR_put$`Put GBM`)

#1 Month
#data
data_asli_1month <- data[c(1452:1473),]
data_1month_1YR <- data[c(1200:1451),]
data_1month_5YR <- data[c(193:1451),]
At_1month <- mean(data_asli_1month$Last)

time_to_maturity_1month <- 31
T_1month <- time_to_maturity_1month/days  
i_1month <- 0.05396 #US Treasury Yield (1 month) 24 OCT
r_1month <- log(1+i_1month)
n_1month <- nrow(data_asli_1month) + 1

call_realprice_1month <- c(30.85, 20.9, 10.98, 2.25, 0.04, 0.01)
put_realprice_1month <- c(0.01, 0.01, 0.04, 1.26, 9, 18.92)
real_price_1month <- as.data.frame(list(Call_RealPrice = call_realprice_1month, Put_RealPrice = put_realprice_1month))

real_payoff_call_1month <- NULL
real_payoff_put_1month <- NULL
for (i in 1:length(K)) {
  real_payoff_call_1month <- c(real_payoff_call_1month, max(At_1month - K[i],0)*exp(-r_1month*T_1month))
  real_payoff_put_1month <- c(real_payoff_put_1month, max(K[i] - At_1month,0)*exp(-r_1month*T_1month))
} 

#1YR
#Menghitung Logreturn
log_return_1month_1YR <- NULL
log_return_1month_1YR[1] <- 0
for(i in 2:nrow(data_1month_1YR)){
  log_return_1month_1YR[i] = log(data_1month_1YR$Last[i]/data_1month_1YR$Last[i-1])
}
data_1month_1YR$`log return` <- log_return_1month_1YR

#Volatilitas
vol_1month_1YR <- sd(data_1month_1YR$`log return`)

#Jump Parameter
is.jump_1month_1YR <- which(abs(data_1month_1YR$`log return`) > 2*sd(data_1month_1YR$`log return`))
num_jump_1month_1YR <- length(is.jump_1month_1YR)

mu_y_1month_1YR <- mean(data_1month_1YR$`log return`[is.jump_1month_1YR])
sigma_y_1month_1YR <- sd(data_1month_1YR$`log return`[is.jump_1month_1YR])
lam_1month_1YR <- num_jump_1month_1YR/nrow(data_1month_1YR)
So_1month_1YR <- data_1month_1YR$Last[nrow(data_1month_1YR)] 
sigma_1month_1YR <- sd(data_1month_1YR$`log return`)*sqrt(time_to_maturity_1month)

p_1month_1YR <- ggplot(data_1month_1YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_1month_1YR[is.jump_1month_1YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 29 Agustus 2022 - 29 Agustus 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_1month_1YR)

cv_price_option_jump_1month_1YR_call <- NULL
cv_price_option_jump_1month_1YR_put <- NULL
cv_price_option_gbm_1month_1YR_call <- NULL
cv_price_option_gbm_1month_1YR_put <- NULL

start_1month_1YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_1month_1YR_call <- rbind(cv_price_option_jump_1month_1YR_call, cv_price_jump_call(w, So = So_1month_1YR, K[i], r = r_1month, delta, sigma = sigma_1month_1YR, n = n_1month, j, 
                                                                                                         T = T_1month, lam = lam_1month_1YR, mu_y=mu_y_1month_1YR, sigma_y=sigma_y_1month_1YR, days, iter))
} 
end_1month_1YR_JUMP_call <- Sys.time()

start_1month_1YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_1month_1YR_put <- rbind(cv_price_option_jump_1month_1YR_put, cv_price_jump_put(w, So = So_1month_1YR, K[i], r = r_1month, delta, sigma = sigma_1month_1YR, n = n_1month, j, 
                                                                                                      T = T_1month, lam = lam_1month_1YR, mu_y=mu_y_1month_1YR, sigma_y=sigma_y_1month_1YR, days, iter))
}
end_1month_1YR_JUMP_put <- Sys.time()

start_1month_1YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_1month_1YR_call <- rbind(cv_price_option_gbm_1month_1YR_call, cv_price_gbm_call(w, So = So_1month_1YR, K[i], r = r_1month, delta, sigma = sigma_1month_1YR, n = n_1month, j, T = T_1month, days, iter))
  
}
end_1month_1YR_GBM_call <- Sys.time()

start_1month_1YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_1month_1YR_put <- rbind(cv_price_option_gbm_1month_1YR_put, cv_price_gbm_put(w, So = So_1month_1YR, K[i], r = r_1month, delta, sigma = sigma_1month_1YR, n = n_1month, j, T = T_1month, days, iter))
  
}
end_1month_1YR_GBM_put <- Sys.time()

time_1month_1YR_JUMP_call <- start_1month_1YR_JUMP_call - end_1month_1YR_JUMP_call
time_1month_1YR_JUMP_put <- start_1month_1YR_JUMP_put - end_1month_1YR_JUMP_put
time_1month_1YR_GBM_call <- start_1month_1YR_GBM_call - end_1month_1YR_GBM_call
time_1month_1YR_GBM_put <- start_1month_1YR_GBM_put - end_1month_1YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_1month_1YR_call <- cbind(cv_price_option_jump_1month_1YR_call$Call_Jump, cv_price_option_gbm_1month_1YR_call$Call_GBM)
cv_price_option_1month_1YR_call <- cbind(strike, cv_price_option_1month_1YR_call, real_price_1month$Call_RealPrice, real_payoff_call_1month)
colnames(cv_price_option_1month_1YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_1month_1YR_put <- cbind(cv_price_option_jump_1month_1YR_put$Put_Jump, cv_price_option_gbm_1month_1YR_put$Put_GBM)
cv_price_option_1month_1YR_put <- cbind(strike, cv_price_option_1month_1YR_put, real_price_1month$Put_RealPrice, real_payoff_put_1month)
colnames(cv_price_option_1month_1YR_put) <- col_names_put

rmse_call_jump_1month_1YR <- rmse(cv_price_option_1month_1YR_call$Barchart, cv_price_option_1month_1YR_call$`Call Jump`)
rmse_call_gbm_1month_1YR <- rmse(cv_price_option_1month_1YR_call$Barchart, cv_price_option_1month_1YR_call$`Call GBM`)
rmse_put_jump_1month_1YR <- rmse(cv_price_option_1month_1YR_put$Barchart, cv_price_option_1month_1YR_put$`Put Jump`)
rmse_put_gbm_1month_1YR <- rmse(cv_price_option_1month_1YR_put$Barchart, cv_price_option_1month_1YR_put$`Put GBM`)

rmse_call_jump_real_1month_1YR <- rmse(cv_price_option_1month_1YR_call$`Payoff Realisasi`, cv_price_option_1month_1YR_call$`Call Jump`)
rmse_call_gbm_real_1month_1YR <- rmse(cv_price_option_1month_1YR_call$`Payoff Realisasi`, cv_price_option_1month_1YR_call$`Call GBM`)
rmse_put_jump_real_1month_1YR <- rmse(cv_price_option_1month_1YR_put$`Payoff Realisasi`, cv_price_option_1month_1YR_put$`Put Jump`)
rmse_put_gbm_real_1month_1YR <- rmse(cv_price_option_1month_1YR_put$`Payoff Realisasi`, cv_price_option_1month_1YR_put$`Put GBM`)

#5YR
#Menghitung Logreturn
log_return_1month_5YR <- NULL
log_return_1month_5YR[1] <- 0
for(i in 2:nrow(data_1month_5YR)){
  log_return_1month_5YR[i] = log(data_1month_5YR$Last[i]/data_1month_5YR$Last[i-1])
}
data_1month_5YR$`log return` <- log_return_1month_5YR

#Volatilitas
vol_1month_5YR <- sd(data_1month_5YR$`log return`)

#Jump Parameter
is.jump_1month_5YR <- which(abs(data_1month_5YR$`log return`) > 2*sd(data_1month_5YR$`log return`))
num_jump_1month_5YR <- length(is.jump_1month_5YR)

mu_y_1month_5YR <- mean(data_1month_5YR$`log return`[is.jump_1month_5YR])
sigma_y_1month_5YR <- sd(data_1month_5YR$`log return`[is.jump_1month_5YR])
lam_1month_5YR <- num_jump_1month_5YR/nrow(data_1month_5YR)
So_1month_5YR <- data_1month_5YR$Last[nrow(data_1month_5YR)] 
sigma_1month_5YR <- sd(data_1month_5YR$`log return`)*sqrt(time_to_maturity_1month)

p_1month_5YR <- ggplot(data_1month_5YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_1month_5YR[is.jump_1month_5YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 29 Agustus 2018 - 29 Agustus 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_1month_5YR)

cv_price_option_jump_1month_5YR_call <- NULL
cv_price_option_jump_1month_5YR_put <- NULL
cv_price_option_gbm_1month_5YR_call <- NULL
cv_price_option_gbm_1month_5YR_put <- NULL

start_1month_5YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_1month_5YR_call <- rbind(cv_price_option_jump_1month_5YR_call, cv_price_jump_call(w, So = So_1month_5YR, K[i], r = r_1month, delta, sigma = sigma_1month_5YR, n = n_1month, j, 
                                                                                                         T = T_1month, lam = lam_1month_5YR, mu_y=mu_y_1month_5YR, sigma_y=sigma_y_1month_5YR, days, iter))
} 
end_1month_5YR_JUMP_call <- Sys.time()

start_1month_5YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_1month_5YR_put <- rbind(cv_price_option_jump_1month_5YR_put, cv_price_jump_put(w, So = So_1month_5YR, K[i], r = r_1month, delta, sigma = sigma_1month_5YR, n = n_1month, j, 
                                                                                                      T = T_1month, lam = lam_1month_5YR, mu_y=mu_y_1month_5YR, sigma_y=sigma_y_1month_5YR, days, iter))
}
end_1month_5YR_JUMP_put <- Sys.time()

start_1month_5YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_1month_5YR_call <- rbind(cv_price_option_gbm_1month_5YR_call, cv_price_gbm_call(w, So = So_1month_5YR, K[i], r = r_1month, delta, sigma = sigma_1month_5YR, n = n_1month, j, T = T_1month, days, iter))
  
}
end_1month_5YR_GBM_call <- Sys.time()

start_1month_5YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_1month_5YR_put <- rbind(cv_price_option_gbm_1month_5YR_put, cv_price_gbm_put(w, So = So_1month_5YR, K[i], r = r_1month, delta, sigma = sigma_1month_5YR, n = n_1month, j, T = T_1month, days, iter))
  
}
end_1month_5YR_GBM_put <- Sys.time()

time_1month_5YR_JUMP_call <- start_1month_5YR_JUMP_call - end_1month_5YR_JUMP_call
time_1month_5YR_JUMP_put <- start_1month_5YR_JUMP_put - end_1month_5YR_JUMP_put
time_1month_5YR_GBM_call <- start_1month_5YR_GBM_call - end_1month_5YR_GBM_call
time_1month_5YR_GBM_put <- start_1month_5YR_GBM_put - end_1month_5YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_1month_5YR_call <- cbind(cv_price_option_jump_1month_5YR_call$Call_Jump, cv_price_option_gbm_1month_5YR_call$Call_GBM)
cv_price_option_1month_5YR_call <- cbind(strike, cv_price_option_1month_5YR_call, real_price_1month$Call_RealPrice, real_payoff_call_1month)
colnames(cv_price_option_1month_5YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_1month_5YR_put <- cbind(cv_price_option_jump_1month_5YR_put$Put_Jump, cv_price_option_gbm_1month_5YR_put$Put_GBM)
cv_price_option_1month_5YR_put <- cbind(strike, cv_price_option_1month_5YR_put, real_price_1month$Put_RealPrice, real_payoff_put_1month)
colnames(cv_price_option_1month_5YR_put) <- col_names_put

rmse_call_jump_1month_5YR <- rmse(cv_price_option_1month_5YR_call$Barchart, cv_price_option_1month_5YR_call$`Call Jump`)
rmse_call_gbm_1month_5YR <- rmse(cv_price_option_1month_5YR_call$Barchart, cv_price_option_1month_5YR_call$`Call GBM`)
rmse_put_jump_1month_5YR <- rmse(cv_price_option_1month_5YR_put$Barchart, cv_price_option_1month_5YR_put$`Put Jump`)
rmse_put_gbm_1month_5YR <- rmse(cv_price_option_1month_5YR_put$Barchart, cv_price_option_1month_5YR_put$`Put GBM`)

rmse_call_jump_real_1month_5YR <- rmse(cv_price_option_1month_5YR_call$`Payoff Realisasi`, cv_price_option_1month_5YR_call$`Call Jump`)
rmse_call_gbm_real_1month_5YR <- rmse(cv_price_option_1month_5YR_call$`Payoff Realisasi`, cv_price_option_1month_5YR_call$`Call GBM`)
rmse_put_jump_real_1month_5YR <- rmse(cv_price_option_1month_5YR_put$`Payoff Realisasi`, cv_price_option_1month_5YR_put$`Put Jump`)
rmse_put_gbm_real_1month_5YR <- rmse(cv_price_option_1month_5YR_put$`Payoff Realisasi`, cv_price_option_1month_5YR_put$`Put GBM`)

#3 months
#data
data_asli_3months <- data[c(1410:1473),]
data_3months_1YR <- data[c(1158:1409),]
data_3months_5YR <- data[c(151:1409),]
At_3months <- mean(data_asli_3months$Last)

time_to_maturity_3months <- 92
T_3months <- time_to_maturity_3months/days  
i_3months <- 0.05456 #US Treasury Yield (3 months)
r_3months <- log(1+i_3months)
n_3months <- nrow(data_asli_3months) + 1

call_realprice_3months <- c(19.98, 11.01, 4.33, 1.06, 0.23, 0.08)
put_realprice_3months <- c(0.23, 1.12, 4.31, 10.91, 19.95, 29.66)
real_price_3months <- as.data.frame(list(Call_RealPrice = call_realprice_3months, Put_RealPrice = put_realprice_3months))

real_payoff_call_3months <- NULL
real_payoff_put_3months <- NULL
for (i in 1:length(K)) {
  real_payoff_call_3months <- c(real_payoff_call_3months, max(At_3months - K[i],0)*exp(-r_3months*T_3months))
  real_payoff_put_3months <- c(real_payoff_put_3months, max(K[i] - At_3months,0)*exp(-r_3months*T_3months))
} 

#1YR
#Menghitung Logreturn
log_return_3months_1YR <- NULL
log_return_3months_1YR[1] <- 0
for(i in 2:nrow(data_3months_1YR)){
  log_return_3months_1YR[i] = log(data_3months_1YR$Last[i]/data_3months_1YR$Last[i-1])
}
data_3months_1YR$`log return` <- log_return_3months_1YR

#Volatilitas
vol_3months_1YR <- sd(data_3months_1YR$`log return`)

#Jump Parameter
is.jump_3months_1YR <- which(abs(data_3months_1YR$`log return`) > 2*sd(data_3months_1YR$`log return`))
num_jump_3months_1YR <- length(is.jump_3months_1YR)

mu_y_3months_1YR <- mean(data_3months_1YR$`log return`[is.jump_3months_1YR])
sigma_y_3months_1YR <- sd(data_3months_1YR$`log return`[is.jump_3months_1YR])
lam_3months_1YR <- num_jump_3months_1YR/nrow(data_3months_1YR)
So_3months_1YR <- data_3months_1YR$Last[nrow(data_3months_1YR)] 
sigma_3months_1YR <- sd(data_3months_1YR$`log return`)*sqrt(time_to_maturity_3months)

p_3months_1YR <- ggplot(data_3months_1YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_3months_1YR[is.jump_3months_1YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 29 Juni 2022 - 29 Juni 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_3months_1YR)

cv_price_option_jump_3months_1YR_call <- NULL
cv_price_option_jump_3months_1YR_put <- NULL
cv_price_option_gbm_3months_1YR_call <- NULL
cv_price_option_gbm_3months_1YR_put <- NULL

start_3months_1YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_3months_1YR_call <- rbind(cv_price_option_jump_3months_1YR_call, cv_price_jump_call(w, So = So_3months_1YR, K[i], r = r_3months, delta, sigma = sigma_3months_1YR, n = n_3months, j, 
                                                                                                           T = T_3months, lam = lam_3months_1YR, mu_y=mu_y_3months_1YR, sigma_y=sigma_y_3months_1YR, days, iter))
} 
end_3months_1YR_JUMP_call <- Sys.time()

start_3months_1YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_3months_1YR_put <- rbind(cv_price_option_jump_3months_1YR_put, cv_price_jump_put(w, So = So_3months_1YR, K[i], r = r_3months, delta, sigma = sigma_3months_1YR, n = n_3months, j, 
                                                                                                        T = T_3months, lam = lam_3months_1YR, mu_y=mu_y_3months_1YR, sigma_y=sigma_y_3months_1YR, days, iter))
}
end_3months_1YR_JUMP_put <- Sys.time()

start_3months_1YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_3months_1YR_call <- rbind(cv_price_option_gbm_3months_1YR_call, cv_price_gbm_call(w, So = So_3months_1YR, K[i], r = r_3months, delta, sigma = sigma_3months_1YR, n = n_3months, j, T = T_3months, days, iter))
  
}
end_3months_1YR_GBM_call <- Sys.time()

start_3months_1YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_3months_1YR_put <- rbind(cv_price_option_gbm_3months_1YR_put, cv_price_gbm_put(w, So = So_3months_1YR, K[i], r = r_3months, delta, sigma = sigma_3months_1YR, n = n_3months, j, T = T_3months, days, iter))
  
}
end_3months_1YR_GBM_put <- Sys.time()

time_3months_1YR_JUMP_call <- start_3months_1YR_JUMP_call - end_3months_1YR_JUMP_call
time_3months_1YR_JUMP_put <- start_3months_1YR_JUMP_put - end_3months_1YR_JUMP_put
time_3months_1YR_GBM_call <- start_3months_1YR_GBM_call - end_3months_1YR_GBM_call
time_3months_1YR_GBM_put <- start_3months_1YR_GBM_put - end_3months_1YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_3months_1YR_call <- cbind(cv_price_option_jump_3months_1YR_call$Call_Jump, cv_price_option_gbm_3months_1YR_call$Call_GBM)
cv_price_option_3months_1YR_call <- cbind(strike, cv_price_option_3months_1YR_call, real_price_3months$Call_RealPrice, real_payoff_call_3months)
colnames(cv_price_option_3months_1YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_3months_1YR_put <- cbind(cv_price_option_jump_3months_1YR_put$Put_Jump, cv_price_option_gbm_3months_1YR_put$Put_GBM)
cv_price_option_3months_1YR_put <- cbind(strike, cv_price_option_3months_1YR_put, real_price_3months$Put_RealPrice, real_payoff_put_3months)
colnames(cv_price_option_3months_1YR_put) <- col_names_put

rmse_call_jump_3months_1YR <- rmse(cv_price_option_3months_1YR_call$Barchart, cv_price_option_3months_1YR_call$`Call Jump`)
rmse_call_gbm_3months_1YR <- rmse(cv_price_option_3months_1YR_call$Barchart, cv_price_option_3months_1YR_call$`Call GBM`)
rmse_put_jump_3months_1YR <- rmse(cv_price_option_3months_1YR_put$Barchart, cv_price_option_3months_1YR_put$`Put Jump`)
rmse_put_gbm_3months_1YR <- rmse(cv_price_option_3months_1YR_put$Barchart, cv_price_option_3months_1YR_put$`Put GBM`)

rmse_call_jump_real_3months_1YR <- rmse(cv_price_option_3months_1YR_call$`Payoff Realisasi`, cv_price_option_3months_1YR_call$`Call Jump`)
rmse_call_gbm_real_3months_1YR <- rmse(cv_price_option_3months_1YR_call$`Payoff Realisasi`, cv_price_option_3months_1YR_call$`Call GBM`)
rmse_put_jump_real_3months_1YR <- rmse(cv_price_option_3months_1YR_put$`Payoff Realisasi`, cv_price_option_3months_1YR_put$`Put Jump`)
rmse_put_gbm_real_3months_1YR <- rmse(cv_price_option_3months_1YR_put$`Payoff Realisasi`, cv_price_option_3months_1YR_put$`Put GBM`)

#5YR
#Menghitung Logreturn
log_return_3months_5YR <- NULL
log_return_3months_5YR[1] <- 0
for(i in 2:nrow(data_3months_5YR)){
  log_return_3months_5YR[i] = log(data_3months_5YR$Last[i]/data_3months_5YR$Last[i-1])
}
data_3months_5YR$`log return` <- log_return_3months_5YR

#Volatilitas
vol_3months_5YR <- sd(data_3months_5YR$`log return`)

#Jump Parameter
is.jump_3months_5YR <- which(abs(data_3months_5YR$`log return`) > 2*sd(data_3months_5YR$`log return`))
num_jump_3months_5YR <- length(is.jump_3months_5YR)

mu_y_3months_5YR <- mean(data_3months_5YR$`log return`[is.jump_3months_5YR])
sigma_y_3months_5YR <- sd(data_3months_5YR$`log return`[is.jump_3months_5YR])
lam_3months_5YR <- num_jump_3months_5YR/nrow(data_3months_5YR)
So_3months_5YR <- data_3months_5YR$Last[nrow(data_3months_5YR)] 
sigma_3months_5YR <- sd(data_3months_5YR$`log return`)*sqrt(time_to_maturity_3months)

p_3months_5YR <- ggplot(data_3months_5YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_3months_5YR[is.jump_3months_5YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 29 Juni 2018 - 29 Juni 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_3months_5YR)

cv_price_option_jump_3months_5YR_call <- NULL
cv_price_option_jump_3months_5YR_put <- NULL
cv_price_option_gbm_3months_5YR_call <- NULL
cv_price_option_gbm_3months_5YR_put <- NULL

start_3months_5YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_3months_5YR_call <- rbind(cv_price_option_jump_3months_5YR_call, cv_price_jump_call(w, So = So_3months_5YR, K[i], r = r_3months, delta, sigma = sigma_3months_5YR, n = n_3months, j, 
                                                                                                           T = T_3months, lam = lam_3months_5YR, mu_y=mu_y_3months_5YR, sigma_y=sigma_y_3months_5YR, days, iter))
} 
end_3months_5YR_JUMP_call <- Sys.time()

start_3months_5YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_3months_5YR_put <- rbind(cv_price_option_jump_3months_5YR_put, cv_price_jump_put(w, So = So_3months_5YR, K[i], r = r_3months, delta, sigma = sigma_3months_5YR, n = n_3months, j, 
                                                                                                        T = T_3months, lam = lam_3months_5YR, mu_y=mu_y_3months_5YR, sigma_y=sigma_y_3months_5YR, days, iter))
}
end_3months_5YR_JUMP_put <- Sys.time()

start_3months_5YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_3months_5YR_call <- rbind(cv_price_option_gbm_3months_5YR_call, cv_price_gbm_call(w, So = So_3months_5YR, K[i], r = r_3months, delta, sigma = sigma_3months_5YR, n = n_3months, j, T = T_3months, days, iter))
  
}
end_3months_5YR_GBM_call <- Sys.time()

start_3months_5YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_3months_5YR_put <- rbind(cv_price_option_gbm_3months_5YR_put, cv_price_gbm_put(w, So = So_3months_5YR, K[i], r = r_3months, delta, sigma = sigma_3months_5YR, n = n_3months, j, T = T_3months, days, iter))
  
}
end_3months_5YR_GBM_put <- Sys.time()

time_3months_5YR_JUMP_call <- start_3months_5YR_JUMP_call - end_3months_5YR_JUMP_call
time_3months_5YR_JUMP_put <- start_3months_5YR_JUMP_put - end_3months_5YR_JUMP_put
time_3months_5YR_GBM_call <- start_3months_5YR_GBM_call - end_3months_5YR_GBM_call
time_3months_5YR_GBM_put <- start_3months_5YR_GBM_put - end_3months_5YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_3months_5YR_call <- cbind(cv_price_option_jump_3months_5YR_call$Call_Jump, cv_price_option_gbm_3months_5YR_call$Call_GBM)
cv_price_option_3months_5YR_call <- cbind(strike, cv_price_option_3months_5YR_call, real_price_3months$Call_RealPrice, real_payoff_call_3months)
colnames(cv_price_option_3months_5YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_3months_5YR_put <- cbind(cv_price_option_jump_3months_5YR_put$Put_Jump, cv_price_option_gbm_3months_5YR_put$Put_GBM)
cv_price_option_3months_5YR_put <- cbind(strike, cv_price_option_3months_5YR_put, real_price_3months$Put_RealPrice, real_payoff_put_3months)
colnames(cv_price_option_3months_5YR_put) <- col_names_put

rmse_call_jump_3months_5YR <- rmse(cv_price_option_3months_5YR_call$Barchart, cv_price_option_3months_5YR_call$`Call Jump`)
rmse_call_gbm_3months_5YR <- rmse(cv_price_option_3months_5YR_call$Barchart, cv_price_option_3months_5YR_call$`Call GBM`)
rmse_put_jump_3months_5YR <- rmse(cv_price_option_3months_5YR_put$Barchart, cv_price_option_3months_5YR_put$`Put Jump`)
rmse_put_gbm_3months_5YR <- rmse(cv_price_option_3months_5YR_put$Barchart, cv_price_option_3months_5YR_put$`Put GBM`)

rmse_call_jump_real_3months_5YR <- rmse(cv_price_option_3months_5YR_call$`Payoff Realisasi`, cv_price_option_3months_5YR_call$`Call Jump`)
rmse_call_gbm_real_3months_5YR <- rmse(cv_price_option_3months_5YR_call$`Payoff Realisasi`, cv_price_option_3months_5YR_call$`Call GBM`)
rmse_put_jump_real_3months_5YR <- rmse(cv_price_option_3months_5YR_put$`Payoff Realisasi`, cv_price_option_3months_5YR_put$`Put Jump`)
rmse_put_gbm_real_3months_5YR <- rmse(cv_price_option_3months_5YR_put$`Payoff Realisasi`, cv_price_option_3months_5YR_put$`Put GBM`)

#6 months
#data
data_asli_6months <- data[c(1347:1473),]
data_6months_1YR <- data[c(1095:1346),]
data_6months_5YR <- data[c(87:1346),]
At_6months <- mean(data_asli_6months$Last)

time_to_maturity_6months <- 184
T_6months <- time_to_maturity_6months/days  
i_6months <- 0.05552 #US Treasury Yield (6 month)
r_6months <- log(1+i_6months)
n_6months <- nrow(data_asli_6months) + 1

call_realprice_6months <- c(22.74, 14.95, 8.75, 4.44, 2.04, 0.95)
put_realprice_6months <- c(1.36, 3.33, 6.88, 12.33, 19.69, 28.36)
real_price_6months <- as.data.frame(list(Call_RealPrice = call_realprice_6months, Put_RealPrice = put_realprice_6months))

real_payoff_call_6months <- NULL
real_payoff_put_6months <- NULL
for (i in 1:length(K)) {
  real_payoff_call_6months <- c(real_payoff_call_6months, max(At_6months - K[i],0)*exp(-r_6months*T_6months))
  real_payoff_put_6months <- c(real_payoff_put_6months, max(K[i] - At_6months,0)*exp(-r_6months*T_6months))
} 

#1YR
#Menghitung Logreturn
log_return_6months_1YR <- NULL
log_return_6months_1YR[1] <- 0
for(i in 2:nrow(data_6months_1YR)){
  log_return_6months_1YR[i] = log(data_6months_1YR$Last[i]/data_6months_1YR$Last[i-1])
}
data_6months_1YR$`log return` <- log_return_6months_1YR

#Volatilitas
vol_6months_1YR <- sd(data_6months_1YR$`log return`)

#Jump Parameter
is.jump_6months_1YR <- which(abs(data_6months_1YR$`log return`) > 2*sd(data_6months_1YR$`log return`))
num_jump_6months_1YR <- length(is.jump_6months_1YR)

mu_y_6months_1YR <- mean(data_6months_1YR$`log return`[is.jump_6months_1YR])
sigma_y_6months_1YR <- sd(data_6months_1YR$`log return`[is.jump_6months_1YR])
lam_6months_1YR <- num_jump_6months_1YR/nrow(data_6months_1YR)
So_6months_1YR <- data_6months_1YR$Last[nrow(data_6months_1YR)] 
sigma_6months_1YR <- sd(data_6months_1YR$`log return`)*sqrt(time_to_maturity_6months)

p_6months_1YR <- ggplot(data_6months_1YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_6months_1YR[is.jump_6months_1YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 29 Maret 2022 - 29 Maret 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_6months_1YR)

cv_price_option_jump_6months_1YR_call <- NULL
cv_price_option_jump_6months_1YR_put <- NULL
cv_price_option_gbm_6months_1YR_call <- NULL
cv_price_option_gbm_6months_1YR_put <- NULL

start_6months_1YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_6months_1YR_call <- rbind(cv_price_option_jump_6months_1YR_call, cv_price_jump_call(w, So = So_6months_1YR, K[i], r = r_6months, delta, sigma = sigma_6months_1YR, n = n_6months, j, 
                                                                                                           T = T_6months, lam = lam_6months_1YR, mu_y=mu_y_6months_1YR, sigma_y=sigma_y_6months_1YR, days, iter))
} 
end_6months_1YR_JUMP_call <- Sys.time()

start_6months_1YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_6months_1YR_put <- rbind(cv_price_option_jump_6months_1YR_put, cv_price_jump_put(w, So = So_6months_1YR, K[i], r = r_6months, delta, sigma = sigma_6months_1YR, n = n_6months, j, 
                                                                                                        T = T_6months, lam = lam_6months_1YR, mu_y=mu_y_6months_1YR, sigma_y=sigma_y_6months_1YR, days, iter))
}
end_6months_1YR_JUMP_put <- Sys.time()

start_6months_1YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_6months_1YR_call <- rbind(cv_price_option_gbm_6months_1YR_call, cv_price_gbm_call(w, So = So_6months_1YR, K[i], r = r_6months, delta, sigma = sigma_6months_1YR, n = n_6months, j, T = T_6months, days, iter))
  
}
end_6months_1YR_GBM_call <- Sys.time()

start_6months_1YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_6months_1YR_put <- rbind(cv_price_option_gbm_6months_1YR_put, cv_price_gbm_put(w, So = So_6months_1YR, K[i], r = r_6months, delta, sigma = sigma_6months_1YR, n = n_6months, j, T = T_6months, days, iter))
  
}
end_6months_1YR_GBM_put <- Sys.time()

time_6months_1YR_JUMP_call <- start_6months_1YR_JUMP_call - end_6months_1YR_JUMP_call
time_6months_1YR_JUMP_put <- start_6months_1YR_JUMP_put - end_6months_1YR_JUMP_put
time_6months_1YR_GBM_call <- start_6months_1YR_GBM_call - end_6months_1YR_GBM_call
time_6months_1YR_GBM_put <- start_6months_1YR_GBM_put - end_6months_1YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_6months_1YR_call <- cbind(cv_price_option_jump_6months_1YR_call$Call_Jump, cv_price_option_gbm_6months_1YR_call$Call_GBM)
cv_price_option_6months_1YR_call <- cbind(strike, cv_price_option_6months_1YR_call, real_price_6months$Call_RealPrice, real_payoff_call_6months)
colnames(cv_price_option_6months_1YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_6months_1YR_put <- cbind(cv_price_option_jump_6months_1YR_put$Put_Jump, cv_price_option_gbm_6months_1YR_put$Put_GBM)
cv_price_option_6months_1YR_put <- cbind(strike, cv_price_option_6months_1YR_put, real_price_6months$Put_RealPrice, real_payoff_put_6months)
colnames(cv_price_option_6months_1YR_put) <- col_names_put

rmse_call_jump_6months_1YR <- rmse(cv_price_option_6months_1YR_call$Barchart, cv_price_option_6months_1YR_call$`Call Jump`)
rmse_call_gbm_6months_1YR <- rmse(cv_price_option_6months_1YR_call$Barchart, cv_price_option_6months_1YR_call$`Call GBM`)
rmse_put_jump_6months_1YR <- rmse(cv_price_option_6months_1YR_put$Barchart, cv_price_option_6months_1YR_put$`Put Jump`)
rmse_put_gbm_6months_1YR <- rmse(cv_price_option_6months_1YR_put$Barchart, cv_price_option_6months_1YR_put$`Put GBM`)

rmse_call_jump_real_6months_1YR <- rmse(cv_price_option_6months_1YR_call$`Payoff Realisasi`, cv_price_option_6months_1YR_call$`Call Jump`)
rmse_call_gbm_real_6months_1YR <- rmse(cv_price_option_6months_1YR_call$`Payoff Realisasi`, cv_price_option_6months_1YR_call$`Call GBM`)
rmse_put_jump_real_6months_1YR <- rmse(cv_price_option_6months_1YR_put$`Payoff Realisasi`, cv_price_option_6months_1YR_put$`Put Jump`)
rmse_put_gbm_real_6months_1YR <- rmse(cv_price_option_6months_1YR_put$`Payoff Realisasi`, cv_price_option_6months_1YR_put$`Put GBM`)

#5YR
#Menghitung Logreturn
log_return_6months_5YR <- NULL
log_return_6months_5YR[1] <- 0
for(i in 2:nrow(data_6months_5YR)){
  log_return_6months_5YR[i] = log(data_6months_5YR$Last[i]/data_6months_5YR$Last[i-1])
}
data_6months_5YR$`log return` <- log_return_6months_5YR

#Volatilitas
vol_6months_5YR <- sd(data_6months_5YR$`log return`)

#Jump Parameter
is.jump_6months_5YR <- which(abs(data_6months_5YR$`log return`) > 2*sd(data_6months_5YR$`log return`))
num_jump_6months_5YR <- length(is.jump_6months_5YR)

mu_y_6months_5YR <- mean(data_6months_5YR$`log return`[is.jump_6months_5YR])
sigma_y_6months_5YR <- sd(data_6months_5YR$`log return`[is.jump_6months_5YR])
lam_6months_5YR <- num_jump_6months_5YR/nrow(data_6months_5YR)
So_6months_5YR <- data_6months_5YR$Last[nrow(data_6months_5YR)] 
sigma_6months_5YR <- sd(data_6months_5YR$`log return`)*sqrt(time_to_maturity_6months)

p_6months_5YR <- ggplot(data_6months_5YR, aes(x = Time, y = Last)) +
  geom_line() +
  geom_point(data = data_6months_5YR[is.jump_6months_5YR, ], aes(x = Time, y = Last), color = "red", size = 2) +
  labs(x = "Tanggal", y = "Harga", title = "Pergerakan Harga Historis 29 Maret 2018 - 29 Maret 2023") +
  theme_minimal() +
  theme(plot.title = element_text(lineheight = 0.7, face = "bold"))
print(p_6months_5YR)

cv_price_option_jump_6months_5YR_call <- NULL
cv_price_option_jump_6months_5YR_put <- NULL
cv_price_option_gbm_6months_5YR_call <- NULL
cv_price_option_gbm_6months_5YR_put <- NULL

start_6months_5YR_JUMP_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_6months_5YR_call <- rbind(cv_price_option_jump_6months_5YR_call, cv_price_jump_call(w, So = So_6months_5YR, K[i], r = r_6months, delta, sigma = sigma_6months_5YR, n = n_6months, j, 
                                                                                                           T = T_6months, lam = lam_6months_5YR, mu_y=mu_y_6months_5YR, sigma_y=sigma_y_6months_5YR, days, iter))
} 
end_6months_5YR_JUMP_call <- Sys.time()

start_6months_5YR_JUMP_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_jump_6months_5YR_put <- rbind(cv_price_option_jump_6months_5YR_put, cv_price_jump_put(w, So = So_6months_5YR, K[i], r = r_6months, delta, sigma = sigma_6months_5YR, n = n_6months, j, 
                                                                                                        T = T_6months, lam = lam_6months_5YR, mu_y=mu_y_6months_5YR, sigma_y=sigma_y_6months_5YR, days, iter))
}
end_6months_5YR_JUMP_put <- Sys.time()

start_6months_5YR_GBM_call <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_6months_5YR_call <- rbind(cv_price_option_gbm_6months_5YR_call, cv_price_gbm_call(w, So = So_6months_5YR, K[i], r = r_6months, delta, sigma = sigma_6months_5YR, n = n_6months, j, T = T_6months, days, iter))
  
}
end_6months_5YR_GBM_call <- Sys.time()

start_6months_5YR_GBM_put <- Sys.time()
for (i in 1:length(K)) {
  cv_price_option_gbm_6months_5YR_put <- rbind(cv_price_option_gbm_6months_5YR_put, cv_price_gbm_put(w, So = So_6months_5YR, K[i], r = r_6months, delta, sigma = sigma_6months_5YR, n = n_6months, j, T = T_6months, days, iter))
  
}
end_6months_5YR_GBM_put <- Sys.time()

time_6months_5YR_JUMP_call <- start_6months_5YR_JUMP_call - end_6months_5YR_JUMP_call
time_6months_5YR_JUMP_put <- start_6months_5YR_JUMP_put - end_6months_5YR_JUMP_put
time_6months_5YR_GBM_call <- start_6months_5YR_GBM_call - end_6months_5YR_GBM_call
time_6months_5YR_GBM_put <- start_6months_5YR_GBM_put - end_6months_5YR_GBM_put

col_names_call <- c("K", "Call Jump", "Call GBM", "Barchart", "Payoff Realisasi")
cv_price_option_6months_5YR_call <- cbind(cv_price_option_jump_6months_5YR_call$Call_Jump, cv_price_option_gbm_6months_5YR_call$Call_GBM)
cv_price_option_6months_5YR_call <- cbind(strike, cv_price_option_6months_5YR_call, real_price_6months$Call_RealPrice, real_payoff_call_6months)
colnames(cv_price_option_6months_5YR_call) <- col_names_call

col_names_put <- c("K", "Put Jump", "Put GBM", "Barchart", "Payoff Realisasi")
cv_price_option_6months_5YR_put <- cbind(cv_price_option_jump_6months_5YR_put$Put_Jump, cv_price_option_gbm_6months_5YR_put$Put_GBM)
cv_price_option_6months_5YR_put <- cbind(strike, cv_price_option_6months_5YR_put, real_price_6months$Put_RealPrice, real_payoff_put_6months)
colnames(cv_price_option_6months_5YR_put) <- col_names_put

rmse_call_jump_6months_5YR <- rmse(cv_price_option_6months_5YR_call$Barchart, cv_price_option_6months_5YR_call$`Call Jump`)
rmse_call_gbm_6months_5YR <- rmse(cv_price_option_6months_5YR_call$Barchart, cv_price_option_6months_5YR_call$`Call GBM`)
rmse_put_jump_6months_5YR <- rmse(cv_price_option_6months_5YR_put$Barchart, cv_price_option_6months_5YR_put$`Put Jump`)
rmse_put_gbm_6months_5YR <- rmse(cv_price_option_6months_5YR_put$Barchart, cv_price_option_6months_5YR_put$`Put GBM`)

rmse_call_jump_real_6months_5YR <- rmse(cv_price_option_6months_5YR_call$`Payoff Realisasi`, cv_price_option_6months_5YR_call$`Call Jump`)
rmse_call_gbm_real_6months_5YR <- rmse(cv_price_option_6months_5YR_call$`Payoff Realisasi`, cv_price_option_6months_5YR_call$`Call GBM`)
rmse_put_jump_real_6months_5YR <- rmse(cv_price_option_6months_5YR_put$`Payoff Realisasi`, cv_price_option_6months_5YR_put$`Put Jump`)
rmse_put_gbm_real_6months_5YR <- rmse(cv_price_option_6months_5YR_put$`Payoff Realisasi`, cv_price_option_6months_5YR_put$`Put GBM`)

#DATA FINAL
simulation_type <- c("1YR_10000_Call_16Hari",
                     "1YR_10000_Put_16Hari",
                     "1YR_10000_Call_1Bulan",
                     "1YR_10000_Put_1Bulan",
                     "1YR_10000_Call_3Bulan",
                     "1YR_10000_Put_3Bulan",
                     "1YR_10000_Call_6Bulan",
                     "1YR_10000_Put_6Bulan",
                     "5YR_10000_Call_16Hari",
                     "5YR_10000_Put_16Hari",
                     "5YR_10000_Call_1Bulan",
                     "5YR_10000_Put_1Bulan",
                     "5YR_10000_Call_3Bulan",
                     "5YR_10000_Put_3Bulan",
                     "5YR_10000_Call_6Bulan",
                     "5YR_10000_Put_6Bulan")

jump_error <- c(rmse_call_jump_16days_1YR,
                rmse_put_jump_16days_1YR,
                rmse_call_jump_1month_1YR,
                rmse_put_jump_1month_1YR,
                rmse_call_jump_3months_1YR,
                rmse_put_jump_3months_1YR,
                rmse_call_jump_6months_1YR,
                rmse_put_jump_6months_1YR,
                rmse_call_jump_16days_5YR,
                rmse_put_jump_16days_5YR,
                rmse_call_jump_1month_5YR,
                rmse_put_jump_1month_5YR,
                rmse_call_jump_3months_5YR,
                rmse_put_jump_3months_5YR,
                rmse_call_jump_6months_5YR,
                rmse_put_jump_6months_5YR)

jump_time <- c(time_16days_1YR_JUMP_call, time_16days_1YR_JUMP_put,
               time_1month_1YR_JUMP_call, time_1month_1YR_JUMP_put,
               time_3months_1YR_JUMP_call, time_3months_1YR_JUMP_put,
               time_6months_1YR_JUMP_call, time_6months_1YR_JUMP_put,
               time_16days_5YR_JUMP_call, time_16days_5YR_JUMP_put,
               time_1month_5YR_JUMP_call, time_1month_5YR_JUMP_put,
               time_3months_5YR_JUMP_call, time_3months_5YR_JUMP_put,
               time_6months_5YR_JUMP_call, time_6months_5YR_JUMP_put)
jump_time <- abs(jump_time)

GBM_error <- c(rmse_call_gbm_16days_1YR,
               rmse_put_gbm_16days_1YR,
               rmse_call_gbm_1month_1YR,
               rmse_put_gbm_1month_1YR,
               rmse_call_gbm_3months_1YR,
               rmse_put_gbm_3months_1YR,
               rmse_call_gbm_6months_1YR,
               rmse_put_gbm_6months_1YR,
               rmse_call_gbm_16days_5YR,
               rmse_put_gbm_16days_5YR,
               rmse_call_gbm_1month_5YR,
               rmse_put_gbm_1month_5YR,
               rmse_call_gbm_3months_5YR,
               rmse_put_gbm_3months_5YR,
               rmse_call_gbm_6months_5YR,
               rmse_put_gbm_6months_5YR)

GBM_time <- c(time_16days_1YR_GBM_call, time_16days_1YR_GBM_put,
              time_1month_1YR_GBM_call, time_1month_1YR_GBM_put,
              time_3months_1YR_GBM_call, time_3months_1YR_GBM_put,
              time_6months_1YR_GBM_call, time_6months_1YR_GBM_put,
              time_16days_5YR_GBM_call, time_16days_5YR_GBM_put,
              time_1month_5YR_GBM_call, time_1month_5YR_GBM_put,
              time_3months_5YR_GBM_call, time_3months_5YR_GBM_put,
              time_6months_5YR_GBM_call, time_6months_5YR_GBM_put)
GBM_time <- abs(GBM_time)

error <- as.data.frame(list("sim_type" = simulation_type, "jump_error" = jump_error,
                            "jump_time" = jump_time, "gbm_error" = GBM_error, "gbm_time" = GBM_time))

jump_error_real <- c(rmse_call_jump_real_16days_1YR,
                     rmse_put_jump_real_16days_1YR,
                     rmse_call_jump_real_1month_1YR,
                     rmse_put_jump_real_1month_1YR,
                     rmse_call_jump_real_3months_1YR,
                     rmse_put_jump_real_3months_1YR,
                     rmse_call_jump_real_6months_1YR,
                     rmse_put_jump_real_6months_1YR,
                     rmse_call_jump_real_16days_5YR,
                     rmse_put_jump_real_16days_5YR,
                     rmse_call_jump_real_1month_5YR,
                     rmse_put_jump_real_1month_5YR,
                     rmse_call_jump_real_3months_5YR,
                     rmse_put_jump_real_3months_5YR,
                     rmse_call_jump_real_6months_5YR,
                     rmse_put_jump_real_6months_5YR)

gbm_error_real <- c(rmse_call_gbm_real_16days_1YR,
                    rmse_put_gbm_real_16days_1YR,
                    rmse_call_gbm_real_1month_1YR,
                    rmse_put_gbm_real_1month_1YR,
                    rmse_call_gbm_real_3months_1YR,
                    rmse_put_gbm_real_3months_1YR,
                    rmse_call_gbm_real_6months_1YR,
                    rmse_put_gbm_real_6months_1YR,
                    rmse_call_gbm_real_16days_5YR,
                    rmse_put_gbm_real_16days_5YR,
                    rmse_call_gbm_real_1month_5YR,
                    rmse_put_gbm_real_1month_5YR,
                    rmse_call_gbm_real_3months_5YR,
                    rmse_put_gbm_real_3months_5YR,
                    rmse_call_gbm_real_6months_5YR,
                    rmse_put_gbm_real_6months_5YR)

real_error <- as.data.frame(list("sim_type" = simulation_type, "jump_error" = jump_error_real,
                                 "jump_time" = jump_time, "gbm_error" = gbm_error_real, "gbm_time" = GBM_time))

#Moneyness
#16 Hari Call
strike_50_16days_call <- rbind(cv_price_option_16days_1YR_call[1,], 
                               setNames(cv_price_option_16days_5YR_call[1,], names(cv_price_option_16days_1YR_call[1,])))
strike_60_16days_call <- rbind(cv_price_option_16days_1YR_call[2,],
                               setNames(cv_price_option_16days_5YR_call[2,], names(cv_price_option_16days_1YR_call[2,])))
strike_70_16days_call <- rbind(cv_price_option_16days_1YR_call[3,],
                               setNames(cv_price_option_16days_5YR_call[3,], names(cv_price_option_16days_1YR_call[3,])))
strike_80_16days_call <- rbind(cv_price_option_16days_1YR_call[4,],
                               setNames(cv_price_option_16days_5YR_call[4,], names(cv_price_option_16days_1YR_call[4,])))
strike_90_16days_call <- rbind(cv_price_option_16days_1YR_call[5,],
                               setNames(cv_price_option_16days_5YR_call[5,], names(cv_price_option_16days_1YR_call[5,])))
strike_100_16days_call <- rbind(cv_price_option_16days_1YR_call[6,],
                                setNames(cv_price_option_16days_5YR_call[6,], names(cv_price_option_16days_1YR_call[6,])))

#16 Hari put
strike_50_16days_put <- rbind(cv_price_option_16days_1YR_put[1,], 
                              setNames(cv_price_option_16days_5YR_put[1,], names(cv_price_option_16days_1YR_put[1,])))
strike_60_16days_put <- rbind(cv_price_option_16days_1YR_put[2,],
                              setNames(cv_price_option_16days_5YR_put[2,], names(cv_price_option_16days_1YR_put[2,])))
strike_70_16days_put <- rbind(cv_price_option_16days_1YR_put[3,],
                              setNames(cv_price_option_16days_5YR_put[3,], names(cv_price_option_16days_1YR_put[3,])))
strike_80_16days_put <- rbind(cv_price_option_16days_1YR_put[4,],
                              setNames(cv_price_option_16days_5YR_put[4,], names(cv_price_option_16days_1YR_put[4,])))
strike_90_16days_put <- rbind(cv_price_option_16days_1YR_put[5,],
                              setNames(cv_price_option_16days_5YR_put[5,], names(cv_price_option_16days_1YR_put[5,])))
strike_100_16days_put <- rbind(cv_price_option_16days_1YR_put[6,],
                               setNames(cv_price_option_16days_5YR_put[6,], names(cv_price_option_16days_1YR_put[6,])))

#1 Month Call
strike_50_1month_call <- rbind(cv_price_option_1month_1YR_call[1,], 
                               setNames(cv_price_option_1month_5YR_call[1,], names(cv_price_option_1month_1YR_call[1,])))
strike_60_1month_call <- rbind(cv_price_option_1month_1YR_call[2,],
                               setNames(cv_price_option_1month_5YR_call[2,], names(cv_price_option_1month_1YR_call[2,])))
strike_70_1month_call <- rbind(cv_price_option_1month_1YR_call[3,],
                               setNames(cv_price_option_1month_5YR_call[3,], names(cv_price_option_1month_1YR_call[3,])))
strike_80_1month_call <- rbind(cv_price_option_1month_1YR_call[4,],
                               setNames(cv_price_option_1month_5YR_call[4,], names(cv_price_option_1month_1YR_call[4,])))
strike_90_1month_call <- rbind(cv_price_option_1month_1YR_call[5,],
                               setNames(cv_price_option_1month_5YR_call[5,], names(cv_price_option_1month_1YR_call[5,])))
strike_100_1month_call <- rbind(cv_price_option_1month_1YR_call[6,],
                                setNames(cv_price_option_1month_5YR_call[6,], names(cv_price_option_1month_1YR_call[6,])))

#1 Month put
strike_50_1month_put <- rbind(cv_price_option_1month_1YR_put[1,], 
                              setNames(cv_price_option_1month_5YR_put[1,], names(cv_price_option_1month_1YR_put[1,])))
strike_60_1month_put <- rbind(cv_price_option_1month_1YR_put[2,],
                              setNames(cv_price_option_1month_5YR_put[2,], names(cv_price_option_1month_1YR_put[2,])))
strike_70_1month_put <- rbind(cv_price_option_1month_1YR_put[3,],
                              setNames(cv_price_option_1month_5YR_put[3,], names(cv_price_option_1month_1YR_put[3,])))
strike_80_1month_put <- rbind(cv_price_option_1month_1YR_put[4,],
                              setNames(cv_price_option_1month_5YR_put[4,], names(cv_price_option_1month_1YR_put[4,])))
strike_90_1month_put <- rbind(cv_price_option_1month_1YR_put[5,],
                              setNames(cv_price_option_1month_5YR_put[5,], names(cv_price_option_1month_1YR_put[5,])))
strike_100_1month_put <- rbind(cv_price_option_1month_1YR_put[6,],
                               setNames(cv_price_option_1month_5YR_put[6,], names(cv_price_option_1month_1YR_put[6,])))

#3 Months Call
strike_50_3months_call <- rbind(cv_price_option_3months_1YR_call[1,], 
                                setNames(cv_price_option_3months_5YR_call[1,], names(cv_price_option_3months_1YR_call[1,])))
strike_60_3months_call <- rbind(cv_price_option_3months_1YR_call[2,],
                                setNames(cv_price_option_3months_5YR_call[2,], names(cv_price_option_3months_1YR_call[2,])))
strike_70_3months_call <- rbind(cv_price_option_3months_1YR_call[3,],
                                setNames(cv_price_option_3months_5YR_call[3,], names(cv_price_option_3months_1YR_call[3,])))
strike_80_3months_call <- rbind(cv_price_option_3months_1YR_call[4,],
                                setNames(cv_price_option_3months_5YR_call[4,], names(cv_price_option_3months_1YR_call[4,])))
strike_90_3months_call <- rbind(cv_price_option_3months_1YR_call[5,],
                                setNames(cv_price_option_3months_5YR_call[5,], names(cv_price_option_3months_1YR_call[5,])))
strike_100_3months_call <- rbind(cv_price_option_3months_1YR_call[6,],
                                 setNames(cv_price_option_3months_5YR_call[6,], names(cv_price_option_3months_1YR_call[6,])))

#3 Months put
strike_50_3months_put <- rbind(cv_price_option_3months_1YR_put[1,], 
                               setNames(cv_price_option_3months_5YR_put[1,], names(cv_price_option_3months_1YR_put[1,])))
strike_60_3months_put <- rbind(cv_price_option_3months_1YR_put[2,],
                               setNames(cv_price_option_3months_5YR_put[2,], names(cv_price_option_3months_1YR_put[2,])))
strike_70_3months_put <- rbind(cv_price_option_3months_1YR_put[3,],
                               setNames(cv_price_option_3months_5YR_put[3,], names(cv_price_option_3months_1YR_put[3,])))
strike_80_3months_put <- rbind(cv_price_option_3months_1YR_put[4,],
                               setNames(cv_price_option_3months_5YR_put[4,], names(cv_price_option_3months_1YR_put[4,])))
strike_90_3months_put <- rbind(cv_price_option_3months_1YR_put[5,],
                               setNames(cv_price_option_3months_5YR_put[5,], names(cv_price_option_3months_1YR_put[5,])))
strike_100_3months_put <- rbind(cv_price_option_3months_1YR_put[6,],
                                setNames(cv_price_option_3months_5YR_put[6,], names(cv_price_option_3months_1YR_put[6,])))

#6 Months Call
strike_50_6months_call <- rbind(cv_price_option_6months_1YR_call[1,], 
                                setNames(cv_price_option_6months_5YR_call[1,], names(cv_price_option_6months_1YR_call[1,])))
strike_60_6months_call <- rbind(cv_price_option_6months_1YR_call[2,],
                                setNames(cv_price_option_6months_5YR_call[2,], names(cv_price_option_6months_1YR_call[2,])))
strike_70_6months_call <- rbind(cv_price_option_6months_1YR_call[3,],
                                setNames(cv_price_option_6months_5YR_call[3,], names(cv_price_option_6months_1YR_call[3,])))
strike_80_6months_call <- rbind(cv_price_option_6months_1YR_call[4,],
                                setNames(cv_price_option_6months_5YR_call[4,], names(cv_price_option_6months_1YR_call[4,])))
strike_90_6months_call <- rbind(cv_price_option_6months_1YR_call[5,],
                                setNames(cv_price_option_6months_5YR_call[5,], names(cv_price_option_6months_1YR_call[5,])))
strike_100_6months_call <- rbind(cv_price_option_6months_1YR_call[6,],
                                 setNames(cv_price_option_6months_5YR_call[6,], names(cv_price_option_6months_1YR_call[6,])))

#6 Months put
strike_50_6months_put <- rbind(cv_price_option_6months_1YR_put[1,], 
                               setNames(cv_price_option_6months_5YR_put[1,], names(cv_price_option_6months_1YR_put[1,])))
strike_60_6months_put <- rbind(cv_price_option_6months_1YR_put[2,],
                               setNames(cv_price_option_6months_5YR_put[2,], names(cv_price_option_6months_1YR_put[2,])))
strike_70_6months_put <- rbind(cv_price_option_6months_1YR_put[3,],
                               setNames(cv_price_option_6months_5YR_put[3,], names(cv_price_option_6months_1YR_put[3,])))
strike_80_6months_put <- rbind(cv_price_option_6months_1YR_put[4,],
                               setNames(cv_price_option_6months_5YR_put[4,], names(cv_price_option_6months_1YR_put[4,])))
strike_90_6months_put <- rbind(cv_price_option_6months_1YR_put[5,],
                               setNames(cv_price_option_6months_5YR_put[5,], names(cv_price_option_6months_1YR_put[5,])))
strike_100_6months_put <- rbind(cv_price_option_6months_1YR_put[6,],
                                setNames(cv_price_option_6months_5YR_put[6,], names(cv_price_option_6months_1YR_put[6,])))

#cor
cor_16days_1YR_call <- cbind(strike, cv_price_option_jump_16days_1YR_call$cor, cv_price_option_gbm_16days_1YR_call$cor)
cor_16days_5YR_call <- cbind(strike, cv_price_option_jump_16days_5YR_call$cor, cv_price_option_gbm_16days_5YR_call$cor)
cor_1month_1YR_call <- cbind(strike, cv_price_option_jump_1month_1YR_call$cor, cv_price_option_gbm_1month_1YR_call$cor)
cor_1month_5YR_call <- cbind(strike, cv_price_option_jump_1month_5YR_call$cor, cv_price_option_gbm_1month_5YR_call$cor)
cor_3months_1YR_call <- cbind(strike, cv_price_option_jump_3months_1YR_call$cor, cv_price_option_gbm_3months_1YR_call$cor)
cor_3months_5YR_call <- cbind(strike, cv_price_option_jump_3months_5YR_call$cor, cv_price_option_gbm_3months_5YR_call$cor)
cor_6months_1YR_call <- cbind(strike, cv_price_option_jump_6months_1YR_call$cor, cv_price_option_gbm_6months_1YR_call$cor)
cor_6months_5YR_call <- cbind(strike, cv_price_option_jump_6months_5YR_call$cor, cv_price_option_gbm_6months_5YR_call$cor)

cor_16days_1YR_put <- cbind(strike, cv_price_option_jump_16days_1YR_put$cor, cv_price_option_gbm_16days_1YR_put$cor)
cor_16days_5YR_put <- cbind(strike, cv_price_option_jump_16days_5YR_put$cor, cv_price_option_gbm_16days_5YR_put$cor)
cor_1month_1YR_put <- cbind(strike, cv_price_option_jump_1month_1YR_put$cor, cv_price_option_gbm_1month_1YR_put$cor)
cor_1month_5YR_put <- cbind(strike, cv_price_option_jump_1month_5YR_put$cor, cv_price_option_gbm_1month_5YR_put$cor)
cor_3months_1YR_put <- cbind(strike, cv_price_option_jump_3months_1YR_put$cor, cv_price_option_gbm_3months_1YR_put$cor)
cor_3months_5YR_put <- cbind(strike, cv_price_option_jump_3months_5YR_put$cor, cv_price_option_gbm_3months_5YR_put$cor)
cor_6months_1YR_put <- cbind(strike, cv_price_option_jump_6months_1YR_put$cor, cv_price_option_gbm_6months_1YR_put$cor)
cor_6months_5YR_put <- cbind(strike, cv_price_option_jump_6months_5YR_put$cor, cv_price_option_gbm_6months_5YR_put$cor)

library(dplyr)
combined_data_cor <- bind_rows(
  cor_16days_1YR_call,
  cor_16days_5YR_call,
  cor_1month_1YR_call,
  cor_1month_5YR_call,
  cor_3months_1YR_call,
  cor_3months_5YR_call,
  cor_6months_1YR_call,
  cor_6months_5YR_call,
  cor_16days_1YR_put,
  cor_16days_5YR_put,
  cor_1month_1YR_put,
  cor_1month_5YR_put,
  cor_3months_1YR_put,
  cor_3months_5YR_put,
  cor_6months_1YR_put,
  cor_6months_5YR_put
)

#Std Errorr
std_err_16days_1YR_call <- cbind(strike, cv_price_option_jump_16days_1YR_call$std_err, cv_price_option_gbm_16days_1YR_call$std_err)
std_err_16days_5YR_call <- cbind(strike, cv_price_option_jump_16days_5YR_call$std_err, cv_price_option_gbm_16days_5YR_call$std_err)
std_err_1month_1YR_call <- cbind(strike, cv_price_option_jump_1month_1YR_call$std_err, cv_price_option_gbm_1month_1YR_call$std_err)
std_err_1month_5YR_call <- cbind(strike, cv_price_option_jump_1month_5YR_call$std_err, cv_price_option_gbm_1month_5YR_call$std_err)
std_err_3months_1YR_call <- cbind(strike, cv_price_option_jump_3months_1YR_call$std_err, cv_price_option_gbm_3months_1YR_call$std_err)
std_err_3months_5YR_call <- cbind(strike, cv_price_option_jump_3months_5YR_call$std_err, cv_price_option_gbm_3months_5YR_call$std_err)
std_err_6months_1YR_call <- cbind(strike, cv_price_option_jump_6months_1YR_call$std_err, cv_price_option_gbm_6months_1YR_call$std_err)
std_err_6months_5YR_call <- cbind(strike, cv_price_option_jump_6months_5YR_call$std_err, cv_price_option_gbm_6months_5YR_call$std_err)

std_err_16days_1YR_put <- cbind(strike, cv_price_option_jump_16days_1YR_put$std_err, cv_price_option_gbm_16days_1YR_put$std_err)
std_err_16days_5YR_put <- cbind(strike, cv_price_option_jump_16days_5YR_put$std_err, cv_price_option_gbm_16days_5YR_put$std_err)
std_err_1month_1YR_put <- cbind(strike, cv_price_option_jump_1month_1YR_put$std_err, cv_price_option_gbm_1month_1YR_put$std_err)
std_err_1month_5YR_put <- cbind(strike, cv_price_option_jump_1month_5YR_put$std_err, cv_price_option_gbm_1month_5YR_put$std_err)
std_err_3months_1YR_put <- cbind(strike, cv_price_option_jump_3months_1YR_put$std_err, cv_price_option_gbm_3months_1YR_put$std_err)
std_err_3months_5YR_put <- cbind(strike, cv_price_option_jump_3months_5YR_put$std_err, cv_price_option_gbm_3months_5YR_put$std_err)
std_err_6months_1YR_put <- cbind(strike, cv_price_option_jump_6months_1YR_put$std_err, cv_price_option_gbm_6months_1YR_put$std_err)
std_err_6months_5YR_put <- cbind(strike, cv_price_option_jump_6months_5YR_put$std_err, cv_price_option_gbm_6months_5YR_put$std_err)

combined_data_stder <- bind_rows(
  std_err_16days_1YR_call,
  std_err_16days_5YR_call,
  std_err_1month_1YR_call,
  std_err_1month_5YR_call,
  std_err_3months_1YR_call,
  std_err_3months_5YR_call,
  std_err_6months_1YR_call,
  std_err_6months_5YR_call,
  std_err_16days_1YR_put,
  std_err_16days_5YR_put,
  std_err_1month_1YR_put,
  std_err_1month_5YR_put,
  std_err_3months_1YR_put,
  std_err_3months_5YR_put,
  std_err_6months_1YR_put,
  std_err_6months_5YR_put
)



