# Load necessary libraries
library(data.table)
library(tidyverse)
library(lubridate)
library(magrittr)
library(dyn)
library(reshape2)
library(formattable)
library(sparkline)
library(xtable)
library(DT)
library(htmltools)
library(dplyr)
library(tidyr)
library(kableExtra)
library(tibble)
library(formattable)
library(DT)
library(htmltools)
library(scales)
library(lmtest)
library(sandwich)

# Import the data
annual <- fread(file.path("/Users/tolgacanbaybes/Desktop/Goyal_welch_replications/PredictorData2023(goyalwebsite).csv"))
table1_2005_g_w <- fread(file.path("/Users/tolgacanbaybes/Desktop/Goyal_welch_replications/annual_2005_goyal_welch.csv"))

# Define variables
### Default yield spread (dfy)
annual <- annual[, dfy := BAA - AAA]

### Default Return Spread (dfr)
annual <- annual[, dfr := corpr - ltr]

### Stock returns (IndexDiv)
annual$Index <- as.numeric(annual$Index)
annual <- annual[, IndexDiv := Index + D12]

### Dividend-price ratio (dp)
annual <- annual[, dp := log(D12) - log(Index)]

### Log 3-month treasury bill (logRfree)
annual <- annual[, logRfree := log(Rfree + 1)]

### Term spread (tms)
annual <- annual[, tms := lty - tbl]

### Dividend yield vector (vec_dy)
vec_dy <- c(NA, annual[2:nrow(annual), log(D12)] - annual[1:(nrow(annual)-1), log(Index)])

### Dividend yield (dy)
annual <- annual[, dy := vec_dy]

### Log return (logret)
annual <- annual[, logret := c(NA, diff(log(Index)))]

### Log return with dividends (logretdiv)
vec_logretdiv <- c(NA, annual[2:nrow(annual), log(IndexDiv)] - annual[1:(nrow(annual)-1), log(Index)])
annual <- annual[, logretdiv := vec_logretdiv]

### Earnings-price ratio (ep)
annual <- annual[, ep := log(E12) - log(Index)]

### Dividend payout ratio (de)
annual <- annual[, de := log(D12) - log(E12)]

### Log equity premium (rp_div)
annual <- annual[, rp_div := logretdiv - logRfree]

### News
annual <- annual[, news := news]

### VIX
annual <- annual[, vix := vix]

# Create time series object
ts_annual <- ts(annual, start=annual[1, yyyy], end=annual[nrow(annual), yyyy])

# Function to get statistics and test for heteroscedasticity
get_statistics2 <- function(ts_df, indep, dep, h=1, start=1872, end=2023, est_periods_OOS = 20) {
  # In Sample
  avg <- mean(window(ts_df, start, end)[, dep], na.rm=TRUE)
  IS_error_N <- (window(ts_df, start, end)[, dep] - avg)
  reg <- dyn$lm(eval(parse(text=dep)) ~ stats::lag(eval(parse(text=indep)), -1), data=window(ts_df, start, end))
  IS_error_A <- reg$residuals
  
  # Out Of Sample
  OOS_error_N <- numeric(end - start - est_periods_OOS)
  OOS_error_A <- numeric(end - start - est_periods_OOS)
  j <- 0
  
  for (i in (start + est_periods_OOS):(end-1)) {
    j <- j + 1
    actual_ERP <- as.numeric(window(ts_df, i+1, i+1)[, dep])
    OOS_error_N[j] <- actual_ERP - mean(window(ts_df, start, i)[, dep], na.rm=TRUE)
    reg_OOS <- dyn$lm(eval(parse(text=dep)) ~ stats::lag(eval(parse(text=indep)), -1), data=window(ts_df, start, i))
    df <- data.frame(x=as.numeric(window(ts_df, i, i)[, indep]))
    names(df) <- indep
    pred_ERP <- predict.lm(reg_OOS, newdata=df)
    OOS_error_A[j] <-  pred_ERP - actual_ERP
  }
  
  MSE_N <- mean(OOS_error_N^2)
  MSE_A <- mean(OOS_error_A^2)
  OOS_R2 <- 1 - MSE_A/MSE_N
  OOS_aR2 <- 1 - (((1-OOS_R2)*(reg_OOS$df.residual))/(reg_OOS$df.residual - 1))
  dRMSE <- sqrt(MSE_N) - sqrt(MSE_A)
  
  # Residual plot
  plot(fitted(reg), resid(reg), main=paste("Residuals vs Fitted for", indep), xlab="Fitted values", ylab="Residuals")
  abline(h = 0, col = "red")
  
  # Breusch-Pagan Test
  bp_test <- bptest(reg)
  
  # Robust standard errors
  robust_se <- coeftest(reg, vcov = vcovHC(reg, type = "HC1"))
  
  return(list(
    IS_error_N = IS_error_N,
    IS_error_A = reg$residuals,
    OOS_error_N = OOS_error_N,
    OOS_error_A = OOS_error_A,
    IS_R2 = summary(reg)$r.squared, 
    IS_aR2 = summary(reg)$adj.r.squared, 
    OOS_R2 = OOS_R2,
    OOS_aR2 = OOS_aR2,
    dRMSE = dRMSE,
    plotGG = plotGG,
    bp_test = bp_test,
    robust_se = robust_se
  ))
}

# Perform analysis for each predictor
dp_stat <- get_statistics2(ts_annual, "dp", "rp_div", start=1872)
dy_stat <- get_statistics2(ts_annual, "dy", "rp_div", start=1872)
ep_stat <- get_statistics2(ts_annual, "ep", "rp_div", start=1872)
de_stat <- get_statistics2(ts_annual, "de", "rp_div", start=1872)
svar_stat <- get_statistics2(ts_annual, "svar", "rp_div", start=1885)
tbl_stat <- get_statistics2(ts_annual, "tbl", "rp_div", start=1920)
lty_stat <- get_statistics2(ts_annual, "lty", "rp_div", start=1919)
ltr_stat <- get_statistics2(ts_annual, "ltr", "rp_div", start=1926)
tms_stat <- get_statistics2(ts_annual, "tms", "rp_div", start=1920)
dfy_stat <- get_statistics2(ts_annual, "dfy", "rp_div", start=1919)
dfr_stat <- get_statistics2(ts_annual, "dfr", "rp_div", start=1926)
infl_stat <- get_statistics2(ts_annual, "infl", "rp_div", start=1919)
bm_stat <- get_statistics2(ts_annual, "bm", "rp_div", start=1921)
ik_stat <- get_statistics2(ts_annual, "ik", "rp_div", start=1947)
ntis_stat <- get_statistics2(ts_annual, "ntis", "rp_div", start=1927)
eqis_stat <- get_statistics2(ts_annual, "eqis", "rp_div", start=1927)
news_stat <- get_statistics2(ts_annual, "news", "rp_div", start=1980)
vix_stat <- get_statistics2(ts_annual, "vix", "rp_div", start=1990)

# Print Breusch-Pagan test results and robust standard errors for each predictor
print(dp_stat$bp_test)
print(dp_stat$robust_se)

print(dy_stat$bp_test)
print(dy_stat$robust_se)

print(ep_stat$bp_test)
print(ep_stat$robust_se)

print(de_stat$bp_test)
print(de_stat$robust_se)

print(svar_stat$bp_test)
print(svar_stat$robust_se)

print(tbl_stat$bp_test)
print(tbl_stat$robust_se)

print(lty_stat$bp_test)
print(lty_stat$robust_se)

print(ltr_stat$bp_test)
print(ltr_stat$robust_se)

print(tms_stat$bp_test)
print(tms_stat$robust_se)

print(dfy_stat$bp_test)
print(dfy_stat$robust_se)

print(dfr_stat$bp_test)
print(dfr_stat$robust_se)

print(infl_stat$bp_test)
print(infl_stat$robust_se)

print(bm_stat$bp_test)
print(bm_stat$robust_se)

print(ik_stat$bp_test)
print(ik_stat$robust_se)

print(ntis_stat$bp_test)
print(ntis_stat$robust_se)

print(eqis_stat$bp_test)
print(eqis_stat$robust_se)

print(news_stat$bp_test)
print(news_stat$robust_se)

print(vix_stat$bp_test)
print(vix_stat$robust_se)
