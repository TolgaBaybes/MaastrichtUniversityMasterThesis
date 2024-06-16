#importing the libraries

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
library(strucchange)


#importing the data

annual <- fread(file.path("/Users/tolgacanbaybes/Desktop/Goyal_welch_replications/PredictorData2023(goyalwebsite).csv"))

table1_2005_g_w <- fread(file.path("/Users/tolgacanbaybes/Desktop/Goyal_welch_replications/annual_2005_goyal_welch.csv"))

#define variables 

#Default yield spread (dfy)

annual <- annual[, dfy := BAA - AAA]

#Default Return Spread (dfr)

annual <- annual[, dfr := corpr - ltr] 

##  Stock returns (IndexDiv)                                                ####
##  ............................................................................
##  continuously compounded S&P 500 returns (Index) + 12-month moving sum of
##  S&P Dividends (D12)

annual$Index <- as.numeric(annual$Index) # In 2018 dataset Index is set to 'chr' datatype 
annual <- annual[, IndexDiv := Index + D12]

##  Dividend-price ratio (dp)      

annual <- annual[, dp := log(D12) - log(Index)]

##  Log 3-month treasury bill (logRfree)     

annual <- annual[, logRfree := log(Rfree + 1)]


##  Term spread (tms)        

annual <- annual[, tms := lty - tbl]

##  Dividend yield vector (vec_dy)    


vec_dy <- c(NA, annual[2:nrow(annual), log(D12)] - annual[1:(nrow(annual)-1), log(Index)])


##  Dividend vector (dy) 

annual <- annual[, dy := vec_dy]

annual <- annual[, logret   := c(NA,diff(log(Index)))]

vec_logretdiv <- c(NA, annual[2:nrow(annual), log(IndexDiv)] - annual[1:(nrow(annual)-1), log(Index)])

vec_logretdiv <- c(NA, log(annual[2:nrow(annual), IndexDiv]/annual[1:(nrow(annual)-1), Index]))

annual <- annual[, logretdiv := vec_logretdiv]


##  Earnings price ratio (ep)  

annual <- annual[, ep := log(E12) - log(Index)]

##  Dividend payout ratio (de)  

annual <- annual[, de := log(D12) - log(E12)] 


#Log equity premium (rp_div)   

annual <- annual[, rp_div := logretdiv - logRfree]

#define news 

annual <- annual[, news := news]

#define vix

annual <- annual[, vix := vix]





#   ____________________________________________________________________________
#   TIME SERIES                                                             ####
#   ____________________________________________________________________________

##  ............................................................................
##  Annual data time series object (ts_annual)                              ####
##  ............................................................................

ts_annual <- ts(annual, start=annual[1, yyyy], end=annual[nrow(annual), yyyy])

##  ............................................................................
##  Plot ts_annual                                                          ####
##  ............................................................................
##  log equity premium, dividend-price ratio, dividend yield

plot(ts_annual[, c("rp_div", "dp", "dy", "news", "vix")])





####################################################################################################################


#######################################################################################################################################
#   ____________________________________________________________________________
#   STATISTICS FUNCTION (2)                                                 ####
#   ____________________________________________________________________________
#   Workaround function to adjust chart scale for stock variance (svar_stat)
get_statistics2 <- function(ts_df, indep, dep, h=1, start=1872, end=2023, est_periods_OOS = 20) {
  
  ##  ............................................................................
  ##  In Sample                                                               ####
  ##  ............................................................................
  
  avg <- mean(window(ts_df, start, end)[, dep], na.rm=TRUE)
  
  IS_error_N <- (window(ts_df, start, end)[, dep] - avg)
  
  reg <- dyn$lm(eval(parse(text=dep)) ~ stats::lag(eval(parse(text=indep)), -1), data=window(ts_df, start, end))
  
  IS_error_A <- reg$residuals
  
  ## Calculate Confidence Intervals
  ci <- confint(reg, level = 0.90)
  
  ##  ............................................................................
  ##  Out Of Sample                                                           ####
  ##  ............................................................................
  
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
  
  T <- length(!is.na(ts_df[, dep]))
  
  OOS_R2  <- 1 - MSE_A/MSE_N
  OOS_aR2 <- 1 - (((1-OOS_R2)*(reg_OOS$df.residual))/(reg_OOS$df.residual - 1))
  dRMSE <- sqrt(MSE_N) - sqrt(MSE_A)
  
  IS  <- cumsum(IS_error_N[2:length(IS_error_N)]^2)-cumsum(IS_error_A^2)
  OOS <- cumsum(OOS_error_N^2)-cumsum(OOS_error_A^2)
  
  df  <- data.frame(x=seq.int(from=start + 1 + est_periods_OOS, to=end), 
                    IS=IS[(1 + est_periods_OOS):length(IS)], 
                    OOS=OOS) 
  
  df$IS <- df$IS - df$IS[1] 
  df  <- melt(df, id.var="x") 
  plotGG <- ggplot(df) + 
    geom_line(aes(x=x, y=value,color=variable)) + 
    geom_rect(
      data=data.frame(),
      aes(xmin=1929, xmax=1941, ymin=-0.3, ymax=0.3),
      fill='darkolivegreen', 
      alpha=0.3
    ) +
    geom_rect(
      data=data.frame(),
      aes(xmin=2020, xmax=2023,ymin=-0.3,ymax=0.3),
      fill='deepskyblue1',
      alpha=0.3
    ) + 
    geom_rect(
      data=data.frame(), 
      aes(xmin=1974, xmax=1975, ymin=-0.3, ymax=0.3), 
      fill='red',
      alpha=0.3
    ) +
    geom_rect(
      data=data.frame(),
      aes(xmin=1987, xmax=1988,ymin=-0.3,ymax=0.3),
      fill='blue',
      alpha=0.3
    ) +
    geom_rect(
      data=data.frame(),
      aes(xmin=2000, xmax=2002,ymin=-0.3, ymax=0.3),
      fill='green',
      alpha=0.3
    ) +
    geom_rect(
      data=data.frame(),
      aes(xmin=2007, xmax=2010,ymin=-0.3, ymax=0.3),
      fill='purple',
      alpha=0.3
    ) +
    annotate(geom = "text", x = 1973, y = 0.15, label = "Oil Shock", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 1976, y = 0.15, label = "1974", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 1986, y = 0.15, label = "Black Monday", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 1989, y = 0.15, label = "1987", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 1999, y = 0.15, label = "Dot-com Burst", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 2003, y = 0.15, label = "2000-2002", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 2006, y = 0.15, label = "Global Financial Crisis", color = "Black",
             angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 2011, y = 0.15, label = "2008-2010", color ="Black", angle = 90, size = 2.5
    ) +
    annotate(geom = "text", x = 1928, y= 0.15, label = "The Great Depression", color = "Black", angle = 90, size = 2.5
    ) +
    annotate (geom = "text", x =1942, y = 0.15, label = "1929-1941", color = "Black", angle = 90, size = 2.5
    ) +
    annotate (geom = 'text', x= 2019, y= 0.15, label = "Covid-19 Low", color = "Black", angle = 90, size = 2.5
    ) + 
    annotate(geom = "text", x = 2024, y= 0.15, label = "2020-2022", color = "Black", angle = 90, size = 2.5
    ) +
    scale_x_continuous('Year') +
    scale_y_continuous('Cumulative SSE Difference', limits=c(-1.0, 0.3)) 
  
  # Conduct Chow Test
  chow_test <- sctest(reg, type = "Chow", point = (start + end) / 2)
  
  return(list(IS_error_N = IS_error_N,
              IS_error_A = reg$residuals,
              OOS_error_N = OOS_error_N,
              OOS_error_A = OOS_error_A,
              IS_R2 = summary(reg)$r.squared, 
              IS_aR2 = summary(reg)$adj.r.squared, 
              OOS_R2  = OOS_R2,
              OOS_aR2 = OOS_aR2,
              dRMSE = dRMSE,
              plotGG = plotGG,
              conf_int = ci,
              chow_test = chow_test
  ))
}

#   ____________________________________________________________________________
#   CREATE PLOTS & STATISTICS                                               ####
#   ____________________________________________________________________________

##  ............................................................................
##  Dividend-price ratio (dp_stat)                                          ####
##  ............................................................................

dp_stat <- get_statistics2(ts_annual, "dp", "rp_div", start=1872)
dp_stat$plotGG + ggtitle("Dividend-Price Ratio (dp)") + theme(plot.title = element_text(hjust = 0.5))
print(dp_stat$conf_int)
print(dp_stat$chow_test)

##  ............................................................................
##  Dividend-yield (dy_stat)                                                ####
##  ............................................................................

dy_stat <- get_statistics2(ts_annual, "dy", "rp_div", start=1872)
dy_stat$plotGG + ggtitle("Dividend Yield (dy)") + theme(plot.title = element_text(hjust = 0.5))
print(dy_stat$conf_int)
print(dy_stat$chow_test)

##  ............................................................................
##  Earnings-price ratio (ep_stat)                                          ####
##  ............................................................................

ep_stat <- get_statistics2(ts_annual, "ep", "rp_div", start=1872)
ep_stat$plotGG + ggtitle("Earnings-Price Ratio (ep)") + theme(plot.title = element_text(hjust = 0.5))
print(ep_stat$conf_int)
print(ep_stat$chow_test)

##  ............................................................................
##  Dividend payout ratio (de_stat)                                         ####
##  ............................................................................

de_stat <- get_statistics2(ts_annual, "de", "rp_div", start=1872)
de_stat$plotGG + ggtitle("Dividend Payout Ratio (de)") + theme(plot.title = element_text(hjust = 0.5))
print(de_stat$conf_int)
print(de_stat$chow_test)

##  ............................................................................
##  Stock variance (svar_stat)                                              ####
##  ............................................................................

svar_stat <- get_statistics2(ts_annual, "svar", "rp_div", start=1885)
svar_stat$plotGG + ggtitle("Stock Variance (svar)") + theme(plot.title = element_text(hjust = 0.5))
print(svar_stat$conf_int)
print(svar_stat$chow_test)

##  ............................................................................
##  Treasury bill rate (tbl_stat)                                           ####
##  ............................................................................

tbl_stat <- get_statistics2(ts_annual, "tbl", "rp_div", start=1920)
tbl_stat$plotGG + ggtitle("Treasury Bill Rate (tbl)") + theme(plot.title = element_text(hjust = 0.5))
print(tbl_stat$conf_int)
print(tbl_stat$chow_test)

##  ............................................................................
##  Long-term yield (lty_stat)                                              ####
##  ............................................................................

lty_stat <- get_statistics2(ts_annual, "lty", "rp_div", start=1919)
lty_stat$plotGG + ggtitle("Long-Term Yield (lty)") + theme(plot.title = element_text(hjust = 0.5))
print(lty_stat$conf_int)
print(lty_stat$chow_test)

##  ............................................................................
##  Long-term return (ltr_stat)                                             ####
##  ............................................................................

ltr_stat <- get_statistics2(ts_annual, "ltr", "rp_div", start=1926)
ltr_stat$plotGG + ggtitle("Long-Term Return (ltr)") + theme(plot.title = element_text(hjust = 0.5))
print(ltr_stat$conf_int)
print(ltr_stat$chow_test)

##  ............................................................................
##  Term spread (tms_stat)                                                  ####
##  ............................................................................

tms_stat <- get_statistics2(ts_annual, "tms", "rp_div", start=1920)
tms_stat$plotGG + ggtitle("Term Spread (tms)") + theme(plot.title = element_text(hjust = 0.5))
print(tms_stat$conf_int)
print(tms_stat$chow_test)

##  ............................................................................
##  Default yield spread (dfy_stat)                                         ####
##  ............................................................................

dfy_stat <- get_statistics2(ts_annual, "dfy", "rp_div", start=1919)
dfy_stat$plotGG + ggtitle("Default Yield Spread (dfy)") + theme(plot.title = element_text(hjust = 0.5))
print(dfy_stat$conf_int)
print(dfy_stat$chow_test)

##  ............................................................................
##  Default return spread (dfr_stat)                                        ####
##  ............................................................................

dfr_stat <- get_statistics2(ts_annual, "dfr", "rp_div", start=1926)
dfr_stat$plotGG + ggtitle("Default Return Spread (dfr)") + theme(plot.title = element_text(hjust = 0.5))
print(dfr_stat$conf_int)
print(dfr_stat$chow_test)

##  ............................................................................
##  Inflation (infl_stat)                                                   ####
##  ............................................................................

infl_stat <- get_statistics2(ts_annual, "infl", "rp_div", start=1919)
infl_stat$plotGG + ggtitle("Inflation (infl)") + theme(plot.title = element_text(hjust = 0.5))
print(infl_stat$conf_int)
print(infl_stat$chow_test)

##  ............................................................................
##  Book to market (bm_stat)                                                ####
##  ............................................................................

bm_stat <- get_statistics2(ts_annual, "bm", "rp_div", start=1921)
bm_stat$plotGG + ggtitle("Book to Market (bm)") + theme(plot.title = element_text(hjust = 0.5))
print(bm_stat$conf_int)
print(bm_stat$chow_test)

##  ............................................................................
##  Investment-capital ratio (ik_stat)                                      ####
##  ............................................................................

ik_stat <- get_statistics2(ts_annual, "ik", "rp_div", start=1947)
ik_stat$plotGG + ggtitle("Investment-Capital Ratio (ik)") + theme(plot.title = element_text(hjust = 0.5))
print(ik_stat$conf_int)
print(ik_stat$chow_test)

##  ............................................................................
##  Net equity expansion (ntis_stat)                                        ####
##  ............................................................................

ntis_stat <- get_statistics2(ts_annual, "ntis", "rp_div", start=1927)
ntis_stat$plotGG + ggtitle("Net Equity Expansion (ntis)") + theme(plot.title = element_text(hjust = 0.5))
print(ntis_stat$conf_int)
print(ntis_stat$chow_test)

##  ............................................................................
##  Percent equity issuing (eqis_stat)                                      ####
##  ............................................................................

eqis_stat <- get_statistics2(ts_annual, "eqis", "rp_div", start=1927)
eqis_stat$plotGG + ggtitle("Percent Equity Issuing (eqis)")
print(eqis_stat$conf_int)
print(eqis_stat$chow_test)

##news 
news_stat <- get_statistics2(ts_annual, "news", "rp_div", start=1980)
news_stat$plotGG + ggtitle("Daily News Sentiment Index")
print(news_stat$conf_int)
print(news_stat$chow_test)

##vix
vix_stat <- get_statistics2(ts_annual, "vix", "rp_div", start=1990)
vix_stat$plotGG + ggtitle("VIX")
print(vix_stat$conf_int)
print(vix_stat$chow_test)

#   ____________________________________________________________________________
#   TABLE 1                                                                 ####
#   ____________________________________________________________________________

##  ............................................................................
##  Create tibble                                                           ####
##  ............................................................................

table1 <- tibble(
  
  Variable = c(
    "dfy",
    "infl",
    "svar",
    "d/e",
    "lty",
    "tms",
    "tbl",
    "dfr",
    "d/p",
    "d/y",
    "ltr",
    "e/p",
    "b/m",
    "i/k",
    "ntis",
    "eqis",
    "news",
    "vix"
  ),
  
  
  IS_error_N = c(
    mean(dfy_stat$IS_error_N),
    mean(infl_stat$IS_error_N),
    mean(svar_stat$IS_error_N),
    mean(de_stat$IS_error_N),
    mean(lty_stat$IS_error_N),
    mean(tms_stat$IS_error_N),
    mean(tbl_stat$IS_error_N),
    mean(dfr_stat$IS_error_N),
    mean(dp_stat$IS_error_N),
    mean(dy_stat$IS_error_N),
    mean(ltr_stat$IS_error_N),
    mean(ep_stat$IS_error_N),
    mean(bm_stat$IS_error_N),
    mean(ik_stat$IS_error_N),
    mean(ntis_stat$IS_error_N),
    mean(eqis_stat$IS_error_N),
    mean(news_stat$IS_error_N),
    mean(vix_stat$IS_error_N)
  ),
  
  IS_error_A = c(
    mean(dfy_stat$IS_error_A),
    mean(infl_stat$IS_error_A),
    mean(svar_stat$IS_error_A),
    mean(de_stat$IS_error_A),
    mean(lty_stat$IS_error_A),
    mean(tms_stat$IS_error_A),
    mean(tbl_stat$IS_error_A),
    mean(dfr_stat$IS_error_A),
    mean(dp_stat$IS_error_A),
    mean(dy_stat$IS_error_A),
    mean(ltr_stat$IS_error_A),
    mean(ep_stat$IS_error_A),
    mean(bm_stat$IS_error_A),
    mean(ik_stat$IS_error_A),
    mean(ntis_stat$IS_error_A),
    mean(eqis_stat$IS_error_A),
    mean(news_stat$IS_error_A),
    mean(vix_stat$IS_error_A)
  ),
  OOS_error_N = c(
    mean(dfy_stat$OOS_error_N),
    mean(infl_stat$OOS_error_N),
    mean(svar_stat$OOS_error_N),
    mean(de_stat$OOS_error_N),
    mean(lty_stat$OOS_error_N),
    mean(tms_stat$OOS_error_N),
    mean(tbl_stat$OOS_error_N),
    mean(dfr_stat$OOS_error_N),
    mean(dp_stat$OOS_error_N),
    mean(dy_stat$OOS_error_N),
    mean(ltr_stat$OOS_error_N),
    mean(ep_stat$OOS_error_N),
    mean(bm_stat$OOS_error_N),
    mean(ik_stat$OOS_error_N),
    mean(ntis_stat$OOS_error_N),
    mean(eqis_stat$OOS_error_N),
    mean(news_stat$OOS_error_N),
    mean(vix_stat$OOS_error_N)
  ),
  OOS_error_A = c(
    mean(dfy_stat$OOS_error_A),
    mean(infl_stat$OOS_error_A),
    mean(svar_stat$OOS_error_A),
    mean(de_stat$OOS_error_A),
    mean(lty_stat$OOS_error_A),
    mean(tms_stat$OOS_error_A),
    mean(tbl_stat$OOS_error_A),
    mean(dfr_stat$OOS_error_A),
    mean(dp_stat$OOS_error_A),
    mean(dy_stat$OOS_error_A),
    mean(ltr_stat$OOS_error_A),
    mean(ep_stat$OOS_error_A),
    mean(bm_stat$OOS_error_A),
    mean(ik_stat$OOS_error_A),
    mean(ntis_stat$OOS_error_A),
    mean(eqis_stat$OOS_error_A),
    mean(news_stat$OOS_error_A),
    mean(vix_stat$OOS_error_A)
  ),
  IS_R2 = c(
    dfy_stat$IS_R2,
    infl_stat$IS_R2,
    svar_stat$IS_R2,
    de_stat$IS_R2,
    lty_stat$IS_R2,
    tms_stat$IS_R2,
    tbl_stat$IS_R2,
    dfr_stat$IS_R2,
    dp_stat$IS_R2,
    dy_stat$IS_R2,
    ltr_stat$IS_R2,
    ep_stat$IS_R2,
    bm_stat$IS_R2,
    ik_stat$IS_R2,
    ntis_stat$IS_R2,
    eqis_stat$IS_R2,
    news_stat$IS_R2,
    vix_stat$IS_R2
  ),
  
  IS_aR2 = c(
    format(round(dfy_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(infl_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(svar_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(de_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(lty_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(tms_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(tbl_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(dfr_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(dp_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(dy_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(ltr_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(ep_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(bm_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(ik_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(ntis_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(eqis_stat$IS_aR2*100,digits=2),nsmall=2), 
    format(round(news_stat$IS_aR2*100,digits=2),nsmall=2),
    format(round(vix_stat$IS_aR2*100,digits=2),nsmall=2)
  ),
  
  OOS_R2 = c(
    dfy_stat$OOS_R2,
    infl_stat$OOS_R2,
    svar_stat$OOS_R2,
    de_stat$OOS_R2,
    lty_stat$OOS_R2,
    tms_stat$OOS_R2,
    tbl_stat$OOS_R2,
    dfr_stat$OOS_R2,
    dp_stat$OOS_R2,
    dy_stat$OOS_R2,
    ltr_stat$OOS_R2,
    ep_stat$OOS_R2,
    bm_stat$OOS_R2,
    ik_stat$OOS_R2,
    ntis_stat$OOS_R2,
    eqis_stat$OOS_R2,
    news_stat$OOS_R2,
    vix_stat$OOS_R2
  ),
  
  OOS_aR2 = c(
    format(round(dfy_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(infl_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(svar_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(de_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(lty_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(tms_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(tbl_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(dfr_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(dp_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(dy_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(ltr_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(ep_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(bm_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(ik_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(ntis_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(eqis_stat$OOS_aR2*100,digits=2),nsmall=2), 
    format(round(news_stat$OOS_aR2*100,digits=2),nsmall=2),
    format(round(vix_stat$OOS_aR2*100,digits=2),nsmall=2)
  ),
  
  dRMSE = c(
    format(round(dfy_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(infl_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(svar_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(de_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(lty_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(tms_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(tbl_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(dfr_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(dp_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(dy_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(ltr_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(ep_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(bm_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(ik_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(ntis_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(eqis_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(news_stat$dRMSE*100,digits=2),nsmall=2),
    format(round(vix_stat$dRMSE*100,digits=2),nsmall=2)
  )
)

# Convert to numeric format, handling the already formatted string columns correctly

########

###DECIMAL FORMATTING DONE HERE

##############
table1 <- table1 %>%
  mutate(
    IS_error_N = as.numeric(IS_error_N),
    IS_error_A = as.numeric(IS_error_A),
    OOS_error_N = as.numeric(OOS_error_N),
    OOS_error_A = as.numeric(OOS_error_A),
    IS_R2 = as.numeric(IS_R2),
    IS_aR2 = as.numeric(IS_aR2),
    OOS_R2 = as.numeric(OOS_R2),
    OOS_aR2 = as.numeric(OOS_aR2),
    dRMSE = as.numeric(dRMSE)
  ) %>%
  mutate(across(c(IS_error_N, IS_error_A, OOS_error_N, OOS_error_A, IS_R2, IS_aR2, OOS_R2, OOS_aR2, dRMSE), round, digits = 4))
####################




##  ............................................................................
##  Create table comparing results to G&W                                   ####
##  ............................................................................

# Create table comparing results to G&W
table1_sum <- tibble(
  Variable = table1$Variable, 
  IS_aR2 = c(dfy_stat$IS_aR2, infl_stat$IS_aR2, svar_stat$IS_aR2, de_stat$IS_aR2, lty_stat$IS_aR2, 
             tms_stat$IS_aR2, tbl_stat$IS_aR2, dfr_stat$IS_aR2, dp_stat$IS_aR2, dy_stat$IS_aR2, 
             ltr_stat$IS_aR2, ep_stat$IS_aR2, bm_stat$IS_aR2, ik_stat$IS_aR2, ntis_stat$IS_aR2, 
             eqis_stat$IS_aR2, news_stat$IS_aR2, vix_stat$IS_aR2),
  OOS_aR2 = c(dfy_stat$OOS_aR2, infl_stat$OOS_aR2, svar_stat$OOS_aR2, de_stat$OOS_aR2, lty_stat$OOS_aR2, 
              tms_stat$OOS_aR2, tbl_stat$OOS_aR2, dfr_stat$OOS_aR2, dp_stat$OOS_aR2, dy_stat$OOS_aR2, 
              ltr_stat$OOS_aR2, ep_stat$OOS_aR2, bm_stat$OOS_aR2, ik_stat$OOS_aR2, ntis_stat$OOS_aR2, 
              eqis_stat$OOS_aR2, news_stat$OOS_aR2, vix_stat$OOS_aR2),
  dRMSE = c(dfy_stat$dRMSE, infl_stat$dRMSE, svar_stat$dRMSE, de_stat$dRMSE, lty_stat$dRMSE, 
            tms_stat$dRMSE, tbl_stat$dRMSE, dfr_stat$dRMSE, dp_stat$dRMSE, dy_stat$dRMSE, 
            ltr_stat$dRMSE, ep_stat$dRMSE, bm_stat$dRMSE, ik_stat$dRMSE, ntis_stat$dRMSE, 
            eqis_stat$dRMSE, news_stat$dRMSE, vix_stat$dRMSE)
)


###################################VIEWER CODE SNIPPET####################################

# Define custom color functions
customGreen0 <- "#DeF7E9"
customGreen <- "#71CA97"
customRed <- "#ff7f7f"


# Define custom color functions based on intuitive thresholds
color_tile_good <- function(value) {
  sapply(value, function(x) style(color = if (x > 0) customGreen else customRed))
}

color_tile_bad <- function(value) {
  sapply(value, function(x) style(color = if (x < 0) customGreen else customRed))
}

# Define colors for IS and OOS errors (positive errors are generally not good)
color_tile_error <- function(value) {
  sapply(value, function(x) style(color = if (x <= 0) customGreen else customRed))
}

# Create a formattable object with better styling
formatted_table <- formattable(table1, align = c("l", rep("r", NCOL(table1) - 1)),
                               list(
                                 `Variable` = formatter("span", style = ~ style(color = "grey", font.weight = "bold")),
                                 IS_error_N = formatter("span", style = ~ color_tile_error(IS_error_N)),
                                 IS_error_A = formatter("span", style = ~ color_tile_error(IS_error_A)),
                                 OOS_error_N = formatter("span", style = ~ color_tile_error(OOS_error_N)),
                                 OOS_error_A = formatter("span", style = ~ color_tile_error(OOS_error_A)),
                                 IS_R2 = formatter("span", style = ~ color_tile_good(IS_R2)),
                                 IS_aR2 = formatter("span", style = ~ color_tile_good(IS_aR2)),
                                 OOS_R2 = formatter("span", style = ~ color_tile_good(OOS_R2)),
                                 OOS_aR2 = formatter("span", style = ~ color_tile_good(OOS_aR2)),
                                 dRMSE = formatter("span", style = ~ color_tile_good(dRMSE))
                               )
)

# Convert formattable to an HTML widget
formatted_table_widget <- as.htmlwidget(formatted_table)

# Add a title to the table
formatted_table_with_title <- htmltools::browsable(htmltools::tagList(
  htmltools::tags$h3("Predictive Metrics Table"),
  formatted_table_widget
))

# Define the legend HTML
legend_html <- '
<div style="position: absolute; bottom: 10px; right: 10px;">
  <h4>Legend</h4>
  <div><span style="background-color: #71CA97; padding: 2px 10px;">Positive</span></div>
  <div><span style="background-color: #ff7f7f; padding: 2px 10px;">Negative</span></div>
  <div><span style="background-color: #DeF7E9; padding: 2px 10px;">Neutral</span></div>
</div>
'

# Convert formattable to an HTML widget with legend
formatted_table_with_legend <- htmltools::browsable(htmltools::tagList(
  htmltools::tags$h3("Predictive Metrics Table"),
  formatted_table_widget,
  htmltools::HTML(legend_html)
))

# Print the formatted table with the title and legend
print(formatted_table_with_legend)

# Use the datatable object for better interactivity in the viewer
datatable_object <- datatable(table1, options = list(paging = FALSE, filter = 'top'))

# Print or view the datatable
datatable_object
