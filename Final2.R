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
library(knitr)
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

annual <- annual[, diffdp := diffdp]

##  Log 3-month treasury bill (logRfree)     

annual <- annual[, logRfree := log(Rfree + 1)]

#define differenced lty

annual <- annual[, difflty := difflty]

#define differenced book-to-market ratio

annual <- annual[, diffbm := diffbm]

#define differenced tbl 

annual <- annual[, difftbl := difftbl]

##  Term spread (tms)        

annual <- annual[, tms := lty - tbl]

# Use diffdy directly instead of calculating dy
annual <- annual[, dy := as.numeric(diffdy)]     # Ensure diffdy is numeric

# Calculate log returns (logret)
annual <- annual[, logret := c(NA, diff(log(Index)))]

# Calculate log returns with dividends (logretdiv)
vec_logretdiv <- c(NA, log(annual[2:nrow(annual), IndexDiv] / annual[1:(nrow(annual)-1), Index]))
annual <- annual[, logretdiv := vec_logretdiv]

##  Earnings price ratio (ep)  

annual <- annual[, ep := log(E12) - log(Index)]

##  Dividend payout ratio (de)  

annual <- annual[, diffde := diffde] 


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

plot(ts_annual[, c("rp_div", "diffdp", "diffdy", "news", "vix", "difftbl")])





####################################################################################################################


#######################################################################################################################################
#   ____________________________________________________________________________
#   STATISTICS FUNCTION (2)                                                 ####
#   ____________________________________________________________________________

# Define a function to perform Bai-Perron test
perform_bai_perron_test <- function(ts_df, indep, dep, start, end) {
  # Subset the data to the relevant period
  data_subset <- window(ts_df, start, end)
  
  # Build the formula for regression
  formula <- as.formula(paste(dep, "~", indep))
  
  # Perform the Bai-Perron test
  bp_result <- breakpoints(formula, data = data_subset)
  
  return(bp_result)
}





get_statistics2 <- function(ts_df, indep, dep, h=1, start=1872, end=2023, est_periods_OOS = 20)
  
{
  
  ##  ............................................................................
  ##  In Sample                                                               ####
  ##  ............................................................................
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### IS Historical mean model                                                ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### na.rm=TRUE - Excludes missing values from analyses                      
  
  avg <- mean(window(ts_df, start, end)[, dep], na.rm=TRUE)
  
  
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### IS historical mean error                                                ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  IS_error_N <- (window(ts_df, start, end)[, dep] - avg)
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### IS OLS Regression - dyn package used to regress using lagged variable   ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  reg <- dyn$lm(eval(parse(text=dep)) ~ stats::lag(eval(parse(text=indep)), -1), data=window(ts_df, start, end))
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### IS OLS error                                                            ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  IS_error_A <- reg$residuals
  
  ## Calculate Confidence Intervals
  ci <- confint(reg, level = 0.90)
  
  
  ##  ............................................................................
  ##  Out Of Sample                                                           ####
  ##  ............................................................................
  ##  Uses recursive estimation window                                        ####
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Create OOS historical mean error as "numeric" data type                 ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  OOS_error_N <- numeric(end - start - est_periods_OOS)
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Create OOS OLS error as "numeric" data type                             ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  OOS_error_A <- numeric(end - start - est_periods_OOS)
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Only use information that is available up to the time of forecast       ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  j <- 0
  
  for (i in (start + est_periods_OOS):(end-1)) 
    
  {
    
    j <- j + 1
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### Actual ERP                                                              ####
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    
    actual_ERP <- as.numeric(window(ts_df, i+1, i+1)[, dep])
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### OOS Historical mean model                                               ####
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    
    OOS_error_N[j] <- actual_ERP - mean(window(ts_df, start, i)[, dep], na.rm=TRUE)
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### OOS OLS model                                                           ####
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    
    reg_OOS <- dyn$lm(eval(parse(text=dep)) ~ stats::lag(eval(parse(text=indep)), -1), data=window(ts_df, start, i))
    
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    ### OOS OLS error                                                           ####
    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
    
    df <- data.frame(x=as.numeric(window(ts_df, i, i)[, indep]))
    names(df) <- indep
    pred_ERP   <- predict.lm(reg_OOS, newdata=df)
    OOS_error_A[j] <-  pred_ERP - actual_ERP
    
  }
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Historical mean model                                                   ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### vector of rolling OOS errors                                            
  
  MSE_N <- mean(OOS_error_N^2)
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### OLS model                                                               ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### vector of rolling OOS errors                                            
  
  MSE_A <- mean(OOS_error_A^2) 
  
  T <- length(!is.na(ts_df[, dep]))
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### OOS R-squared                                                           ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  OOS_R2  <- 1 - MSE_A/MSE_N
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### OOS Adjusted R-squared (OOS_aR2)                                        ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  
  OOS_aR2 <- 1 - (((1-OOS_R2)*(reg_OOS$df.residual))/(reg_OOS$df.residual - 1))
  
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Delta root mean squared error (dRMSE)                                   ####
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### difference between the historical mean model RMSE and the OLS model RMSE
  
  dRMSE <- sqrt(MSE_N) - sqrt(MSE_A)
  
  #   ____________________________________________________________________________
  #   CREATE PLOT                                                             ####
  #   ____________________________________________________________________________
  
  ##  ............................................................................
  ##  IS Cummulative SSE Difference                                           ####
  ##  ............................................................................
  
  IS  <- cumsum(IS_error_N[2:length(IS_error_N)]^2)-cumsum(IS_error_A^2)
  
  ##  ............................................................................
  ##  OOS Cummulative SSE Difference                                          ####
  ##  ............................................................................
  
  OOS <- cumsum(OOS_error_N^2)-cumsum(OOS_error_A^2)
  
  df  <- data.frame(x=seq.int(from=start + 1 + est_periods_OOS, to=end), 
                    IS=IS[(1 + est_periods_OOS):length(IS)], 
                    OOS=OOS) # One observation lost due to the lag
  
  ##  ............................................................................
  ##  Shift IS errors vertically                                              ####
  ##  ............................................................................
  ##  Sets IS line to begin at zero on the date of first OOS prediction       
  
  df$IS <- df$IS - df$IS[1] 
  df  <- melt(df, id.var="x") 
  plotGG <- ggplot(df) + 
    geom_line(aes(x=x, y=value,color=variable)) + 
    
    
    
    #HERE COME THE PLOTS ############################################################################################
  
  
  ##  ............................................................................
  ##  Highlight oil shock of 1974                                             ####
  ##  ............................................................................
  
  
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
      data=data.frame(),# Needed by ggplot2 for transparency
      aes(xmin=1987, xmax=1988,ymin=-0.3,ymax=0.3),
      fill='blue',
      alpha=0.3
    ) +
    
    geom_rect(
      data=data.frame(),# Needed by ggplot2 for transparency
      aes(xmin=2000, xmax=2002,ymin=-0.3,ymax=0.3),
      fill='green',
      alpha=0.3
    ) +
    
    geom_rect(
      data=data.frame(),# Needed by ggplot2 for transparency
      aes(xmin=2007, xmax=2010,ymin=-0.3,ymax=0.3),
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
    
    
    ##  ............................................................................
    ##  Label x-axis and y-axis                                                 ####
  ##  ............................................................................
  
  scale_x_continuous('Year') +
    scale_y_continuous('Cumulative SSE Difference', limits=c(-1.0, 0.3)) 
  
  # Conduct Chow Test
  chow_test <- sctest(reg, type = "Chow", point = (start + end) / 2)
  # Perform Bai-Perron test for structural breaks
  bp_result <- perform_bai_perron_test(ts_df, indep, dep, start, end)
  
  ##  ............................................................................
  ##  Return function values                                                  ####
  ##  ............................................................................
  
  return(list(IS_error_N = IS_error_N,
              IS_error_A = reg$residuals,
              OOS_error_N = OOS_error_N,
              OOS_error_A = OOS_error_A,
              IS_R2 = summary(reg)$r.squared, 
              IS_aR2 = summary(reg)$adj.r.squared, 
              OOS_R2  = OOS_R2,
              OOS_aR2 = OOS_aR2,
              dRMSE = dRMSE,
              conf_int = ci,
              chow_test = chow_test,
              bai_perron_test = bp_result,
              plotGG = plotGG
  ))
}

#   ____________________________________________________________________________
#   CREATE PLOTS & STATISTICS                                               ####
#   ____________________________________________________________________________

##  ............................................................................
##  Dividend-price ratio (dp_stat)                                          ####
##  ............................................................................

dp_stat <- get_statistics2(ts_annual, "diffdp", "rp_div", start=1873)
dp_stat$plotGG + ggtitle("Dividend-Price Ratio (differenced)") + theme(plot.title = element_text(hjust = 0.5))
print(dp_stat$conf_int)
print(dp_stat$chow_test)
print(summary(dp_stat$bai_perron_test))


##  ............................................................................
##  Dividend-yield (dy_stat)                                                ####
##  ............................................................................

dy_stat <- get_statistics2(ts_annual, "diffdy", "rp_div", start=1874)
dy_stat$plotGG + ggtitle("Dividend Yield (differenced)") + theme(plot.title = element_text(hjust = 0.5))
print(dy_stat$conf_int)
print(dy_stat$chow_test)
print(summary(dy_stat$bai_perron_test))

##  ............................................................................
##  Earnings-price ratio (ep_stat)                                          ####
##  ............................................................................

ep_stat <- get_statistics2(ts_annual, "ep", "rp_div", start=1872)
ep_stat$plotGG + ggtitle("Earnings-Price Ratio (ep)") + theme(plot.title = element_text(hjust = 0.5))
print(ep_stat$conf_int)
print(ep_stat$chow_test)
print(summary(ep_stat$bai_perron_test))

##  ............................................................................
##  Dividend payout ratio (de_stat)                                         ####
##  ............................................................................

de_stat <- get_statistics2(ts_annual, "diffde", "rp_div", start=1873)
de_stat$plotGG + ggtitle("Dividend Payout Ratio (differenced)") + theme(plot.title = element_text(hjust = 0.5))
print(de_stat$conf_int)
print(de_stat$chow_test)
print(summary(de_stat$bai_perron_test))


##  ............................................................................
##  Stock variance (svar_stat)                                              ####
##  ............................................................................

svar_stat <- get_statistics2(ts_annual, "svar", "rp_div", start=1885)
svar_stat$plotGG + ggtitle("Stock Variance (svar)") + theme(plot.title = element_text(hjust = 0.5))
print(svar_stat$conf_int)
print(svar_stat$chow_test)
print(summary(svar_stat$bai_perron_test))

##  ............................................................................
##  Treasury bill rate (tbl_stat)                                           ####
##  ............................................................................

tbl_stat <- get_statistics2(ts_annual, "difftbl", "rp_div", start=1921)
tbl_stat$plotGG + ggtitle("Treasury Bill Rate (differenced)") + theme(plot.title = element_text(hjust = 0.5))
print(tbl_stat$conf_int)
print(tbl_stat$chow_test)
print(summary(tbl_stat$bai_perron_test))


..............................................
##  Long-term yield (lty_stat)                                              ####
##  ............................................................................

lty_stat <- get_statistics2(ts_annual, "difflty", "rp_div", start=1920)
lty_stat$plotGG + ggtitle("Long-Term Yield (differenced)") + theme(plot.title = element_text(hjust = 0.5))
print(lty_stat$conf_int)
print(lty_stat$chow_test)
print(summary(lty_stat$bai_perron_test))


##  ............................................................................
##  Long-term return (ltr_stat)                                             ####
##  ............................................................................

ltr_stat <- get_statistics2(ts_annual, "ltr", "rp_div", start=1926)
ltr_stat$plotGG + ggtitle("Long-Term Return (ltr)") + theme(plot.title = element_text(hjust = 0.5))
print(ltr_stat$conf_int)
print(ltr_stat$chow_test)
print(summary(ltr_stat$bai_perron_test))

##  ............................................................................
##  Term spread (tms_stat)                                                  ####
##  ............................................................................

tms_stat <- get_statistics2(ts_annual, "tms", "rp_div", start=1920)
tms_stat$plotGG + ggtitle("Term Spread (tms)") + theme(plot.title = element_text(hjust = 0.5))
print(tms_stat$conf_int)
print(tms_stat$chow_test)
print(summary(tms_stat$bai_perron_test))

##  ............................................................................
##  Default yield spread (dfy_stat)                                         ####
##  ............................................................................

dfy_stat <- get_statistics2(ts_annual, "dfy", "rp_div", start=1919)
dfy_stat$plotGG + ggtitle("Default Yield Spread (dfy)") + theme(plot.title = element_text(hjust = 0.5))
print(dfy_stat$conf_int)
print(dfy_stat$chow_test)
print(summary(dfy_stat$bai_perron_test))

##  ............................................................................
##  Default return spread (dfr_stat)                                        ####
##  ............................................................................

dfr_stat <- get_statistics2(ts_annual, "dfr", "rp_div", start=1926)
dfr_stat$plotGG + ggtitle("Default Return Spread (dfr)") + theme(plot.title = element_text(hjust = 0.5))
print(dfr_stat$conf_int)
print(dfr_stat$chow_test)
print(summary(dfr_stat$bai_perron_test))

##  ............................................................................
##  Inflation (infl_stat)                                                   ####
##  ............................................................................

infl_stat <- get_statistics2(ts_annual, "infl", "rp_div", start=1919)
infl_stat$plotGG + ggtitle("Inflation (infl)") + theme(plot.title = element_text(hjust = 0.5))
print(infl_stat$conf_int)
print(infl_stat$chow_test)
print(summary(infl_stat$bai_perron_test))

##  ............................................................................
##  Book to market (bm_stat)                                                ####
##  ............................................................................

bm_stat <- get_statistics2(ts_annual, "diffbm", "rp_div", start=1922)
bm_stat$plotGG + ggtitle("Book to Market (differenced)") + theme(plot.title = element_text(hjust = 0.5))
print(bm_stat$conf_int)
print(bm_stat$chow_test)
print(summary(bm_stat$bai_perron_test))


##  ............................................................................
##  Investment-capital ratio (ik_stat)                                      ####
##  ............................................................................

ik_stat <- get_statistics2(ts_annual, "ik", "rp_div", start=1947)
ik_stat$plotGG + ggtitle("Investment-Capital Ratio (ik)") + theme(plot.title = element_text(hjust = 0.5))
print(ik_stat$conf_int)
print(ik_stat$chow_test)
print(summary(ik_stat$bai_perron_test))

##  ............................................................................
##  Net equity expansion (ntis_stat)                                        ####
##  ............................................................................

ntis_stat <- get_statistics2(ts_annual, "ntis", "rp_div", start=1927)
ntis_stat$plotGG + ggtitle("Net Equity Expansion (ntis)") + theme(plot.title = element_text(hjust = 0.5))
print(ntis_stat$conf_int)
print(ntis_stat$chow_test)
print(summary(ntis_stat$bai_perron_test))


##  ............................................................................
##  Percent equity issuing (eqis_stat)                                      ####
##  ............................................................................

eqis_stat <- get_statistics2(ts_annual, "eqis", "rp_div", start=1927)
eqis_stat$plotGG + ggtitle("Percent Equity Issuing (eqis)") 
print(eqis_stat$conf_int)
print(eqis_stat$chow_test)
print(summary(eqis_stat$bai_perron_test))


##news 
news_stat <- get_statistics2(ts_annual, "news", "rp_div", start=1980)
news_stat$plotGG + ggtitle("Daily News Sentiment Index") 
print(news_stat$conf_int)
print(news_stat$chow_test)
print(summary(news_stat$bai_perron_test))

##vix

vix_stat <- get_statistics2(ts_annual, "vix", "rp_div", start=1990)
vix_stat$plotGG + ggtitle("VIX") 
print(vix_stat$conf_int)
print(vix_stat$chow_test)
print(summary(vix_stat$bai_perron_test))



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
    "d/e (differenced)",
    "lty (differenced)",
    "tms",
    "tbl (differenced)",
    "dfr",
    "d/p(differenced)",
    "d/y (differenced)",
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


print(table1)

##  ............................................................................
##  Create table comparing results to G&W                                   ####
##  ............................................................................

print(table1_2005_g_w)
