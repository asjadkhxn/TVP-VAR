#FRAM G11

library(readxl)
library(zoo)
library(moments)
library(tseries)
library(urca)
library(openxlsx)
library(utils)  
library(ConnectednessApproach)
library(grDevices)

data_path <- "~/Downloads/Final_Data.xlsx"
data <- read_excel(data_path)
data$Date <- as.Date(data$Date)

# Summary statistics
calculate_summary_stats <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  var_val <- var(x, na.rm = TRUE)
  skew <- skewness(x, na.rm = TRUE)
  kurt <- kurtosis(x, na.rm = TRUE)
  jb_test <- jarque.bera.test(x)
  ers_test <- ur.ers(x, type = "P-test", model = "constant")
  q_test <- Box.test(x, lag = 20, type = "Ljung-Box")
  q_squared_test <- Box.test(x^2, lag = 20, type = "Ljung-Box")
  
  return(c(
    Mean = round(mean_val, 3),
    Variance = round(var_val, 3),
    Skewness = round(skew, 3),
    Kurtosis = round(kurt, 3),
    JB = round(jb_test$statistic, 3),
    ERS = round(ers_test@teststat[1], 3),
    "Q(20)" = round(q_test$statistic, 3),
    "Q²(20)" = round(q_squared_test$statistic, 3)
  ))
}

# Statistics for each commodity
commodities <- names(data)[-1]  # exclude Date column
stats_matrix <- sapply(data[commodities], calculate_summary_stats)

# Correlation matrix
cor_matrix <- cor(data[commodities], use = "complete.obs")


wb <- createWorkbook()
addWorksheet(wb, "Summary Statistics")
addWorksheet(wb, "Correlation Matrix")
writeData(wb, "Summary Statistics", t(stats_matrix), rowNames = TRUE)
writeData(wb, "Correlation Matrix", cor_matrix, rowNames = TRUE)

headerStyle <- createStyle(
  textDecoration = "bold",
  halign = "center",
  border = "bottom",
  borderColour = "black"
)

numberStyle <- createStyle(
  halign = "right",
  numFmt = "0.000"
)
addStyle(wb, "Summary Statistics", headerStyle, rows = 1, cols = 1:(ncol(stats_matrix)+1))
addStyle(wb, "Correlation Matrix", headerStyle, rows = 1, cols = 1:(ncol(cor_matrix)+1))
saveWorkbook(wb, "Table1_Summary_Statistics.xlsx", overwrite = TRUE)


print("Tables have been saved to Table1_Summary_Statistics.xlsx")





data_path <- "~/Downloads/Final_Data.xlsx"
data <- read_excel(data_path)


data$Date <- as.Date(data$Date)
data_zoo <- zoo(data[,-1], order.by=data$Date)

# COVID period
covid_start_date <- as.Date("2020-03-01")

# Split into pre-COVID and COVID periods
data_precovid <- window(data_zoo, end=covid_start_date-1)
data_covid <- window(data_zoo, start=covid_start_date)

# TVP-VAR connectedness for pre-COVID period
tvp_var_precovid <- ConnectednessApproach(
  x = data_precovid,
  model = "TVP-VAR",
  connectedness = "Time",
  nlag = 1,
  nfore = 20,
  VAR_config = list(
    TVPVAR = list(
      kappa1 = 0.99,
      kappa2 = 0.99,
      prior = "BayesPrior"
    )
  )
)

# TVP-VAR connectedness for COVID period
tvp_var_covid <- ConnectednessApproach(
  x = data_covid,
  model = "TVP-VAR",
  connectedness = "Time",
  nlag = 1,
  nfore = 20,
  VAR_config = list(
    TVPVAR = list(
      kappa1 = 0.99,
      kappa2 = 0.99,
      prior = "BayesPrior"
    )
  )
)


format_connectedness_table <- function(pre, covid) {
  
  n <- ncol(pre$CT[,,1])  # Using first time slice of CT array
  result <- matrix(NA, nrow = n + 2, ncol = n + 1)
  pre_means <- apply(pre$CT, c(1,2), mean, na.rm = TRUE)
  covid_means <- apply(covid$CT, c(1,2), mean, na.rm = TRUE)
  for(i in 1:n) {
    for(j in 1:n) {
      result[i,j] <- sprintf("%.2f (%.2f)", 
                             pre_means[i,j], 
                             covid_means[i,j])
    }
  }
  
  # FROM others
  pre_from <- rowSums(pre_means) - diag(pre_means)
  covid_from <- rowSums(covid_means) - diag(covid_means)
  
  # Add FROM others column
  for(i in 1:n) {
    result[i, n+1] <- sprintf("%.2f (%.2f)",
                              pre_from[i],
                              covid_from[i])
  }
  
  # Calculate and add TO others row
  pre_to <- colSums(pre_means) - diag(pre_means)
  covid_to <- colSums(covid_means) - diag(covid_means)
  for(j in 1:n) {
    result[n+1, j] <- sprintf("%.2f (%.2f)",
                              pre_to[j],
                              covid_to[j])
  }
  
  # Calculate and add NET row
  pre_net <- pre_to - pre_from
  covid_net <- covid_to - covid_from
  for(j in 1:n) {
    result[n+2, j] <- sprintf("%.2f (%.2f)",
                              pre_net[j],
                              covid_net[j])
  }
  
  # Add TCI (Total Connectedness Index)
  pre_tci <- mean(pre$TCI, na.rm = TRUE)
  covid_tci <- mean(covid$TCI, na.rm = TRUE)
  result[n+1, n+1] <- sprintf("%.2f (%.2f)",
                              pre_tci,
                              covid_tci)
  
  return(result)
}


result_table <- format_connectedness_table(tvp_var_precovid, tvp_var_covid)
colnames(result_table) <- c(colnames(data_zoo), "FROM others")
rownames(result_table) <- c(colnames(data_zoo), "TO others", "NET")


wb <- createWorkbook()
addWorksheet(wb, "Dynamic Connectedness")
writeData(wb, "Dynamic Connectedness", result_table, 
          rowNames = TRUE, colNames = TRUE)

writeData(wb, "Dynamic Connectedness", 
          "Results are based on a TVP-VAR model with a lag length of order one (BIC), κ1 = 0.99, κ2 = 0.99, and a 20-step-ahead generalized forecast error variance decomposition. Values in parentheses represent connectedness measures during the COVID-19 pandemic while others stand for the connectedness measures prior to the COVID-19 period",
          startRow = nrow(result_table) + 2, startCol = 1)
saveWorkbook(wb, "Table2_Dynamic_Connectedness.xlsx", overwrite = TRUE)




# We already have tvp_var_precovid and tvp_var_covid from previous code
# Pairwise connectedness

calculate_pairwise_connectedness <- function(tvp_var) {
  ct_array <- tvp_var$CT
  ct_mean <- apply(ct_array, c(1,2), mean, na.rm = TRUE)
  return(ct_mean * 100) 
}

pairwise_precovid <- calculate_pairwise_connectedness(tvp_var_precovid)
pairwise_covid <- calculate_pairwise_connectedness(tvp_var_covid)

wb <- createWorkbook()
addWorksheet(wb, "Pairwise Connectedness")
result_matrix <- matrix(NA, 
                        nrow = ncol(data_zoo), 
                        ncol = ncol(data_zoo))

for(i in 1:ncol(data_zoo)) {
  for(j in 1:ncol(data_zoo)) {
    result_matrix[i,j] <- sprintf("%.2f (%.2f)", 
                                  pairwise_precovid[i,j], 
                                  pairwise_covid[i,j])
  }
}


colnames(result_matrix) <- colnames(data_zoo)
rownames(result_matrix) <- colnames(data_zoo)
writeData(wb, "Pairwise Connectedness", result_matrix, 
          rowNames = TRUE, colNames = TRUE)

headerStyle <- createStyle(
  textDecoration = "bold",
  halign = "center",
  border = "bottom"
)

addStyle(wb, "Pairwise Connectedness", headerStyle, 
         rows = 1, cols = 1:(ncol(result_matrix)+1))
addStyle(wb, "Pairwise Connectedness", headerStyle, 
         rows = 1:(nrow(result_matrix)+1), cols = 1)


writeData(wb, "Pairwise Connectedness", 
          "Results are based on a TVP-VAR model with a lag length of order one (BIC), κ1 = 0.99, κ2 = 0.99, and a 20-step-ahead generalized forecast error variance decomposition. Values in parentheses represent connectedness measures during the COVID-19 pandemic while others stand for the connectedness measures prior to the COVID-19 period",
          startRow = nrow(result_matrix) + 2, startCol = 1)
saveWorkbook(wb, "Table3_Pairwise_Connectedness.xlsx", overwrite = TRUE)





# Set working directory
setwd("~/Downloads")

# 1. Dynamic Total Connectedness
pdf("Figure1_Dynamic_TCI.pdf", width = 12, height = 8)
PlotTCI(tvp_var_precovid, 
        ylim = c(0,80))  # Removed main parameter
dev.off()

# 2. Net Total Directional Connectedness
pdf("Figure2_NET.pdf", width = 12, height = 8)
PlotNET(tvp_var_precovid, 
        ylim = c(-40,20))
dev.off()

# 3. Net Pairwise Directional Connectedness
pdf("Figure3_NPDC.pdf", width = 12, height = 8)
PlotNPDC(tvp_var_precovid, 
         ylim = c(-20,10))
dev.off()