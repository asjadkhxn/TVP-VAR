#TVP-VAR

# Install and load required packages
library(readxl)
library(tidyverse)
library(zoo)

# Function to read and process a single sheet
process_sheet <- function(filename, sheet_name) {
  # Read the sheet
  df <- read_excel(filename, sheet = sheet_name)
  
  # Extract date and volatility columns (4th and 5th columns)
  data.frame(
    date = as.Date(as.numeric(df[[4]]), origin = "1899-12-30"),
    volatility = as.numeric(df[[5]]),
    commodity = sheet_name
  ) %>%
    # Handle duplicate dates by taking the mean
    group_by(date, commodity) %>%
    summarize(volatility = mean(volatility, na.rm = TRUE), .groups = "drop")
}

# Main processing function
process_data <- function(filename) {
  # Get sheet names
  sheets <- excel_sheets(filename)
  sheets <- sheets[sheets != "INDEX"]  # Remove INDEX sheet
  
  # Process each sheet and combine
  all_data <- map_df(sheets, ~process_sheet(filename, .x))
  
  # Create wide format
  wide_data <- all_data %>%
    pivot_wider(
      id_cols = date,
      names_from = commodity,
      values_from = volatility
    )
  
  return(wide_data)
}

# Read and process the data
wide_data <- process_data("Raw_Data.xlsx")

# Basic data checks
print("Data Structure:")
str(wide_data)

print("\nDate Range:")
print(paste("Date range:", min(wide_data$date), "to", max(wide_data$date)))

print("\nMissing Values by Column:")
print(colSums(is.na(wide_data)))

print("\nFirst Few Rows:")
print(head(wide_data))

# Save processed data
write.csv(wide_data, "processed_volatility_data.csv", row.names = FALSE)

# Calculate correlation matrix
cor_matrix <- cor(wide_data[,-1], use = "pairwise.complete.obs")
print("\nCorrelation Matrix (first few columns):")
print(round(head(cor_matrix), 2))




# Part 1: Initial Setup and Data Preparation
selected_commodities <- c("CL", "HO", "NG",  # Energy
                          "C", "W", "S", "CC", "CT", "KC", "SB", "BO")  # Agricultural

clean_data <- wide_data %>%
  select(date, all_of(selected_commodities)) %>%
  na.omit()

Y <- as.matrix(clean_data[,-1])  # Remove date column
N <- ncol(Y)
T <- nrow(Y)

# Print dimensions and first few rows to verify
print("Data Dimensions:")
print(dim(Y))
print("\nFirst few rows of data:")
print(head(Y))



# Parameters
p <- 1  # VAR lag order
h <- 10 # Forecast horizon
N <- ncol(Y)  # Number of variables
T <- nrow(Y)  # Number of observations

# Create lagged data matrix (Yt-1)
Y_lag <- Y[-T,]  # Remove last row
Y_current <- Y[-1,]  # Remove first row

# Initialize arrays for time-varying parameters
beta_t <- array(0, dim = c(N, N, T-1))  # β_t coefficients
S_t <- array(0, dim = c(N, N, T-1))     # Σ_t covariance matrices

# Rolling window estimation
window_size <- 252  # One year of daily data
for(t in window_size:(T-1)) {
  # Get window data
  window_start <- t - window_size + 1
  Y_window <- Y[window_start:t,]
  Y_lag_window <- Y_window[-window_size,]
  Y_current_window <- Y_window[-1,]
  
  # Estimate β_t using OLS for each window
  beta_t[,,t] <- solve(t(Y_lag_window) %*% Y_lag_window) %*% 
    t(Y_lag_window) %*% Y_current_window
  
  # Calculate residuals and Σ_t
  residuals <- Y_current_window - Y_lag_window %*% beta_t[,,t]
  S_t[,,t] <- t(residuals) %*% residuals / (window_size - 1)
}

# Print sample of results to verify
print("Sample of beta_t coefficients (first window, first few elements):")
print(round(beta_t[1:3, 1:3, window_size], 4))

print("\nSample of S_t covariance matrix (first window, first few elements):")
print(round(S_t[1:3, 1:3, window_size], 4))



# Parameters from previous step remain the same
N <- ncol(Y)
h <- 10  # forecast horizon

# Function to calculate h-step error variance (ψ²ij,t) and GFEVD (θ̃ij(h))
calculate_variance_gfevd <- function(beta, Sigma, h) {
  # Initialize arrays
  psi_sq <- array(0, dim = c(N, N, h-1))  # ψ²ij,t
  theta <- array(0, dim = c(N, N, h))      # θ̃ij(h)
  
  # Calculate ψ²ij,t following the formula in the image
  for(l in 1:(h-1)) {
    for(i in 1:N) {
      for(j in 1:N) {
        # Calculate S^(-1/2)
        S_inv_sqrt <- Sigma[i,j]^(-1/2)
        
        # Calculate Ah,t * εj,t
        if(l == 1) {
          A_h <- beta
        } else {
          A_h <- A_h %*% beta
        }
        
        # Calculate ψ²ij,t using the formula from the image
        psi_sq[i,j,l] <- (S_inv_sqrt * sum(A_h[i,] * A_h[j,]))^2
      }
    }
  }
  
  # Calculate GFEVD (θ̃ij(h)) using the formula from the image
  for(i in 1:N) {
    for(j in 1:N) {
      numerator <- sum(psi_sq[i,j,])
      denominator <- sum(sapply(1:N, function(k) sum(psi_sq[i,k,])))
      theta[i,j,] <- rep(numerator/denominator, h)
    }
  }
  
  return(list(psi_sq = psi_sq, theta = theta))
}

# Calculate for the last window
variance_gfevd <- calculate_variance_gfevd(
  beta_t[,,dim(beta_t)[3]], 
  S_t[,,dim(S_t)[3]], 
  h
)

# Print sample results
print("Sample of ψ²ij,t (first few elements):")
print(round(variance_gfevd$psi_sq[1:3, 1:3, 1], 4))

print("\nSample of θ̃ij(h) (first few elements):")
print(round(variance_gfevd$theta[1:3, 1:3, 1], 4))



#Data Analysis Part



# Calculate directional connectedness measures
calculate_connectedness <- function(theta) {
    N <- dim(theta)[1]
    
    # Initialize matrices for directional connectedness
    C_to <- matrix(0, N, N)    # C^g_{i←j,t}(h)
    C_from <- matrix(0, N, N)  # C^g_{i→j,t}(h)
    C_net <- matrix(0, N, N)   # C^n_{i,j}(h)
    
    # Calculate directional connectedness following the formulas in the image
    for(i in 1:N) {
        for(j in 1:N) {
            if(i != j) {
                # To others (C^g_{i←j,t}(h))
                C_to[i,j] <- (sum(theta[i,j,]) / sum(theta[i,,])) * 100
                
                # From others (C^g_{i→j,t}(h))
                C_from[i,j] <- (sum(theta[j,i,]) / sum(theta[j,,])) * 100
            }
        }
    }
    
    # Calculate net directional connectedness (C^n_{i,j}(h))
    for(i in 1:N) {
        for(j in 1:N) {
            if(i != j) {
                C_net[i,j] <- C_to[i,j] - C_from[i,j]
            }
        }
    }
    
    return(list(to = C_to, from = C_from, net = C_net))
}

# Calculate NPDC following the final equation in the image
calculate_npdc <- function(theta) {
    N <- dim(theta)[1]
    npdc <- matrix(0, N, N)
    
    for(i in 1:N) {
        for(j in 1:N) {
            if(i != j) {
                # Following the NPDC formula: (θ̃ij(h) - θ̃ji(h))/N × 100
                npdc[i,j] <- (sum(theta[i,j,]) - sum(theta[j,i,])) * 100 / N
            }
        }
    }
    
    return(npdc)
}

# Calculate all measures
connectedness <- calculate_connectedness(variance_gfevd$theta)
npdc_result <- calculate_npdc(variance_gfevd$theta)

# Add labels
colnames(npdc_result) <- colnames(Y)
rownames(npdc_result) <- colnames(Y)

# Print results
print("Directional Connectedness TO others (first few elements):")
print(round(connectedness$to[1:3, 1:3], 4))

print("\nDirectional Connectedness FROM others (first few elements):")
print(round(connectedness$from[1:3, 1:3], 4))

print("\nNet Directional Connectedness (first few elements):")
print(round(connectedness$net[1:3, 1:3], 4))

print("\nNet Pairwise Directional Connectedness (NPDC) (first few elements):")
print(round(npdc_result[1:3, 1:3], 4))

# Create visualization of NPDC
library(ggplot2)
npdc_df <- as.data.frame(as.table(npdc_result))
names(npdc_df) <- c("From", "To", "Value")

ggplot(npdc_df, aes(x = From, y = To, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0) +
    theme_minimal() +
    labs(title = "Net Pairwise Directional Connectedness",
         x = "From", y = "To") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("npdc_heatmap.png", width = 10, height = 8)


# Install and load required packages
install.packages("gridExtra")
library(tidyverse)
library(ggplot2)
library(gridExtra)

# 1. Create Enhanced Summary Statistics
summary_stats <- clean_data %>%
  select(-date) %>%
  summarise(across(everything(), 
                   list(
                     Mean = ~mean(., na.rm = TRUE),
                     SD = ~sd(., na.rm = TRUE),
                     Min = ~min(., na.rm = TRUE),
                     Max = ~max(., na.rm = TRUE),
                     CV = ~sd(., na.rm = TRUE)/mean(., na.rm = TRUE)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = "Commodity",
               values_to = "Value") %>%
  separate(Commodity, into = c("Commodity", "Statistic"), sep = "_") %>%
  pivot_wider(names_from = Statistic, values_from = Value)

# 2. Create Time Series Plot with Market Classification
ts_plot <- clean_data %>%
  pivot_longer(-date, names_to = "Commodity", values_to = "Volatility") %>%
  mutate(Market = case_when(
    Commodity %in% c("CL", "HO", "NG") ~ "Energy",
    TRUE ~ "Agricultural"
  )) %>%
  ggplot(aes(x = date, y = Volatility, color = Market)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~Commodity, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Time Series of Commodity Volatilities",
       x = "Date",
       y = "Volatility") +
  theme(legend.position = "bottom")

# 3. Create Correlation Matrix Heatmap
cor_matrix <- cor(clean_data[,-1], use = "complete.obs")
cor_df <- as.data.frame(as.table(cor_matrix))
names(cor_df) <- c("Var1", "Var2", "Correlation")

cor_plot <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix of Commodity Volatilities")

# 4. Create Box Plot with Market Classification
box_plot <- clean_data %>%
  pivot_longer(-date, names_to = "Commodity", values_to = "Volatility") %>%
  mutate(Market = case_when(
    Commodity %in% c("CL", "HO", "NG") ~ "Energy",
    TRUE ~ "Agricultural"
  )) %>%
  ggplot(aes(x = Commodity, y = Volatility, fill = Market)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Commodity Volatilities")

# Save all plots
ggsave("time_series_plot.png", ts_plot, width = 15, height = 10)
ggsave("correlation_matrix.png", cor_plot, width = 10, height = 8)
ggsave("boxplot.png", box_plot, width = 10, height = 6)

# Print summary statistics
print("Summary Statistics:")
print(summary_stats)

# Save summary statistics
write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)