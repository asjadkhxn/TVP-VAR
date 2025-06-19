# README: Generating Results from the Gabauer_Repl.R Code File

Steps:

1. Loading data, ensuring proper date formatting, and transforming variables as needed (e.g., calculating log returns for stationarity).

2. Computing descriptive statistics (mean, variance, skewness, kurtosis) and running JB, ERS, and Ljung-Box tests.

3. Using the ConnectednessApproach package for setting up a TVP-VAR model and conducting dynamic spillover analysis.

4. Calculating TCI (overall connectedness), directional spillovers (identifying transmitters/receivers), and pairwise spillovers (measuring bilateral effects).

5. Applying a rolling window for tracking changes over time, especially during key events.

6. Generating plots for TCI, directional spillovers, and pairwise spillovers using the package’s visualization tools.


# Generating Results from the spillover_analysis.R Code File

Steps:

1. Loading data from an Excel file, extracting date and volatility columns from each sheet, handling duplicate dates by averaging, and combining all sheets into a single dataset. Pivoting the data into a wide format for analysis.

2. Computing basic statistics like mean, variance, and skewness, and running tests for normality (Jarque-Bera) and autocorrelation (Ljung-Box).

3. Using the TVP-VAR framework for estimating dynamic relationships between commodities, capturing spillover effects.

4. Calculating total connectedness, directional spillovers for identifying transmitters and receivers, and pairwise spillovers for measuring bilateral relationships.

5. Applying a rolling window for tracking changes over time and analyzing responses to external shocks.

6. Generating time-series plots for connectedness measures and saving results as Excel files or visualizations in PDF/PNG formats.


