import functions
import Wavelets


# # Retrieve Historical Data
ohlc = functions.get_data(ticker="AAPL", interval="1d", start_date="2018-01-01", end_date="2020-01-01")
train_signal = ohlc['Close']
ohlc = functions.get_data(ticker="AAPL", interval="1d", start_date="2020-01-01", end_date="2021-01-01")
test_signal = ohlc['Close']

# Compute Returns and Lagged Returns
train_returns = functions.compute_returns(train_signal)
[train_returns, train_returns_lagged] = functions.lag_returns(train_returns)

# Compute Bivariate Density Estimation
density_kernel_x, density_kernel_y, density_kernel_z = functions.kde2D(train_returns, train_returns_lagged, 0.01)
density_kernel_x, density_kernel_y = functions.update_axis_arrays(density_kernel_x, density_kernel_y)

# Backtesting
test_returns = functions.compute_returns(test_signal)
cpt = 0
for i in range(0, len(test_returns)-1):
    predicted_return = functions.compute_VaR_2D(density_kernel_x, density_kernel_y, density_kernel_z, test_returns[i], 0.5)
    if predicted_return > 0 and test_returns[i+1] > 0 or predicted_return <= 0 and test_returns[i+1] <= 0:
        cpt = cpt + 1
print("\n[INFO] {} - Backtesting Results: {}%".format(functions.get_now(), round(100*cpt/(len(test_returns)-1), 2)))
