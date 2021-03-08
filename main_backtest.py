import functions
import WaveletsMultiDim

import numpy as np
import matplotlib.pyplot as plt


# Retrieve Historical Data
ticker_id = "AAPL"
ohlc = functions.get_data(ticker=ticker_id, interval="1d", start_date="2019-01-01", end_date="2020-01-01")
train_signal = np.flip(ohlc['Close'].values)
ohlc = functions.get_data(ticker=ticker_id, interval="1d", start_date="2020-01-01", end_date="2021-01-01")
test_signal = np.flip(ohlc['Close'].values)
test_dates = np.flip(ohlc['Date'].values)

# Compute Returns and Lagged Returns
train_returns = functions.compute_returns(train_signal)
[train_returns, train_returns_lagged] = functions.lag_returns(train_returns)
test_returns = functions.compute_returns(test_signal)


# Compute Bivariate Gaussian Kernel Density Estimation
density_kernel_x, density_kernel_y, density_kernel_z = functions.kde2D(train_returns, train_returns_lagged, 0.01)
density_kernel_x, density_kernel_y = functions.update_axis_arrays(density_kernel_x, density_kernel_y)

# Backtesting (Gaussian Kernel)
cpt = 0
portfolio_value = np.ones(len(test_returns))
for i in range(0, len(test_returns)-1):
    predicted_return = functions.compute_VaR_2D(density_kernel_x, density_kernel_y, density_kernel_z, test_returns[i], 0.5)
    if predicted_return >= 0:
        portfolio_value[i + 1] = portfolio_value[i] * (1 + test_returns[i + 1])
    else:
        portfolio_value[i + 1] = portfolio_value[i] * (1 - test_returns[i + 1])
    if (predicted_return > 0 and test_returns[i+1] > 0) or (predicted_return <= 0 and test_returns[i+1] <= 0):
        cpt = cpt + 1
print("\n[INFO] {} - Gaussian Kernel Backtesting Results: {}%".format(functions.get_now(), round(100*cpt/(len(test_returns)-1), 2)))
plt.plot(test_dates[:-1], portfolio_value, label="Kernel Gaussian")


# Wavelet Parameters
nb_moments = 6
j = 6
wavelet = WaveletsMultiDim.Daubechie(signals=[train_returns, train_returns_lagged], nb_vanishing_moments=nb_moments)

# Compute Bivariate Wavelet Density Estimation
[density_wav_x, density_wav_y, density_wav_z] = wavelet.density(j, "linear")
density_wav_x, density_wav_y = functions.update_axis_arrays(density_wav_x, density_wav_y, type="normal")

# Backtesting (Wavelet)
cpt = 0
error_cpt = 0
portfolio_value = np.ones(len(test_returns))
for i in range(0, len(test_returns)-1):
    try:
        predicted_return = functions.compute_VaR_2D(density_wav_x, density_wav_y, density_wav_z, test_returns[i], 0.5)
        if predicted_return >= 0:
            portfolio_value[i + 1] = portfolio_value[i] * (1 + test_returns[i + 1])
        else:
            portfolio_value[i + 1] = portfolio_value[i] * (1 - test_returns[i + 1])
        if (predicted_return > 0 and test_returns[i+1] > 0) or (predicted_return <= 0 and test_returns[i+1] <= 0):
            cpt = cpt + 1
    except TypeError:
        error_cpt = error_cpt + 1
        portfolio_value[i + 1] = portfolio_value[i]
print("\n[INFO] {} - Wavelet Backtesting Results: {}%".format(functions.get_now(), round(100*cpt/(len(test_returns)-1), 2)))
print("NB ERROR : {}".format(error_cpt))
plt.plot(test_dates[:-1], portfolio_value, label="Wavelet")

# Show Results
plt.title("Portfolio Returns Evolution - " + ticker_id)
plt.plot(test_dates[:-1], test_signal[:-1]/test_signal[0], label="Stock")
plt.xlabel("Date")
plt.ylabel("Returns")
plt.legend()
plt.grid()
plt.show()
