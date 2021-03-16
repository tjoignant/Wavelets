import functions
import WaveletsMultiDim

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad


# Retrieve Historical Data
ticker_id = "^GSPC"
ohlc = functions.get_data(ticker=ticker_id, interval="1d", start_date="2011-09-02", end_date="2013-09-02")
train_signal = np.flip(ohlc['Close'].values)
ohlc = functions.get_data(ticker=ticker_id, interval="1d", start_date="2013-09-03", end_date="2015-04-17")
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
pKernel = cpt
plt.plot(test_dates[:-1], test_signal[:-1]/test_signal[0], label="Ticker")
plt.plot(test_dates[:-1], portfolio_value, label="Kernel Portfolio")

returns = functions.compute_returns(portfolio_value)
risk_free_rate = 0.02
annualized_returns = np.mean(returns) * 252 * 100
volatility = returns.std()
annualized_vol = volatility * np.sqrt(252) * 100
sharpe = (annualized_returns - risk_free_rate) / annualized_vol
VaR_95 = functions.compute_empirical_VaR(returns, 0.95) * 100
VaR_99 = functions.compute_empirical_VaR(returns, 0.99) * 100
lower_vol = functions.compute_semi_deviation(returns)
annualized_lower_vol = lower_vol * np.sqrt(252) * 100
sortino = (annualized_returns - risk_free_rate) / annualized_lower_vol
[i, j] = functions.compute_MDD(portfolio_value)
MDD = (portfolio_value[i] - portfolio_value[j]) * 100
skewness = functions.compute_skewness(returns)
kurtosis = functions.compute_kurtosis(returns)
n = len(portfolio_value)
p0 = 0.5
pKernel /= n
z0Kernel = (pKernel - p0) / np.sqrt(p0*(1-p0) / n)
integrateKernel, err1 = quad(functions.normalfunction, -1000, z0Kernel)
p_valueKernel = 1 - integrateKernel
print("\nKERNEL GAUSSIAN")
print("Backtesting Results: {}%".format(round(100*cpt/(len(test_returns)-1), 2)))
print("P-value: {}".format(round(p_valueKernel, 2)))
print("Returns : {}%".format(round(annualized_returns, 2)))
print("Vol : {}%".format(round(annualized_vol, 2)))
print("Lower Vol : {}%".format(round(annualized_lower_vol, 2)))
print("Sharpe : {}".format(round(sharpe, 2)))
print("Sortino : {}".format(round(sortino, 2)))
print("VaR 95% : {}%".format(round(VaR_95, 2)))
print("VaR 99% : {}%".format(round(VaR_99, 2)))
print("MDD : {}%".format(round(MDD, 2)))
print("Skewness : {}".format(round(skewness, 3)))
print("Kurtosis: {}".format(round(kurtosis, 3)))
plt.plot([test_dates[i], test_dates[j]], [portfolio_value[i], portfolio_value[j]], 'o', color='Red', markersize=5, label="Kernel MDD")


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

pBivariate = cpt
print("NB ERROR : {}".format(error_cpt))
plt.plot(test_dates[:-1], portfolio_value, label="Wavelet Portfolio")

# p-value
n = len(portfolio_value)
p0 = 0.5
pBivariate /= n
z0Bivariate = (pBivariate - p0) / np.sqrt(p0*(1-p0) / n)
integrateBivariate, err2 = quad(functions.normalfunction, -1000, z0Bivariate)
p_valueBivariate = 1 - integrateBivariate
print("\nWAVELET")
print("Backtesting Results: {}%".format(round(100*cpt/(len(test_returns)-1), 2)))
print("P-value: {}".format(round(p_valueBivariate, 2)))
returns = functions.compute_returns(portfolio_value)
risk_free_rate = 0.02
annualized_returns = np.mean(returns) * 252 * 100
volatility = returns.std()
annualized_vol = volatility * np.sqrt(252) * 100
sharpe = (annualized_returns - risk_free_rate) / annualized_vol
VaR_95 = functions.compute_empirical_VaR(returns, 0.95) * 100
VaR_99 = functions.compute_empirical_VaR(returns, 0.99) * 100
lower_vol = functions.compute_semi_deviation(returns)
annualized_lower_vol = lower_vol * np.sqrt(252) * 100
sortino = (annualized_returns - risk_free_rate) / annualized_lower_vol
[i, j] = functions.compute_MDD(portfolio_value)
MDD = (portfolio_value[i] - portfolio_value[j]) * 100
skewness = functions.compute_skewness(returns)
kurtosis = functions.compute_kurtosis(returns)
print("Returns : {}%".format(round(annualized_returns, 2)))
print("Vol : {}%".format(round(annualized_vol, 2)))
print("Lower Vol : {}%".format(round(annualized_lower_vol, 2)))
print("Sharpe : {}".format(round(sharpe, 2)))
print("Sortino : {}".format(round(sortino, 2)))
print("VaR 95% : {}%".format(round(VaR_95, 2)))
print("VaR 99% : {}%".format(round(VaR_99, 2)))
print("MDD : {}%".format(round(MDD, 2)))
print("Skewness : {}".format(round(skewness, 3)))
print("Kurtosis: {}".format(round(kurtosis, 3)))
plt.plot([test_dates[i], test_dates[j]], [portfolio_value[i], portfolio_value[j]], 'o', color='Green', markersize=5, label="Wavelet MDD")

# Show Results
plt.title("Portfolio Returns Evolution - " + ticker_id)
plt.xlabel("Date")
plt.ylabel("Returns")
plt.legend()
plt.grid()

plt.show()
