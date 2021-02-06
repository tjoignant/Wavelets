import math
import numpy as np
import matplotlib.pyplot as plt

import functions
import Wavelets


# VaR Parameters
alpha_array = [0.01, 0.05, 0.10]

# Gaussian Kernel Smoothing Parameters
smoothing_array = [0.002, 0.005, 0.008]

# Retrieve Historical Data
ohlc = functions.get_data(ticker="AAPL", interval="1d", start_date="2019-01-01", end_date="2021-01-01")
signal = ohlc['Close']

# Compute Real & Estimated Returns
real_returns = functions.compute_returns(signal)
simulated_returns = np.random.normal(0, 1, 1750) / 100

# Compute Normal Density
density_normal_x = np.linspace(min(simulated_returns)*100, max(simulated_returns)*100, len(simulated_returns))
density_normal_y = np.empty(len(density_normal_x))
for i in range(0, len(density_normal_x)):
    density_normal_y[i] = functions.normal_pdf(density_normal_x[i], 0, 1)
density_normal_y = density_normal_y / sum(density_normal_y)
density_normal_x = density_normal_x/100

# Initializing Graph Params
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
fig.set_size_inches(16, 7.6)
fig.set_tight_layout(True)
fig.suptitle(str("Returns Density Estimation"), fontsize=16)

# Simulated Data - Non Parametric Density Estimation (Gaussian Kernel)
print("\n[INFO] {} - Computing the simulated returns density using the gaussian kernel estimator".format(functions.get_now()))
ax1.plot(density_normal_x, density_normal_y, label="Normal Density", linestyle='dashed')
for smoothing in smoothing_array:
    [density_kernel_x, density_kernel_y] = functions.compute_kernel_density(simulated_returns, smoothing=smoothing)
    functions.initialize_subplot(ax=ax1, title='Simulated Returns - Gaussian Kernel Estimator', xlabel="Returns")
    ax1.plot(density_kernel_x, density_kernel_y, label="Smoothing h={}".format(smoothing))
    ax1.legend(loc="upper right")

# Simulated Data - Wavelet Parameters
nb_moments = 6
j = 7
# Simulated Data - Wavelet Density Estimations (Wavelet Kernel & Thresholded)
wavelet = Wavelets.Daubechie(signal=simulated_returns, nb_vanishing_moments=nb_moments)
[density_wav_x, density_wav_linear_y] = wavelet.density(j, "linear")
[density_wav_x, density_wav_donoho_y] = wavelet.density(j, "donoho")
functions.initialize_subplot(ax=ax2, title='Simulated Returns - Wavelet Estimators (db_{} / j={})'.format(nb_moments, j), xlabel="Returns")
ax2.plot(density_normal_x, density_normal_y, label="Normal Density", linestyle='dashed')
ax2.plot(density_wav_x, density_wav_linear_y, label="Linear")
ax2.plot(density_wav_x, density_wav_donoho_y, label="Donoho")
ax2.legend(loc="upper right")

# Simulated Data - VaR
print("\nSIMULATED RETURNS:")
smoothing = 0.002
[density_kernel_x, density_kernel_y] = functions.compute_kernel_density(simulated_returns, smoothing=smoothing)
for alpha in alpha_array:
    print("    - Normal density VaR {}%: {}%".format(int((1-alpha)*100), functions.compute_VaR(density_normal_x, density_normal_y, alpha)))
    print("    - Kernel density (h: {}) VaR {}%: {}%".format(smoothing, int((1 - alpha)*100), functions.compute_VaR(density_kernel_x, density_kernel_y, alpha)))
    print("    - Wavelet Linear density VaR {}%: {}%".format(int((1 - alpha)*100), functions.compute_VaR(density_wav_x, density_wav_linear_y, alpha)))
    print("    - Wavelet Donoho density VaR {}%: {}%\n".format(int((1 - alpha)*100), functions.compute_VaR(density_wav_x, density_wav_donoho_y, alpha)))

print("---------------------------------------------------\n")

# Real Data - Non Parametric Density Estimation (Gaussian Kernel)
print("[INFO] {} - Computing the real returns density using the gaussian kernel estimator".format(functions.get_now()))
for smoothing in smoothing_array:
    [density_kernel_x, density_kernel_y] = functions.compute_kernel_density(real_returns, smoothing=smoothing)
    functions.initialize_subplot(ax=ax3, title='Real Returns - Gaussian Kernel Estimator', xlabel="Returns")
    ax3.plot(density_kernel_x, density_kernel_y, label="Smoothing h={}".format(smoothing))
    ax3.legend(loc="upper right")

# Wavelet Parameters
nb_moments = 6
j = 5
# Real Data - Wavelet Density Estimations (Wavelet Kernel & Thresholded)
wavelet = Wavelets.Daubechie(signal=real_returns, nb_vanishing_moments=nb_moments)
[density_wav_x, density_wav_linear_y] = wavelet.density(j, "linear")
[density_wav_x, density_wav_donoho_y] = wavelet.density(j, "donoho")
functions.initialize_subplot(ax=ax4, title='Real Returns - Wavelet Estimators (db_{} / j={})'.format(nb_moments, j), xlabel="Returns")
ax4.plot(density_wav_x, density_wav_linear_y, label="Linear")
ax4.plot(density_wav_x, density_wav_donoho_y, label="Donoho")
ax4.legend(loc="upper right")

# Real Data - VaR
print("\nREAL RETURNS:")
smoothing = 0.008
[density_kernel_x, density_kernel_y] = functions.compute_kernel_density(real_returns, smoothing=smoothing)
for alpha in alpha_array:
    print("    - Kernel density (h: {}) VaR {}%: {}%".format(smoothing, int((1 - alpha)*100), functions.compute_VaR(density_kernel_x, density_kernel_y, alpha)))
    print("    - Wavelet Linear density VaR {}%: {}%".format(int((1 - alpha)*100), functions.compute_VaR(density_wav_x, density_wav_linear_y, alpha)))
    print("    - Wavelet Donoho density VaR {}%: {}%\n".format(int((1 - alpha)*100), functions.compute_VaR(density_wav_x, density_wav_donoho_y, alpha)))

# Show Graph
plt.show()
