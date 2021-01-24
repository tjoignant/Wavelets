import numpy as np
import matplotlib.pyplot as plt

import functions
import Wavelets


# Retrieve Historical Data
ohlc = functions.get_data("AAPL", "1d")
signal = ohlc['Close']

# Choosing Analysis Parameters
j_max = -4
nb_moments = 4
threshold = "SURE"

# Initializing Graph Params
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
fig.set_size_inches(16, 7.6)
fig.set_tight_layout(True)
if nb_moments == 1:
    wavelet_name = "Haar"
else:
    wavelet_name = "Daubechie_{}".format(nb_moments)
fig.suptitle(str("Signal Wavelet Analysis (" + wavelet_name + ")"), fontsize=16)

# Wavelet Analysis
y = np.empty(len(signal))
y_denoised = np.empty(len(signal))
wavelet = Wavelets.Daubechie(signal=signal, nb_vanishing_moments=nb_moments)
for level in range(0, abs(j_max)):
    j = -1 - level
    # Compute Approximation & Details Coefficients (Discrete Wavelet Transformation)
    [cA, cD] = wavelet.dwt(j)
    # Compute Approximation & Details Values (Inverse Discrete Wavelet Transformation)
    [yA, yD] = wavelet.idwt(cA, cD, j)
    # Compute Denoised Approximation & Details Coefficients
    cD_denoised = wavelet.denoise(cD, threshold=threshold)
    # Compute Denoised Approximation & Details Values (Inverse Discrete Wavelet Transformation)
    [yA_denoised, yD_denoised] = wavelet.idwt(cA, cD_denoised, j)
    # Compute Reconstructed Signals
    y = y + yD
    y_denoised = y_denoised + yD_denoised
    if j == j_max:
        y = y + yA
        y_denoised = y_denoised + yA_denoised
    # Plot Details Values
    functions.initialize_subplot(ax=ax1, title='Details', xlabel="Time Series", ylabel="cD Values")
    ax1.plot(yD, label="Scale j={}".format(j))
    ax1.legend(loc="upper right")
    # Plot Denoised Details Values
    functions.initialize_subplot(ax=ax3, title='Denoised Details', xlabel="Time Series", ylabel="cD_denoised Values")
    ax3.plot(yD_denoised, label="Scale j={}".format(j))
    ax3.legend(loc="upper right")

# Plot Reconstructed Signal
functions.initialize_subplot(ax=ax2, title='Signals', xlabel="Time Series", ylabel='Stock Price')
ax2.plot(signal, label="Original Signal")
ax2.plot(y, label="Reconstructed Signal")
ax2.legend(loc="upper right")

# Plot Reconstructed Denoised Signal
functions.initialize_subplot(ax=ax4, title='Denoised Signals', xlabel="Time Series", ylabel='Stock Price')
ax4.plot(signal, label="Original Signal")
ax4.plot(y_denoised, label="Reconstructed Denoised Signal")
ax4.legend(loc="upper right")

# Show Graph
plt.show()
