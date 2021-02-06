import numpy as np
import matplotlib.pyplot as plt

import functions
import Wavelets


# Retrieve Historical Data
ohlc = functions.get_data(ticker="AAPL", interval="1d", start_date="2019-01-01", end_date="2021-01-01")
signal = ohlc['Close']

# Choosing Analysis Parameters
j_max = -4
nb_moments = 4

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
wavelet = Wavelets.Daubechie(signal=signal, nb_vanishing_moments=nb_moments)
for level in range(0, abs(j_max)):
    j = -1 - level
    # Compute Approximation & Details coefficients (Discrete Wavelet Transformation)
    [cA, cD] = wavelet.dwt(j)
    # Compute Approximation & Details values (Inverse Discrete Wavelet Transformation)
    [yA, yD] = wavelet.idwt(cA, cD, j)
    # Compute Reconstructed Signals
    y = y + yD
    if j == j_max:
        y = y + yA
    # Plot Details Values
    functions.initialize_subplot(ax=ax1, title='Details', xlabel="Time Series", ylabel="cD Values")
    ax1.plot(yD, label="Scale j={}".format(j))
    ax1.legend(loc="upper right")
    # Plot Approximations Values
    functions.initialize_subplot(ax=ax3, title='Approximations', xlabel="Time Series", ylabel="Approximation Values")
    ax3.plot(yA, label="Scale j={}".format(j))
    ax3.legend(loc="upper right")

# Plot Reconstructed Signal
functions.initialize_subplot(ax=ax2, title='Signals', xlabel="Time Series", ylabel='Stock Price')
ax2.plot(signal, label="Original Signal")
ax2.plot(y, label="Reconstructed Signal")
ax2.legend(loc="upper right")

# Plot Father Wavelet
father_wavelet = wavelet.get_father_wavelet()
functions.initialize_subplot(ax=ax4, title='Father Wavelet', xlabel='Support')
ax4.plot(father_wavelet[0], father_wavelet[1])

# Show Graph
plt.show()
