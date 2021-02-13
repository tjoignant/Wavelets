import numpy as np
import matplotlib.pyplot as plt

import functions
import WaveletsMultiDim


# VaR Parameters
alpha_array = [0.01, 0.05, 0.10]

# Retrieve Historical Data
ohlc = functions.get_data(ticker="AAPL", interval="1d", start_date="2020-01-01", end_date="2021-01-01")
signal = ohlc['Close']

# Compute Real & Simulated Returns
real_returns = functions.compute_returns(signal)
real_returns_lagged = real_returns[1:len(real_returns)]
real_returns = real_returns[0:-1]
simulated_returns = np.random.normal(0, 1, 250) / 100
simulated_returns_lagged = simulated_returns[1:len(simulated_returns)]
simulated_returns = simulated_returns[0:-1]

# Initializing Graph Params
fig = plt.figure(figsize=plt.figaspect(0.5))
fig.set_size_inches(16, 7.6)
fig.set_tight_layout(True)
fig.suptitle(str("Returns Density Estimation"), fontsize=16)

# Compute Normal Bivariate Density
density_normal_X = np.linspace(min(simulated_returns)*100, max(simulated_returns)*100, len(simulated_returns))
density_normal_Y = np.linspace(min(simulated_returns_lagged)*100, max(simulated_returns_lagged)*100, len(simulated_returns_lagged))
density_normal_x, density_normal_y = np.meshgrid(density_normal_X, density_normal_Y)
density_normal_z = functions.normal_pdf_2D(density_normal_x, density_normal_y)
density_normal_x = density_normal_x/100
density_normal_y = density_normal_y/100
ax1 = fig.add_subplot(2, 3, 1, projection='3d')
functions.initialize_3D_subplot(ax=ax1, title="Simulated Returns", xlabel="Returns", ylabel="Lagged Returns")
ax1.plot_surface(density_normal_x, density_normal_y, density_normal_z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
density_normal_x, density_normal_y = functions.update_axis_arrays(density_normal_x, density_normal_y, type="normal")

# Simulated Data - Non Parametric Density Estimation (Gaussian Kernel)
print("\n[INFO] {} - Computing the simulated returns bivariate density using the gaussian kernel estimator".format(functions.get_now()))
density_kernel_x, density_kernel_y, density_kernel_z = functions.kde2D(simulated_returns, simulated_returns_lagged, 0.01)
ax2 = fig.add_subplot(2, 3, 2, projection='3d')
functions.initialize_3D_subplot(ax=ax2, title="Gaussian Kernel Estimator", xlabel="Returns", ylabel="Lagged Returns")
ax2.plot_surface(density_kernel_x, density_kernel_y, density_kernel_z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
density_kernel_x, density_kernel_y = functions.update_axis_arrays(density_kernel_x, density_kernel_y, type="kernel")

# Simulated Data - Wavelet Parameters
nb_moments = 6
j = 6
# Simulated Data - Wavelet Density Estimations (Wavelet Kernel & Thresholded)
wavelet = WaveletsMultiDim.Daubechie(signals=[simulated_returns, simulated_returns_lagged], nb_vanishing_moments=nb_moments)
[density_wav_x, density_wav_y, density_wav_z] = wavelet.density(j, "linear")
ax3 = fig.add_subplot(2, 3, 3, projection='3d')
functions.initialize_3D_subplot(ax=ax3, title='Wavelet Estimators (db_{} / j={})'.format(nb_moments, j), xlabel="Returns", ylabel="Lagged Returns")
ax3.plot_surface(density_wav_x, density_wav_y, density_wav_z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
density_wav_x, density_wav_y = functions.update_axis_arrays(density_wav_x, density_wav_y, type="normal")

# Simulated Data - VaR
print("\nSIMULATED RETURNS:")
smoothing = 0.002
last_return = 0
for alpha in alpha_array:
    print("    - Normal density VaR {}%: {}%".format(int((1-alpha)*100), functions.compute_VaR_2D(density_normal_x, density_normal_y, density_normal_z, last_return, alpha)))
    print("    - Kernel density (h: {}) VaR {}%: {}%".format(smoothing, int((1 - alpha)*100), functions.compute_VaR_2D(density_kernel_x, density_kernel_y, density_kernel_z, last_return, alpha)))
    print("    - Wavelet Linear density VaR {}%: {}%\n".format(int((1 - alpha)*100), functions.compute_VaR_2D(density_wav_x, density_wav_y, density_wav_z, last_return, alpha)))

print("---------------------------------------------------\n")

# Real Data - Non Parametric Density Estimation (Gaussian Kernel)
print("[INFO] {} - Computing the real returns bivariate density using the gaussian kernel estimator".format(functions.get_now()))
density_kernel_x, density_kernel_y, density_kernel_z = functions.kde2D(real_returns, real_returns_lagged, 0.01)
ax3 = fig.add_subplot(2, 3, 5, projection='3d')
functions.initialize_3D_subplot(ax=ax3, title="Real Returns - Gaussian Kernel Estimator", xlabel="Returns", ylabel="Lagged Returns")
ax3.plot_surface(density_kernel_x, density_kernel_y, density_kernel_z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
density_kernel_x, density_kernel_y = functions.update_axis_arrays(density_kernel_x, density_kernel_y, type="kernel")

# Wavelet Parameters
nb_moments = 6
j = 6
# Real Data - Wavelet Density Estimations (Wavelet Kernel & Thresholded)
wavelet = WaveletsMultiDim.Daubechie(signals=[real_returns, real_returns_lagged], nb_vanishing_moments=nb_moments)
[density_wav_x, density_wav_y, density_wav_z] = wavelet.density(j, "linear")
ax6 = fig.add_subplot(2, 3, 6, projection='3d')
functions.initialize_3D_subplot(ax=ax6, title='Wavelet Estimators (db_{} / j={})'.format(nb_moments, j), xlabel="Returns", ylabel="Lagged Returns")
ax6.plot_surface(density_wav_x, density_wav_y, density_wav_z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
density_wav_x, density_wav_y = functions.update_axis_arrays(density_wav_x, density_wav_y, type="normal")

# Real Data - VaR
print("\nREAL RETURNS:")
smoothing = 0.008
last_return = 0
for alpha in alpha_array:
    print("    - Kernel density (h: {}) VaR {}%: {}%".format(smoothing, int((1 - alpha)*100), functions.compute_VaR_2D(density_kernel_x, density_kernel_y, density_kernel_z, last_return, alpha)))
    print("    - Wavelet Linear density VaR {}%: {}%".format(int((1 - alpha)*100), functions.compute_VaR_2D(density_wav_x, density_wav_y, density_wav_z, last_return, alpha)))

# Show Graph
plt.show()
