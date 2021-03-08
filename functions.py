from sklearn.neighbors import KernelDensity

import datetime as dt
import yfinance as yf
import numpy as np
import math


def get_data(ticker, interval, start_date, end_date):
    """
    Download ticker's ohlcv data for a chosen time period and time interval (either "1d", "1h", "30m", "15m", "5m")
    :param ticker       ticker id [str]
    :param interval     interval type [str]
    :param start_date   interval first day [str]
    :param end_date     interval last day [str]
    :return ohlcv       open high low close volume dataframe [df]
    """
    # Display indication
    print('[INFO] {} - Retrieving {}_{} historical data'.format(get_now(), ticker, interval))
    # Download ticker's ohlcv
    ohlcv = yf.download(tickers=ticker, start=start_date, end=end_date, interval=interval)
    # Modify dataframe
    ohlcv.drop(columns=['Adj Close'], inplace=True)
    ohlcv.sort_index(axis=0, ascending=False, inplace=True)
    ohlcv.reset_index(inplace=True)
    if "Datetime" in ohlcv.columns:
        ohlcv['Datetime'] = ohlcv['Datetime'].astype(str).str[:-9]
    return ohlcv


def get_now():
    """
    Retrieve the current date and time
    :return:                current date and time [str]
    """
    now = dt.datetime.now()
    now_str = now.strftime("%d/%m %H:%M")
    return now_str


def initialize_subplot(ax, title="", xlabel="", ylabel=""):
    """
    Initialize graph subplot
    :param ax:              subplot's axis [2D matplotlib axis]
    :param title:           subplot's title [str]
    :param xlabel:          subplot's horizontal axis's label [str]
    :param ylabel:          subplot's vertical axis's label [str]
    """
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def initialize_3D_subplot(ax, title="", xlabel="", ylabel="", zlabel=""):
    """
    Initialize graph subplot
    :param ax:              subplot's axis [3D matplotlib axis]
    :param title:           subplot's title [str]
    :param xlabel:          subplot's horizontal 1 axis' label [str]
    :param ylabel:          subplot's horizontal 2 axis' label [str]
    :param zlabel:          subplot's vertical axis's label [str]
    """
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)


def compute_returns(prices, method="frac"):
    """
    Compute returns of given close prices using either the "frac" or "log" method
    :param prices:              close prices [1D array]
    :param method:              method used [str]
    :return:                    returns [1D array]
    """
    returns = np.zeros(len(prices))
    for i in range(1, len(prices)):
        if method == "log":
            returns[i] = math.log(prices[i]) - math.log(prices[i - 1])
        elif method == "frac":
            returns[i] = (prices[i] - prices[i - 1]) / prices[i - 1]
        else:
            return "Method undefined"
    returns[0] = None
    return returns[1:len(returns)]


def lag_returns(returns, lag=1):
    """
    Compute
    :param returns:         returns [1D array]
    :param lag:             number of days lagged [int]
    :return:                returns & lagged returns [2D array]
    """
    train_returns_lagged = returns[lag:len(returns)]
    train_returns = returns[0:-lag]
    return [train_returns, train_returns_lagged]


def compute_histogram_density(returns, nb_intervals):
    """
    Compute return's histogram density
    :param returns:                     returns array [1D array]
    :param nb_intervals:                number of histogram [int]
    :return: density index & values     returns density [2D array]
    """
    lowest_value = round(min(returns), 3)
    highest_value = round(max(returns), 3)
    hist_array_x = np.linspace(lowest_value, highest_value, nb_intervals)
    hist_array_y = np.zeros(len(hist_array_x))
    for i in range(0, len(returns)):
        for j in range(0, len(hist_array_x) - 1):
            if hist_array_x[j] < returns[i] <= hist_array_x[j + 1]:
                hist_array_y[j] = hist_array_y[j] + 1
                break
        if returns[i] > hist_array_x[len(hist_array_x) - 1]:
            hist_array_y[len(hist_array_x) - 1] = hist_array_y[len(hist_array_x) - 1] + 1
    hist_array_y = hist_array_y / sum(hist_array_y)
    return [hist_array_x, hist_array_y]


def compute_kernel_density(returns, smoothing):
    """
    Compute return's histogram density
    :param returns:                     returns array [1D array]
    :param smoothing:                   smoothing parameter [float]
    :return:                            density index & values [2D array]
    """
    lowest_value = round(min(returns), 3)
    highest_value = round(max(returns), 3)
    kernel_array_x = np.linspace(lowest_value, highest_value, len(returns))
    kernel_array_y = np.zeros(len(kernel_array_x))
    for i in range(0, len(kernel_array_x)):
        kernel_array_y[i] = kernel_pdf(returns, kernel_array_x[i], smoothing)
    kernel_array_y = kernel_array_y / sum(kernel_array_y)
    return [kernel_array_x, kernel_array_y]


def kernel_cdf(returns, x, smoothing):
    """
    Compute Gaussian Kernel Cumulative Distribution Function
    :param returns:           returns [1D array]
    :param x:                 index [float]
    :param smoothing:         smoothing parameter [float]
    :return:                  cumulative distribution [1D array]
    """
    sum_ = 0
    for i in range(1, len(returns)):
        sum_ = sum_ + normal_cdf((x - returns[i]) / smoothing)
    return sum_ / len(returns)


def kernel_pdf(returns, x, smoothing):
    """
    Compute Gaussian Kernel Probability Distribution Function
    :param returns:           returns [1D array]
    :param x:                 index [float]
    :param smoothing:         smoothing parameter [float]
    :return:                  probability distribution [1D array]
    """
    sum_ = 0
    for i in range(1, len(returns)):
        sum_ = sum_ + normal_pdf((x - returns[i]) / smoothing)
    return sum_ / (len(returns) * smoothing)


def normal_cdf(x, mean=0, std_dev=1):
    """
    Retrieve value from the normal Cumulative Distribution Function
    :param x:                 index [float]
    :param mean:              center of the density distribution [float]
    :param std_dev:           standard deviation of the density distribution [float]
    :return:                  cumulative distribution [1D array]
    """
    t = x - mean
    y = 0.5 * math.erfc(-t / (std_dev * math.sqrt(2.0)))
    if y > 1.0:
        y = 1.0
    return y


def normal_pdf(x, mean=0, std_dev=1):
    """
    Retrieve value from the normal Probability Distribution Function
    :param x:                 index [float]
    :param mean:              center of the density distribution [float]
    :param std_dev:           standard deviation of the density distribution [float]
    :return:                  probability distribution [1D array]
    """
    t = (x - mean) / abs(std_dev)
    y = math.exp(-(t ** 2) / 2) / math.sqrt(2 * math.pi * std_dev ** 2)
    return y


def normal_pdf_2D(x, y):
    """
    Retrieve value from the normal Probability Distribution Function
    :param x:                 first dimension indexes [float]
    :param y:                 second dimension indexes [float]
    :return:                  probability distribution [1D array]
    """
    y = np.exp(-(x ** 2 + y ** 2) / 2) / np.sqrt(2 * math.pi)
    return y


def compute_VaR(density_x, density_y, alpha):
    """
    Compute VaR from density distribution
    :param density_x:       density indexes [1D array]
    :param density_y:       density values [1D array]
    :param alpha:           confidence level [%]
    :return:                VaR [%]
    """
    sum_ = 0
    for i in range(0, len(density_x)):
        sum_ = sum_ + density_y[i]
        if sum_ >= alpha:
            return round(density_x[i] * 100, 3)


def compute_VaR_2D(density_x, density_y, density_z, last_return, alpha):
    """
    Compute VaR from density distribution
    :param density_x:       first dimension density indexes [1D array]
    :param density_y:       second dimension density indexes [1D array]
    :param density_z:       density values [1D array]
    :param last_return:     lagged return [float]
    :param alpha:           confidence level [%]
    :return:                VaR [%]
    """
    index = 0
    for i in range(0, len(density_y)):
        if density_y[i] >= last_return:
            index = i
            break
    density = density_z[:, index] / sum(density_z[:, index])
    sum_ = 0
    for i in range(0, len(density)):
        sum_ = sum_ + density[i]
        if sum_ >= alpha:
            return round(density_x[i] * 100, 3)


def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs):
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins,
             y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    zz = np.reshape(z, xx.shape)
    z = 100 * zz / (len(x) * len(y))

    return xx, yy, z


def update_axis_arrays(xx, yy, type="kernel"):
    if type == "normal":
        x = xx[0]
        y = np.empty(len(yy))
        for i in range(0, len(yy)):
            y[i] = yy[i][0]
    else:
        y = yy[0]
        x = np.empty(len(xx))
        for i in range(0, len(xx)):
            x[i] = xx[i][0]
    return x, y
