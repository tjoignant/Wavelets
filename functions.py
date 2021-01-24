import datetime as dt
import yfinance as yf
import numpy as np
import math


def get_data(ticker, interval):
    """
    Download ticker's ohlcv data for a chosen time period and time interval (either "1d", "1h", "30m", "15m", "5m")
    :param ticker       ticker id [str]
    :param interval     interval type [str]
    :return ohlcv       open high low close volume dataframe [df]
    """
    # Display indication
    print('[INFO] {} - Retrieving {}_{} historical data'.format(get_now(), ticker, interval))
    # Download ticker's ohlcv
    end_date = dt.datetime.now()
    if interval == '1d':
        start_date = end_date - dt.timedelta(days=365)
    elif interval == '1h':
        start_date = end_date - dt.timedelta(days=100)
    else:
        start_date = end_date - dt.timedelta(days=30)
    ohlcv = yf.download(tickers=ticker, start=start_date, end=end_date, interval=interval)
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
    :param ax:              subplot's axis [matplotlib axis]
    :param title:           subplot's title [str]
    :param xlabel:          subplot's horizontal axis's label [str]
    :param ylabel:          subplot's vertical axis's label [str]
    """
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

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
            returns[i] = math.log(prices[i]) - math.log(prices[i-1])
        elif method == "frac":
            returns[i] = (prices[i] - prices[i-1]) / prices[i-1]
        else:
            return "Method undefined"
    returns[0] = None
    return returns[1:len(returns)]

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
        for j in range(0, len(hist_array_x)-1):
            if hist_array_x[j] < returns[i] <= hist_array_x[j + 1]:
                hist_array_y[j] = hist_array_y[j] + 1
                break
        if returns[i] > hist_array_x[len(hist_array_x)-1]:
            hist_array_y[len(hist_array_x)-1] = hist_array_y[len(hist_array_x)-1] + 1
    hist_array_y = hist_array_y / sum(hist_array_y)
    return [hist_array_x, hist_array_y]

def compute_kernel_density(returns, smoothing):
    """
    Compute return's histogram density
    :param returns:                     returns array [1D array]
    :param smoothing:                   smoothing parameter [float]
    :return: density index & values     returns density [2D array]
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
        sum_ = sum_ + normal_cdf((x-returns[i])/smoothing)
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
        sum_ = sum_ + normal_pdf((x-returns[i])/smoothing)
    return sum_ / (len(returns)*smoothing)

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
    t = (x-mean)/abs(std_dev)
    y = math.exp(-(t**2)/2) / math.sqrt(2 * math.pi * std_dev ** 2)
    return y

def compute_VaR(density_x, density_y, alpha):
    """
    Compute VaR from density distribution
    :param density_x:       density index [1D array]
    :param density_y:       density value [1D array]
    :param alpha:           confidence level [%]
    :return:                VaR [%]
    """
    sum_ = 0
    for i in range(0, len(density_x)):
        sum_ = sum_ + density_y[i]
        if sum_ >= alpha:
            return round(density_x[i]*100, 3)
