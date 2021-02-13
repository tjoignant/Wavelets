import functions

import math
import numpy as np


class Daubechie:
    # CONSTRUCTOR (public)
    def __init__(self, signals, nb_vanishing_moments):
        """
        Initialize an class instance
        :param signals:                     signal to analyse [2D array]
        :param nb_vanishing_moments:        amount of vanishing moments to use for the wavelet, 1 --> Haar & 2,3,...,10 --> Daubechie [int]
        """
        self.signal = signals[0]
        self.signal_lagged = signals[1]
        self.len_signal = len(self.signal)
        self.nb_moments = nb_vanishing_moments
        self.scaling_function_matrix = self._build_father_wavelet()
        self.white_noise_std = 1

    # METHODS (public)
    def density(self, j, estimator="linear"):
        """
        Compute the signal's density
        :param j:               scaling parameter [int]
        :param estimator:       estimator type [str]
        :return:                density function estimation values [1D array]
        """
        [c, d] = self._scaling_coefficients(j)
        lowest_value = round(min(self.signal), 3)
        highest_value = round(max(self.signal), 3)
        density_x = np.linspace(lowest_value, highest_value, self.len_signal)
        density_y = np.linspace(lowest_value, highest_value, self.len_signal)
        density_z = np.empty(shape=(len(density_x), len(density_y)))
        print("[INFO] {} - Computing the signal's density using the wavelet {} estimator (resolution level: {})".format(
            functions.get_now(), estimator, j))
        if estimator == "donoho" or estimator == "thresholded":
            density_z = self._density_donoho(j, c, d, density_x, density_y)
        elif estimator == "linear":
            density_z = self._density_linear(j, c, density_x, density_y)
        else:
            print('\n[ERROR] {} - Undefined estimator)'.format(functions.get_now()))
            quit()
        density_x, density_y = np.meshgrid(density_x, density_y)
        return [density_x, density_y, density_z]

    def _density_linear(self, j, c_coef, x_array, y_array):
        """
        Compute the signal's density (Linear estimator)
        :param j:           scaling parameter [int]
        :param c_coef:      scaling coefficients C [1D array]
        :param x_array:     density index [1D array]
        :return:            density values [1D array]
        """
        k_lim_ = int(pow(2, -j) * self.len_signal)
        z_array = np.empty(shape=(len(x_array), len(y_array)))
        cpt = 0
        perc = 0
        for i in range(0, len(x_array)):
            for ii in range(0, len(y_array)):
                sum_ = 0
                for k in range(-k_lim_, k_lim_):
                    for k2 in range(-k_lim_, k_lim_):
                        sum_ = sum_ + c_coef[k + k_lim_, k2 + k_lim_] * self._derived_father_wavelet(x_array[i], j, k) * self._derived_father_wavelet(y_array[ii], j, k2)
                z_array[i, ii] = sum_
                if 100 * cpt / (len(x_array) * len(y_array)) > perc + 1:
                    perc = int(100 * cpt / (len(x_array) * len(y_array)))
                    print("{}%".format(perc))
                cpt = cpt + 1
        z_array[z_array < 0] = 0
        # z_array = z_array / sum(z_array)
        return z_array

    def _density_donoho(self, j, c_coef, d_coef, x_array, y_array):
        """
        Compute the signal's density (Linear estimator)
        :param j:           scaling parameter [int]
        :param c_coef:      scaling coefficients C [1D array]
        :param x_array:     density index [1D array]
        :return:            density values [1D array]
        """
        k_lim_ = int(pow(2, -j) * self.len_signal)
        z_array = np.empty(shape=(len(x_array), len(y_array)))
        cpt = 0
        perc = 0
        for i in range(0, len(x_array)):
            for ii in range(0, len(y_array)):
                sum_ = 0
                for k in range(-k_lim_, k_lim_):
                    for k2 in range(-k_lim_, k_lim_):
                        sum_ = sum_ + c_coef[k + k_lim_, k2 + k_lim_] * self._derived_father_wavelet(x_array[i], j, k) * self._derived_father_wavelet(y_array[ii], j, k2)
                for j_ in range(j, j+1):
                    for k_ in range(-k_lim_, k_lim_):
                        for q in range(0, 3):
                            sum_ = sum_ + d_coef[q, k_ + k_lim_] * self._derived_mother_wavelet_multidim(x_array[i], y_array[ii], j, k_, q)
                if 100 * cpt / (len(x_array) * len(y_array)) > perc + 1:
                    perc = int(100 * cpt / (len(x_array) * len(y_array)))
                    print("{}%".format(perc))
                cpt = cpt + 1
                z_array[i, ii] = sum_
        z_array[z_array < 0] = 0
        # z_array = z_array / sum(z_array)
        return z_array

    def _scaling_coefficients(self, j):
        """
        Compute the wavelet's scaling coefficients C & D
        :param j:       scaling parameter [int]
        :return:        scaling coefficients [2D array]
        """
        c = self._scaling_coef_c(j)
        d = self._scaling_coef_d(j)
        return [c, d]

    def _scaling_coef_c(self, j):
        """
        Compute the wavelet's scaling coefficients C
        :param j:       scaling parameter [int]
        :return:        scaling coefficients D [1D array]
        """
        k_lim_ = int(pow(2, -j) * self.len_signal)
        c = np.empty(shape=(2 * k_lim_, 2 * k_lim_))
        for k in range(-k_lim_, k_lim_):
            for k2 in range(-k_lim_, k_lim_):
                sum_ = 0
                for i in range(0, self.len_signal):
                    sum_ = sum_ + self._derived_father_wavelet(self.signal[i], j, k) * self._derived_father_wavelet(self.signal_lagged[i], j, k2)
                c[k + k_lim_, k2 + k_lim_] = sum_ / self.len_signal
        return c

    def _scaling_coef_d(self, j):
        """
        Compute the wavelet's scaling coefficients C
        :param j:       scaling parameter [int]
        :return:        scaling coefficients D [1D array]
        """
        k_lim_ = int(pow(2, -j) * self.len_signal)
        d = np.empty(shape=(3, 2 * k_lim_))
        for q in range(0, 3):
            for k in range(-k_lim_, k_lim_):
                sum_ = 0
                for i in range(0, self.len_signal):
                    sum_ = sum_ + self._derived_mother_wavelet_multidim(self.signal[i], self.signal_lagged[i], j, k, q)
                d[q, k + k_lim_] = sum_ / self.len_signal
        return d

    def _derived_mother_wavelet_multidim(self, t1, t2, j, k, q):
        if q == 0:
            return self._derived_father_wavelet(t1, j, k) * self._derived_mother_wavelet(t2, j, k)
        elif q == 1:
            return self._derived_mother_wavelet(t1, j, k) * self._derived_father_wavelet(t2, j, k)
        else:
            return self._derived_father_wavelet(t1, j, k) * self._derived_father_wavelet(t2, j, k)

    def _derived_mother_wavelet(self, t, j, k):
        """
        Retrieve derived mother wavelet
        :param t:       time index [float]
        :param j:       scaling parameter [int]
        :param k:       translation parameter [int]
        :return:        derived mother wavelet value [float]
        """
        x = pow(2, j) * t - k
        mu = pow(2, j / 2) * self._mother_wavelet(x)
        return mu

    def _derived_father_wavelet(self, t, j, k):
        """
        Retrieve derived father wavelet
        :param t:       time index [float]
        :param j:       scaling parameter [int]
        :param k:       translation parameter [int]
        :return:        derived father wavelet value [float]
        """
        x = pow(2, j) * t - k
        mu = pow(2, j / 2) * self._father_wavelet(x)
        return mu

    def _mother_wavelet(self, t):
        """
        Retrieve mother wavelet
        :param t:       time index [float]
        :return:        mother wavelet value [float]
        """
        sum_ = 0
        for k in range(0, 2 * self.nb_moments):
            sum_ += pow(-1, k) * coefficients_dict[self.nb_moments][2 * self.nb_moments - 1 - k] * \
                    self._father_wavelet(2 * t - k)
        return sum_

    def _father_wavelet(self, x):
        """
        Retrieve father wavelet
        :param x:       time index [float]
        :return:        father wavelet value [float]
        """
        # Initialization
        [t, phi_t] = self.scaling_function_matrix
        phi_ = 0
        # Only if x belongs to the wavelet's support
        if t[0] <= x <= t[len(t) - 1]:
            # Retrieve corresponding phi value
            for i in range(0, len(t) - 1):
                if x == t[i]:
                    phi_ = phi_t[i]
                    break
                elif x < t[i]:
                    phi_ = self._interpolate(x, t[i], t[i + 1], phi_t[i], phi_t[i + 1])
                    break
        return phi_

    def _build_father_wavelet(self):
        """
        Build the father wavelet
        :return:        father wavelet matrix [2D array]
        """
        # Display indication to user
        if self.nb_moments == 1:
            nb_iterations = 10
            wavelet_name = "haar"
        else:
            nb_iterations = 6
            wavelet_name = "daubechie_{}".format(self.nb_moments)
        print('\n[INFO] {} - Building the {} father wavelet ({} iterations)'.format(functions.get_now(), wavelet_name,
                                                                                    nb_iterations))
        # Define support
        if self.nb_moments == 1:
            support = [0, 2]
        else:
            support = [0, 20]
        # Initialization
        step = 1
        t = np.zeros(support[1] - support[0] + 1)
        phi_t = np.zeros(support[1] - support[0] + 1)
        for i in range(0, len(t)):
            t[i] = support[0] + i
            if t[i] == 1:
                phi_t[i] = 1
        # Cascade algorithm
        for j in range(2, nb_iterations + 1):
            step = step / 2
            t_tampon = t
            phi_t_tampon = phi_t
            t = np.zeros(2 * len(t_tampon) - 1)
            phi_t = np.zeros(2 * len(t_tampon) - 1)
            for i in range(0, len(t), 2):
                t[i] = t_tampon[int(i / 2)]
                phi_t[i] = phi_t_tampon[int(i / 2)]
            for i in range(1, len(t) - 1, 2):
                t[i] = (t_tampon[int(i / 2)] + t_tampon[int(i / 2) + 1]) / 2
                sum_ = 0
                for k in range(0, 2 * self.nb_moments):
                    sum_ = sum_ + coefficients_dict[self.nb_moments][k] * \
                           self._support_value(2 * t[i] - k, t_tampon, phi_t_tampon)
                phi_t[i] = sum_
        # Interpolate Daubechie Wavelet to smoothen the function
        if self.nb_moments != 1:
            for i in range(0, nb_iterations + 1):
                for int_index in range(0, support[1]):
                    if i == 0:
                        index = self._support_index(int_index, t)
                        phi_t[index] = (phi_t[index - 1] + phi_t[index + 1]) / 2
                    else:
                        for j in range(1, pow(2, i), 2):
                            index = self._support_index(int_index + j / pow(2, i), t)
                            phi_t[index] = (phi_t[index - 1] + phi_t[index + 1]) / 2
        # Find last used index
        check = False
        index = len(t) - 1
        for i in range(0, len(t)):
            if phi_t[i] == 0 and check is False:
                index = i
                check = True
            if phi_t[i] != 0 and check is True:
                check = False
        # Remove useless indexes
        t = t[0:index]
        phi_t = phi_t[0:index]
        return [t, phi_t]

    @staticmethod
    def _support_value(t, t_tampon, phi_t_tampon):
        """
        Retrieve value from tampon array
        :param t:               time index [int]
        :param t_tampon:        time indexes array [1D array]
        :param phi_t_tampon:    support values array [1D array]
        :return:                support value [float]
        """
        phi = 0
        for i in range(0, len(t_tampon)):
            if t == t_tampon[i]:
                phi = phi_t_tampon[i]
                break
        return phi

    @staticmethod
    def _support_index(time, t):
        """
        Retrieve time's index from father wavelet's support
        :param time:     time value [float]
        :param t:        time indexes array [1D array]
        :return:         time value index [int]
        """
        time_index = 0
        for i in range(0, len(t)):
            if t[i] == time:
                time_index = i
                break
        return time_index

    @staticmethod
    def _interpolate(x, x1, x2, y1, y2):
        """
        Interpolate between two points
        :param x:       index to interpolate [float]
        :param x1:      first point index [float]
        :param x2:      second point index [float]
        :param y1:      first point value [float]
        :param y2:      second point value [float]
        :return:        value interpolated [float]
        """
        m = (y1 - y2) / (x1 - x2)
        y = (x - x2) * m + y2
        return y


# DICTIONARY (private)
coefficients_dict = dict([
    # (nb_moments, [ak_values])
    (1, [1, 1]),
    (2, [0.6830127, 1.1830127, 0.317069873, -0.1830127]),
    (3, [0.47046721, 1.14111692, 0.650365, -0.19093442, -0.12083221, 0.0498175]),
    (4, [0.32580343, 1.01094572, 0.89220014, -0.03957503, -0.26450717, 0.0436163, 0.0465036, -0.01498699]),
    (5, [0.22641898, 0.85394354, 1.02432694, 0.19576696, -0.34265671, -0.04560113, 0.10970265, -0.00882680,
         -0.01779187,
         0.00471742793]),
    (6,
     [0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.31998660, -0.18351806, 0.13788809, 0.03892321, -0.04466375,
      0.000783251152, 0.00675606236, -0.00152353381]),
    (7,
     [0.11009943, 0.56079128, 1.03114849, 0.66437248, -0.20351382, -0.31683501, 0.1008467, 0.11400345, -0.05378245,
      -0.02343994, 0.01774979, 6.07514995 * 10 - 4, -2.54790472 * 10 - 3, 5.00226853 * 10 - 4]),
    (8, [0.07695562, 0.44246725, 0.95548615, 0.95548615, -0.02238574, -0.40165863, 6.68194092 * 10 - 4, 0.18207636,
         -0.02456390, -0.06235021, 0.01977216, 0.01236884, -6.88771926 * 10 - 3, -5.54004549 * 10 - 4,
         9.55229711 * 10 - 4,
         -1.66137261 * 10 - 4]),
    (9,
     [0.05385035, 0.34483430, 0.85534906, 0.92954571, 0.18836955, -0.41475176, -0.13695355, 0.21006834, 0.043452675,
      -0.09564726, 3.54892813 * 10 - 4, 0.03162417, -6.67962023 * 10 - 3, -6.05496058 * 10 - 3, 2.61296728 * 10 - 3,
      3.25814671 * 10 - 4,
      -3.56329759 * 10 - 4, 5.5645514 * 10 - 5]),
    (10,
     [0.03771716, 0.26612218, 0.74557507, 0.97362811, 0.39763774, -0.35333620, -0.27710988, 0.18012745, 0.13160299,
      -0.10096657, -0.04165925, 0.04696981, 5.10043697 * 10 - 3, -0.01517900, 1.97332536 * 10 - 3,
      2.81768659 * 10 - 3,
      -9.69947840 * 10 - 4, -1.64709006 * 10 - 4, -1.64709006 * 10 - 4, -1.875841 * 10 - 5]),
])
