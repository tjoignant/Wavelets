# Wavelets

## 1 - Introduction

#### What is a wavelet ?

Small wave with a defined start and end. A signal can be reconstructed as the sum of small dilated and translated wavelets.

<br/>

<p align="center">
  <img src="./Images/Wavelets_Examples.png" />
  <br/>
  <b>Fig.1</b> - Wavelets Examples
</p>

<br/>

#### Why do we use wavelets ? <br/>

Break down a signal into 2 components of various scales and frequencies. Financial applications:
* Denoising
* Singularity prediction
* Detection of frequency patterns
* ...

<br/>

## 2 - Building the Wavelets

### Daubechie

The daubechie wavelet's are orthognal with a compact support. Their are the more commonly used for the discrete wavelet transformation.

#### Father Wavelet
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A+++++%5CPhi%28x%29+%3D+%5Csum_%7Bi%3D0%7D%5E%7B2N-1%7D+%5Calpha_i+%5CPhi%282x-i%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
       \Phi(x) = \sum_{i=0}^{2N-1} \alpha_i \Phi(2x-i)
  \end{align*}
  ">
</p>

*Tip: To construct the father wavelet we use a cascade algorithm*

<br/>

#### Mother Wavelet
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A+++++%5CPsi%28x%29+%3D+%5Csum_%7Bi%3D0%7D%5E%7B2N-1%7D+%28-1%29%5Ei+%5Calpha_%7B2N-1-i%7D+%5CPhi%282x-i%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
       \Psi(x) = \sum_{i=0}^{2N-1} (-1)^i \alpha_{2N-1-i} \Phi(2x-i)
  \end{align*}
  ">
</p>

#### Derived Wavelets
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A+++++%5Cphi_%7Bj%2Ck%7D%28u_k%29+%3D+2%5E%7Bj%2F2%7D%5CPhi%282%5Ej+u_t+-+k%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
       \phi_{j,k}(u_k) = 2^{j/2}\Phi(2^j u_t - k)
  \end{align*}
  ">
</p>

<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A+++++%5Cpsi_%7Bj%2Ck%7D%28u_k%29+%3D+2%5E%7Bj%2F2%7D%5CPsi%282%5Ej+u_t+-+k%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
       \psi_{j,k}(u_k) = 2^{j/2}\Psi(2^j u_t - k)
  \end{align*}
  ">
</p>

*Tip: To reduce computing time one can take* 
<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A++++k+%5Cin+%5B-2%5Ej+T%3B2%5Ej+T%5D%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
    k \in [-2^j T;2^j T]
\end{align*}
">

<br/>

### Mexican Hat

Yet to be implemented ...


<br/>

## 3 - Decomposition/Recomposition

Using a discrete Wavelet (Daubechie for example) one can easily decompose and recompose perfectly the signal by doing the following steps:

<br/>

<p align="center">
  <img src="./Images/Decomposition_Recomposition.png" />
  <br/>
  <b>Fig.2</b> - Signal Decomposition/Recomposition
</p>

#### Discrete Wavelet Transform (dwt)
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AcD_%7Bj%2Ck%7D+%3D+%5Csum_%7Bt%3D1%7D%5ET+z%28u_t%29+%5Cpsi_%7Bj%2Ck%7D%28u_t%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  cD_{j,k} = \sum_{t=1}^T z(u_t) \psi_{j,k}(u_t)
  \end{align*}
  ">
</p>

<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AcA_%7Bj%2Ck%7D+%3D+%5Csum_%7Bt%3D1%7D%5ET+z%28u_t%29+%5Cphi_%7Bj%2Ck%7D%28u_t%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  cA_{j,k} = \sum_{t=1}^T z(u_t) \phi_{j,k}(u_t)
  \end{align*}
  ">
</p>


#### Inverse Discrete Wavelet Transform (idwt)
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AyD_j%28u_t%29+%3D+%5Csum_%7Bk%7DcD_%7Bj%2Ck%7D+%5Cpsi_%7Bj%2Ck%7D%28u_t%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  yD_j(u_t) = \sum_{k}cD_{j,k} \psi_{j,k}(u_t)
  \end{align*}
  ">
</p>

<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AyA_j%28u_t%29+%3D+%5Csum_%7Bk%7DcA_%7Bj%2Ck%7D+%5Cphi_%7Bj%2Ck%7D%28u_t%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  yA_j(u_t) = \sum_{k}cA_{j,k} \phi_{j,k}(u_t)
  \end{align*}
  ">
</p>


#### Signal Recomposition
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Chat+z+%28u_t%29+%3D+%5Csum_j+yD_j+%28u_t%29+%2B+yA_%7Bj0%7D+%28u_t%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  \hat z (u_t) = \sum_j yD_j (u_t) + yA_{j0} (u_t)
  \end{align*}
  ">
</p>


#### Examples

<p align="center">
  <img src="./Images/Haar_Decomposition_Recomposition.png" />
  <br/>
  <b>Fig.3</b> - Haar Decomposition/Recomposition
</p>

<br/>

<p align="center">
  <img src="./Images/Daubechie_Decomposition_Recomposition.png" />
  <br/>
  <b>Fig.4</b> - Daubechie_4 Decomposition/Recomposition
</p>

*Tip: Use main_signal.py to plot thoses graphs*


<br/>

## 4 - Signal Denoising

To denoise a signal one simply needs to denoise the details coefficients:

<p align="center">
  <img src="./Images/Denoised_Signal.png" />
  <br/>
  <b>Fig.5</b> - Signal Denoising
</p>


#### Thresholds

Their are multiples ways of computing the threshold to use. The most two famous are:
* Universal threshold (default)
* SURE threshold

<ins>Universal threshold</ins>
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Clambda_%7Buniversal%7D+%3D+%5Csqrt%7B2%2Alog%28N%29%7D%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  \lambda_{universal} = \sqrt{2*log(N)}
  \end{align*}
  ">
</p>


<ins>SURE threshold</ins>
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Clambda_%7BSURE%7D+%3D+argmin_%7B0%3C%5Clambda%3C%5Clambda_%7Buniversal%7D%7DSURE%28%5Clambda%2Cx%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  \lambda_{SURE} = argmin_{0<\lambda<\lambda_{universal}}SURE(\lambda,x)
  \end{align*}
">
  
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0ASURE%28%5Clambda%2Cx%29+%3D+s+%2B+%5Csum_%7Bi%3D1%7D%5Es+%5Bmin%28%5Cmid+X_i%5Cmid%2C+%5Clambda%29%5D%5E2-2card%28i%3A%5Cmid+X_i%5Cmid+%3C%5Clambda%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
  SURE(\lambda,x) = s + \sum_{i=1}^s [min(\mid X_i\mid, \lambda)]^2-2card(i:\mid X_i\mid <\lambda)
  \end{align*}
  ">
</p>


#### Soft/Hard Thresholding

Two types of thresholding exists:
* Soft - Every values above the treshold is shrinken by the threshold value (default)
* Hard - Every values below the treshold is shrinken to 0

<p align="center">
  <img src="./Images/Threshold_Types.png" />
  <br/>
  <b>Fig.6</b> - Threshold_Types
</p>


#### Examples

<p align="center">
  <img src="./Images/Denoising_Examples.png" />
  <br/>
  <b>Fig.7</b> - Daubechie_4 Denoising
</p>

*Tip: Use main_denoising.py to plot this graph*


