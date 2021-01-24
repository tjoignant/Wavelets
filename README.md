# Wavelets

## Introduction

*What is a wavelet ?*

Small wave with a defined start and end. A signal can be reconstructed as the sum of small dilated and translated wavelets.

<br/>

<p align="center">
  <img src="./Images/Wavelets_Examples.png" alt="Test"/>
  <br/>
  <b>Fig.1</b> Wavelets Examples
</p>

<br/>

*Why do we use wavelets ?* <br/>

Break down a signal into 2 components of various scales and frequencies. Financial applications:
* Denoising
* Singularity prediction
* Detection of frequency patterns
* ...

<br/>

## Building the Wavelets

### Daubechie

The daubechie wavelet's are orthognal with a compact support. Their are the more commonly used for the discrete wavelet transformation.

<br/>

#### Father Wavelet

<img src="https://render.githubusercontent.com/render/math?math=\phi_{j, k}(u_t} = -1">

#### Mother Wavelet

#### Derived Wavelets

<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A+++++%5Cphi_%7Bj%2Ck%7D%28u_k%29+%3D+2%5E%7Bj%2F2%7D%5CPhi%282%5Ej+u_t+-+k%29%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
     \phi_{j,k}(u_k) = 2^{j/2}\Phi(2^j u_t - k)
\end{align*}
">

### Mexican Hat

Not implemented yet ...



#### Discrete Wavelet Transformation (dwt)

#### Inverse Discrete Transformation (idwt)
