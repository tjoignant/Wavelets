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

#### Father Wavelet
<p align="center">
  <img src=
  "https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A+++++%5CPhi%28x%29+%3D+%5Csum_%7Bi%3D0%7D%5E%7B2N-1%7D+%5Calpha_i+%5CPhi%282x-i%29%0A%5Cend%7Balign%2A%7D%0A" 
  alt="\begin{align*}
       \Phi(x) = \sum_{i=0}^{2N-1} \alpha_i \Phi(2x-i)
  \end{align*}
  ">
</p>

Tip: To construct the father wavelet we use a cascade algorithm

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

Tip: To reduce computing time one can take 

<br/>

### Mexican Hat

Not implemented yet ...



#### Discrete Wavelet Transformation (dwt)

#### Inverse Discrete Transformation (idwt)
