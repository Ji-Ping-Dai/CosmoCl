# CosmoCl
This package is based on [CosmoMC](https://github.com/cmbant/CosmoMC) (version: Jun2016), a Fortran 2008 parallelized MCMC sampler. It can be used to calculate the 2-D auto-correlation power spectrum of large scale structure (LSS) tracers and the cross-correlation with CMB. Thus we can constrain  the cosmological parameters such as the primordial non-Gaussianity [(J. P. Dai, J. Q. Xia)](https://academic.oup.com/mnrasl/article/491/1/L61/5614501) and ISW effect [(X. H. Tan, J. P. Dai, J. Q. Xia)](https://doi.org/10.1016/j.dark.2020.100585).<br>

## How to use
You can install or run our package following CosmoMC. Here we provide a [guide](https://arxiv.org/pdf/1808.05080).<br>
The main part to calculate the correlation power spectrum is in **source/mpk_cl.f90** .<br>

## How to extend
In this package, it does not consider any complicated theories. Actually, we only calculate the auto-correlation of LSS number density perturbations and its cross-correlation with CMB temperature fluctuations. The only free parameter is the bias term. The corresponding formulae are:<br>
![](https://github.com/Ji-Ping-Dai/CosmoCl/blob/master/docs/Cl1.PNG)<br>
![](https://github.com/Ji-Ping-Dai/CosmoCl/blob/master/docs/Cl2.PNG)<br>
You can extend this package to other measurements such as weak lensing by changing the transfer functions and window functions.
