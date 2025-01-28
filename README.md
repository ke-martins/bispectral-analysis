# bispectral-analysis: a MATLAB toolbox for ocean wave bispectral analysis 🌊

Welcome to this short documentation of bispectral-analysis, a MATLAB library designed for spectral and bispectral analyses of free surface elevation timeseries associated with ocean waves. Initially created for computing the bispectrum of the free surface elevation, it has now evolved to incorporate other relevant functions, broadening the range of applications (wave dispersive properties, non-linear energy transfers between triads etc). A list of relevant publications that used this toolbox is given at the bottom of this page. Although this toolbox was intended for surface elevation datasets, it can be used for other signals bearing in mind the provided units would be wrong.

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
bispectral-analysis v2: new example provided; updates on edof computation; remove (unphysical) merging option; fix bugs in function computing krms terms.  
<sub><sup>:arrow_forward:</sup></sub> *(May. 2021)*
bispectral-analysis v1: first release of the library, and tag version; stable and widely tested (example to be added).

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

## List of functions

<details>
  <summary>📄 fun_compute_spectrum.m</summary>  
  <br>  

  **Description**:  
  Direct FFT-based estimation of surface elevation spectral densities [m²/Hz].  

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `x`       | double | Detrended free surface elevation signal [m] |  
  | `fs`      | int    | Sampling frequency [Hz]                                         |  
  | `nfft`    | int    | Block length for the FFT (default = 256)                        |
  | `overlap` | int    | Percentage overlap (typical is 50%, 75% optimises edof)         |
  | `wind`    | char   | Type of window for tapering ('rectangular', 'hann', or 'kaiser') |

  **Outputs**:  
  &nbsp;&nbsp;Returns `data`, a self-explanatory MATLAB data structure containing spectral products.

</details>

<details>
  <summary>📄 fun_compute_bispectrum.m</summary>  
  <br>  

  **Description**:  
  Direct FFT-based estimation of the surface elevation bispectrum [m^3].  
  We here use the definition by Kim and Powers (1986):  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$`B(f_1,f_2) = E\left[ A(f_1) A(f_2) A^*(f_1 + f_2) \right]`$,  
  where $`A`$ are the complex Fourier coefficients, $`A^*`$ denotes the conjugate of $`A`$ and $`E`$ is the expected value. In this function, the bicoherence is not computed, as we wish to keep a light version, and it is not always needed. If needed, please use e.g. fun_compute_bispectrum_H2001. The option to frequency-merge bispectral estimates is no longer optional as we consider it not appropriate, since it leads to unrealistic bicoherence >1.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `x`       | double | Detrended free surface elevation signal [m] |  
  | `fs`      | int    | Sampling frequency [Hz]                                         |  
  | `nfft`    | int    | Block length for the FFT (default = 256)                        |
  | `overlap` | int    | Percentage overlap (typical is 50%, 75% optimises edof)                     |
  | `wind`    | char   | Type of window for tapering ('rectangular', 'hann', or 'kaiser') |

  **Outputs**:  
  &nbsp;&nbsp;Returns `data`, a self-explanatory MATLAB data structure containing spectral and bispectral products.

</details>

<details>
  <summary>📄 fun_compute_bispectrum_H1982.m</summary>  
  <br>  

  **Description**:  
  Direct FFT-based estimation of the surface elevation bispectrum [m^3].  
  We here use the definition by Kim and Powers (1979):  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$`B(f_1,f_2) = E\left[ A(f_1) A(f_2) A^*(f_1 + f_2) \right]`$,  
  where $`A`$ are the complex Fourier coefficients, $`A^*`$ denotes the conjugate of $`A`$ and $`E`$ is the expected value. The normalisation used here for the bicoherence $`b`$ is slightly modified from Hinich (1982):  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$`b(f1,f2) = |B(f1,f2)| / \sqrt{E[|A(f1)|^2] E[|A(f2)|^2] E[|A(f1 + f2)|^2]}`$  
  Compared to the original version, the absolute value is applied. This is used by Haubrich (1965) and Herbers et al. (1994). The option to frequency-merge bispectral estimates is no longer optional as we consider it not appropriate, since it leads to unrealistic bicoherence >1.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `x`       | double | Detrended free surface elevation signal [m] |  
  | `fs`      | int    | Sampling frequency [Hz]                                         |  
  | `nfft`    | int    | Block length for the FFT (default = 256)                        |
  | `overlap` | int    | Percentage overlap (typical is 50%, 75% optimises edof)                     |
  | `wind`    | char   | Type of window for tapering ('rectangular', 'hann', or 'kaiser') |

  **Outputs**:  
  &nbsp;&nbsp;Returns `data`, a self-explanatory MATLAB data structure containing spectral and bispectral products.

</details>

<details>
  <summary>📄 fun_compute_bispectrum_H2001.m</summary>  

  **Description**:  
  Direct FFT-based estimation of the surface elevation bispectrum [m^3].  
  We here use the definition by Kim and Powers (1986):  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$`B(f_1,f_2) = E\left[ A(f_1) A(f_2) A^*(f_1 + f_2) \right]`$,  
  where $`A`$ are the complex Fourier coefficients, $`A^*`$ denotes the conjugate of $`A`$ and $`E`$ is the expected value. The normalisation used here for the bicoherence $`b`$ is taken from Hagihira et al. (2001):  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$`b(f_1,f_2) = |B(f_1,f_2)| / E\left[ \sqrt{|A(f_1)|^2 |A(f_2)|^2 |A(f_1 + f_2)|^2}\, \right]`$  
  This is the most correct normalisation as it equals one for fully coupled triads (all in phase). The option to frequency-merge bispectral estimates is no longer optional as we consider it not appropriate, since it leads to unrealistic bicoherence >1.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `x`       | double | Detrended free surface elevation signal [m] |  
  | `fs`      | int    | Sampling frequency [Hz]                                         |  
  | `nfft`    | int    | Block length for the FFT (default = 256)                        |
  | `overlap` | int    | Percentage overlap (typical is 50%, 75% optimises edof)                     |
  | `wind`    | char   | Type of window for tapering ('rectangular', 'hann', or 'kaiser') |

  **Outputs**:  
  &nbsp;&nbsp;Returns `data`, a self-explanatory MATLAB data structure containing spectral and bispectral products.

</details>

<details>
  <summary>📄 fun_compute_krms.m</summary>  
  <br>  

  **Description**:  
  Compute the root-mean square wavenumker $`\kappa_{rms}`$ following Herbers et al. (2000, Eq. 12), based on the Boussinesq theory of Herbers and Burton (1997). In this work, the authors neglect fourth-order frequency terms, which can improve the linear predictions in deeper water depths. This is an option here.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `h0`      | double | local water depth [m] |  
  | `f`       | double | frequency array [Hz]                                            |  
  | `P`       | double | power spectrum [m^2], not a density                       |
  | `B`       | complex | power bispectrum [m^3], not a density       |
  | `option`  | char   | optional, 'second' or 'fourth' order frequency dispersion term |

  **Outputs**:  
  &nbsp;&nbsp;Returns `krms`, the root-mean square wavenumker k [1/m], whose size is that of input `f`.

</details>

<details>
  <summary>📄 fun_compute_krms_terms.m</summary>  
  <br>  

  **Description**:  
  Compute the non-linear frequency and amplitude dispersion terms contributing to the root-mean square wavenumker $`\kappa_{rms}`$ following Herbers et al. (2000). Compared to their Eq. 12, we retain frequency terms at order $`\mu^2`$, in order to improve the linear dispersive properties in deeper water depth. This function was written for Boussinesq-based depth inversion applications described in Martins et al. (2023), see the list of reference for more details, especially on notations employed here.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `f`       | double | frequency array [Hz]                                            |  
  | `P`       | double | power spectrum [m^2], not a density                       |
  | `B`       | complex | power bispectrum [m^3], not a density       |
  | `fc`      | double | optional, cutoff frequency [Hz], in case the raw timeseries was noisy     |

  **Outputs**:  
  &nbsp;&nbsp;`gamma_fr` - frequency dispersion term, corresponding to the original term $`\beta_{fr}/h`$, size of input `f`.  
  &nbsp;&nbsp;`gamma_fr2`- frequency dispersion term, corresponding to the original term $`\beta_{fr,2}/h^2`$, size of input `f`.  
  &nbsp;&nbsp;`gamma_am` - amplitude dispersion term, corresponding to the original term $`h\beta_{am}`$, size of input `f`.

</details>

<details>
  <summary>📄 fun_compute_Snl.m</summary>  
  <br>  

  **Description**:  
  Compute $`S_{nl}`$, the non-linear term for the energy transfer between triads following Herbers and Burton (1997), see their Eq. 14.  

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `h0`      | double | local water depth [m] |  
  | `f`       | double | frequency array [Hz]                                            |  
  | `B`       | complex | power bispectrum [m^3], not a density       |

  **Outputs**:  
  &nbsp;&nbsp;Returns `Snl`, the source term for non-linear energy transfers between triads [m^2], whose size is that of input `f`.

</details>

<details>
  <summary>📄 fun_compute_edof.m</summary>  
  <br>  

  **Description**:  
  Computes the spectral estimate effective degrees of freedom, following Percival and Walden (1993, their Eq. 292b).  

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `w`       | double | Window (windowed timeseries of FFT length)
  | `Ns`      | int    | Block length for the FFT                                         |  
  | `N`       | int    | Total number of points                        |
  | `overlap` | int    | Percentage overlap         |

  **Outputs**:  
  &nbsp;&nbsp;Returns `v`, effective degrees of freedom in a PSD estimate, following Percival and Walden (1993).

</details>

# Acknowledgements  

This work was undertaken in 2019 during my postdoctoral fellowship at the University of Bordeaux, France, while working on the manuscript: "Dispersive characteristics of non-linear waves propagating and breaking over a mildly sloping laboratory beach", published in Coastal Engineering in 2021 (see list of references). The financial support from the University, through an International Postdoctoral Grant (Idex, nb. 1024R-5030), is therefore greatly acknowledged. Fruitful exchanges at the time with Steve Elgar regarding the bispectral analysis were also appreciated. Since this work, continuous development has been pursued for various applications, including for the reconstruction of non-linear wave fields from pressure transducers, or for nearshore depth inversion based on Boussinesq theory, through the support of the European Union’s Horizon 2020 research and innovation program under the Marie Skłodowska-Curie Grant Agreement 887867 (lidBathy). In the **References** section below, you will find a list of published work where we heavily relied on this toolbox.


# References

 - Martins, K., Bonneton, P., & Michallet, H. (2021). Dispersive characteristics of non-linear waves propagating and breaking over a mildly sloping laboratory beach. *Coastal Engineering* <strong>167</strong>, 103917. https://doi.org/10.1016/j.coastaleng.2021.103917
 
 - Martins, K., Bonneton, P., Lannes, D., & Michallet, H. (2021). Relation between orbital velocities, pressure, and surface elevation in nonlinear nearshore water waves. *Journal of Physical Oceanography* <strong>51</strong>(11), 3539-3556. https://doi.org/10.1175/JPO-D-21-0061.1
 
 - Martins, K., Bonneton, P., de Viron, O., Turner, I. L., Harley, M. D., & Splinter, K. (2023). New Perspectives for Nonlinear Depth‐Inversion of the Nearshore Using Boussinesq Theory. *Geophysical Research Letters* <strong>50</strong>(2), e2022GL100498. https://doi.org/10.1029/2022GL100498
 
 - Sous, D., Martins, K., Tissier, M., Bouchette, F., & Meulé, S. (2023). Spectral wave dissipation over a roughness‐varying barrier reef. *Geophysical Research Letters* <strong>50</strong>(5), e2022GL102104. https://doi.org/10.1029/2022GL102104
 
 
 
 
 
 
