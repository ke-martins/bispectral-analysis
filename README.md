# Bispectral-Analysis
A MATLAB library designed for spectral and bispectral analyses of surface elevation timeseries.

Author: Kévin Martins

Affiliation: UMR 5805 EPOC, CNRS - University of Bordeaux, Allée Geoffroy Saint-Hilaire, F-33615 Pessac, France

Email : kevin.martins@u-bordeaux.fr

Last updated on 1 July 2022

Included files:

    fun_compute_spectrum.m         --> Compute the power spectral density (psd) of a surface elevation timeseries

    fun_compute_bispectrum.m       --> Compute the power bispectrum of a surface elevation timeseries

    fun_compute_bispectrum_H2001.m --> Compute the power bispectrum of a surface elevation timeseries + bicoherence (Hagihira et al., 2001) and biphase

    fun_compute_bispectrum_H1982.m --> Compute the power bispectrum of a surface elevation timeseries + bicoherence (Hinich, 1982) and biphase

    fun_compute_edof.m             --> Compute the equivalent number of degrees of freedom of a psd

    fun_compute_Snl.m              --> Compute the source term for non-linear energy transfer (Herbers and Burton, 1997)
    
    fun_compute_krms.m             --> Compute the root-mean square wavenumber (Herbers et al., 2002)
    
# Notes
Each function is commented and referenced when found necessary or appropriate. Currently, no running examples are provided but the information provided in the files header should be sufficient. Feel free to contact me for any question regarding the use of these routines or parts of their implementation.

# Acknowledgements
This work was undertaken during my postdoctoral fellowship at the University of Bordeaux, France. Their financial support, through an International Postdoctoral Grant (Idex, nb. 1024R-5030), is therefore greatly acknowledged. Fruitful exchanges with Steve Elgar regarding the bispectral analysis were greatly appreciated. The development of this library was initiated for the manuscript: "Dispersive characteristics of non-linear waves propagating and breaking over a mildly sloping laboratory beach", published in Coastal Engineering:
Martins, K., Bonneton, P., & Michallet, H. (2021). Dispersive characteristics of non-linear waves propagating and breaking over a mildly sloping laboratory beach. Coastal Engineering, 167, 103917. 
Ongoing development is pursued, through the support of the European Union’s Horizon 2020 research and innovation program under the Marie Skłodowska-Curie Grant Agreement 887867 (lidBathy).

