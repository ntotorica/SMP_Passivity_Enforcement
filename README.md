## Author: Nathan Totorica
## Date: 5/14/2021

# Singularity Matrix Pertubation (SMP)
This code was written for a class project in the course entitled ECE 504: "ST: Passive Electromagnetic Systems" taught in Spring 2021 by Dr. Ata Zadehgol at the University of Idaho in Moscow. 

This code was developed, in part, based on the code developed by Jennifer Houle in ECE 504 "ST: Modern Circuit Synthesis Algorithms" taught in Spring 2020 by Dr. Ata Zadehgol. Jennifer's code is available online at https://github.com/JenniferEHoule/Circuit_Synthesis.

## Overview
Singular Matrix Perturbation (SMP) as introduced in [11], is a passivity enforcement algorithm for use on fitted models. This robust method is fast, computationally inexpensive, and accurate in enforcing passivity. This implementation using Python can easily be used with Vector Fitting algorithm implemented in [12]. Example cases to demonstrate passivity enforcmment were taken from [7][15], as well as custom examples designed in simulation software based off of multiport circuit synthesis as described in [6].

## I/O and Flow Description
 - Instructions for input and output and flow description can be found in flow_diagram.pdf. 
 
## Files
- SMP.py: Singularity matrix perturbation implementation of [11]
- s_compare.py: Compare ADS simulated S matrices. Calculate RMS error and plot
- eig_plot.py: Plotting functions. Eigenvalue plots based off of plot.py in [12]
- smp_ex.py: Example of how to pull .sp file in, parse parameters, send to vector fitting, enforce passivity, and generate new passive circuit.

Imported from [12]
- Ex2_y.py 
- vectfit3.py
- intercheig.py
- rot.py
- pass_check.py
- fitcalc.py
- FRPY.py
- intercheig.py
- violextrema.py
- plots.py
- quadprog.py
- pr2ss.py
- utils.py

## Licensing
GPL-3.0 License 

In addition to licensing:

        Embedding any of (or parts from) the routines of the Matrix Fitting Toolbox in a commercial software, or a software requiring licensing, is strictly prohibited. This applies to all routines, see Section 2.1.
        If the code is used in a scientific work, then reference should me made as follows:
        VFdriver.m and/or vectfit3.m: References [1],[2],[3]
        RPdriver.m and/or FRPY.m applied to Y-parameters: [8],[9]

## References
```
[1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
    domain responses by Vector Fitting", IEEE Trans. Power Delivery,
    vol. 14, no. 3, pp. 1052-1061, July 1999.

[2] B. Gustavsen, "Improving the pole relocating properties of vector
    fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
    July 2006.

[3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
    "Macromodeling of Multiport Systems Using a Fast Implementation of
    the Vector Fitting Method", IEEE Microwave and Wireless Components
    Letters, vol. 18, no. 6, pp. 383-385, June 2008.

[4] B. Gustavsen, VFIT3, The Vector Fitting Website. March 20, 2013. Accessed on:
    Jan. 21, 2020. [Online]. Available: 
    https://www.sintef.no/projectweb/vectfit/downloads/vfit3/.

[5] A. Zadehgol, "A semi-analytic and cellular approach to rational system characterization 
    through equivalent circuits", Wiley IJNM, 2015. [Online]. https://doi.org/10.1002/jnm.2119

[6] V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit Synthesis from 
    Frequency Response of Multi-port Networks", EMC EUR, pp. 79-84, 2016. [Online]. 
    Available: ://WOS:000392194100012.

[7] B. Gustavsen, Matrix Fitting Toolbox, The Vector Fitting Website.
    March 20, 2013. Accessed on: Feb. 25, 2020. [Online]. Available:
    https://www.sintef.no/projectweb/vectorfitting/downloads/matrix-fitting-toolbox/.

[8] B. Gustavsen, "Fast passivity enforcement for S-parameter models by perturbation
    of residue matrix eigenvalues",
    IEEE Trans. Advanced Packaging, vol. 33, no. 1, pp. 257-265, Feb. 2010.

[9] B. Gustavsen, "Fast Passivity Enforcement for Pole-Residue Models by Perturbation
    of Residue Matrix Eigenvalues", IEEE Trans. Power Delivery, vol. 23, no. 4,
    pp. 2278-2285, Oct. 2008.

[10] A. Semlyen, B. Gustavsen, "A Half-Size Singularity Test Matrix for Fast and Reliable
    Passivity Assessment of Rational Models," IEEE Trans. Power Delivery, vol. 24, no. 1,
    pp. 345-351, Jan. 2009.

[11] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions
     on Power Delivery, vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

[12] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: February 3, 2021. [Online]. Available: https://github.com/jenniferEhoule/circuit_synthesis
```
