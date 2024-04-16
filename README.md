This repository contains all the codes needed to reproduce the results presented in our recent paper- "Rashid, H. U., Olorode, O. Use of controlled fractures in enhanced geothermal systems. Advances in Geo-Energy Research, 2024, 12(1): 35-51. https://doi.org/10.46690/ager.2024.04.04". This paper is available via open access.

The simulations performed using these codes are based on the "geothermal" and "upr" modules in the MATLAB Reservoir Simulation Toolbox (MRST), which is freely available at [mrst.no](https://www.sintef.no/projectweb/mrst/).

Please note that the codes were based on MRST 2021b. The major change from MRST 2022a+ is the use of the TestCase() function instead of the MRSTExample() function. We haven't tested the idea of simply updating the code to use TestCase() instead, but this might be something you might want to try out. 

To run the cases, you need to replace the lineSites2D.m function in modules\upr\pebi2D\lineSites2D.m with lineSites2D.m available in this repository. 

Please submit an issue if you have any problems running this code.
