References
==========

Model Development
~~~~~~~~~~~~~~~~~

The MNiShed model is described on the Community Surface Dynamics Modeling
System (CSDMS) model page:

* `CSDMS Model Page: MNiShed <https://csdms.colorado.edu/wiki/Model:MNiShed>`_

Snowpack and Rain-on-Snow
~~~~~~~~~~~~~~~~~~~~~~~~~

Hock, R. (2003). Temperature index melt modelling in mountain areas.
*Journal of Hydrology*, 282(1-4), 104--115.
https://doi.org/10.1016/S0022-1694(03)00257-9

McCabe, G. J., Hay, L. E., & Clark, M. P. (2007). Rain-on-snow events in the
western United States.
*Bulletin of the American Meteorological Society*, 88(3), 319--328.
https://doi.org/10.1175/BAMS-88-3-319

Würzer, S., Jonas, T., Wever, N., & Lehning, M. (2016). Influence of initial
snowpack properties on runoff formation during rain-on-snow events.
*Journal of Hydrometeorology*, 17(6), 1801--1815.
https://doi.org/10.1175/JHM-D-15-0181.1

Frozen Ground
~~~~~~~~~~~~~

Dunne, T., & Black, R. D. (1971). Runoff processes during snowmelt.
*Water Resources Research*, 7(5), 1160--1172.
https://doi.org/10.1029/WR007i005p01160

Molnau, M., & Bissell, V. C. (1983). A continuous frozen ground index for
flood forecasting. *Proceedings of the 51st Annual Western Snow Conference*,
pp. 109--119.
https://westernsnowconference.org/bibliography/1983Molnau.pdf
(Original source for the continuous frozen ground index formulation,
including the exponential snow-depth insulation factor.)

Shanley, J. B., & Chalmers, A. (1999). The effect of frozen soil on snowmelt
runoff at Sleepers River, Vermont.
*Hydrological Processes*, 13(12--13), 1843--1857.
https://doi.org/10.1002/(SICI)1099-1085(199909)13:12/13<1843::AID-HYP879>3.0.CO;2-G

Evapotranspiration
~~~~~~~~~~~~~~~~~~

Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
*Crop Evapotranspiration: Guidelines for Computing Crop Water Requirements.*
FAO Irrigation and Drainage Paper 56. FAO, Rome.
https://www.fao.org/3/x0490e/x0490e00.htm
(Equations 2--4 used for photoperiod computation in the Thornthwaite--Chang method.)

Thornthwaite, C. W. (1948). An approach toward a rational classification of
climate. *Geographical Review*, 38(1), 55--94.
https://doi.org/10.2307/210739
(Original temperature-index ET method; the Chang (2019) variant used here is a
daily modification of it.)

Camargo, A. P., Marin, F. R., Sentelhas, P. C., & Picini, A. G. (1999).
Adjust of the Thornthwaite's method to estimate the potential
evapotranspiration for arid and superhumid climates based on daily temperature
amplitude.
*Revista Brasileira de Agrometeorologia*, 7(2), 251--257.
(Source for the k = 0.72 coefficient for monthly ET₀.)

Chang, X., Wang, S., Gao, Z., Luo, Y., & Chen, H. (2019). Forecast of daily
reference evapotranspiration using a modified daily Thornthwaite equation and
temperature forecasts.
*Irrigation and Drainage*, 68(2), 297--317.
https://doi.org/10.1002/ird.2309

Pereira, A. R., & Pruitt, W. O. (2004). Adaptation of the Thom and Oliver
evaporation equation to daily totals.
*Agricultural and Forest Meteorology*, 124(1-2), 26--32.
https://doi.org/10.1016/j.agrformet.2004.01.005

Schwartz, M. D., Ault, T. R., & Betancourt, J. L. (2013). Spring onset
variations and trends in the continental United States: past and regional
assessment using temperature-based indices.
*International Journal of Climatology*, 33(13), 2917--2922.
https://doi.org/10.1002/joc.3625
(Temperature-based spring-index phenology, background for the growing-degree-day
vegetation coefficient.)

Reservoir Cascade Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bergström, S. (1976). *Development and Application of a Conceptual Runoff
Model for Scandinavian Catchments.* SMHI Report RHO 7. Swedish
Meteorological and Hydrological Institute, Norrköping.
(Original HBV model; multi-component runoff cascade structure.)

Moore, R. J. (1985). The probability-distributed principle and runoff
production at point and basin scales.
*Hydrological Sciences Journal*, 30(2), 273--297.
https://doi.org/10.1080/02626668509490989
(Source for the PDM saturation-excess runoff option,
:math:`f_\text{sat} = 1 - e^{-H/H_0}`.)

Recession Analysis
~~~~~~~~~~~~~~~~~~

Brutsaert, W., & Nieber, J. L. (1977). Regionalized drought flow hydrographs
from a mature glaciated plateau.
*Water Resources Research*, 13(3), 637--643.
https://doi.org/10.1029/WR013i003p00637
(Basis for the :class:`~mnished.BrutsaertNieber` recession analysis and the
data-driven recession-exponent prior.)

Channel Routing
~~~~~~~~~~~~~~~

Nash, J. E. (1957). The form of the instantaneous unit hydrograph.
*IAHS Publication*, 45, 114--121.
(Original source for the N-reservoir Nash cascade and its gamma IUH.)

Dooge, J. C. I. (1959). A general theory of the unit hydrograph.
*Journal of Geophysical Research*, 64(2), 241--256.
https://doi.org/10.1029/JZ064i002p00241

Rodriguez-Iturbe, I., & Valdés, J. B. (1979). The geomorphologic structure
of hydrologic response.
*Water Resources Research*, 15(6), 1409--1420.
https://doi.org/10.1029/WR015i006p01409

Model Calibration and Evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eckhardt, K. (2005). How to construct recursive digital baseflow separation
filters: A discussion of two alternative filtering approaches.
*Hydrological Processes*, 19(2), 507--515.
https://doi.org/10.1002/hyp.5675

Yilmaz, K. K., Gupta, H. V., & Wagener, T. (2008). A process-based
diagnostic approach to model evaluation: Application to the NWS
distributed hydrological model.
*Water Resources Research*, 44(9), W09417.
https://doi.org/10.1029/2007WR006716
(Basis for the ``KGE_logKGE`` composite metric: no single efficiency
measure captures both high- and low-flow behaviour.)

Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009).
Decomposition of the mean squared error and NSE: Implications for improving
hydrological modelling.
*Journal of Hydrology*, 377(1-2), 80--91.
https://doi.org/10.1016/j.jhydrol.2009.08.003

Nash, J. E., & Sutcliffe, J. V. (1970). River flow forecasting through
conceptual models part I -- A discussion of principles.
*Journal of Hydrology*, 10(3), 282--290.
https://doi.org/10.1016/0022-1694(70)90255-6

Brun, R., Reichert, P., & Künsch, H. R. (2001). Practical identifiability
analysis of large environmental simulation models.
*Water Resources Research*, 37(4), 1015--1030.
https://doi.org/10.1029/2000WR900350
(Identifiability via parameter sensitivities and collinearity — the basis for
the profile and curvature-eigenspectrum diagnostics.)

Gutenkunst, R. N., Waterfall, J. J., Casey, F. P., Brown, K. S., Myers, C. R.,
& Sethna, J. P. (2007). Universally sloppy parameter sensitivities in systems
biology models. *PLoS Computational Biology*, 3(10), e189.
https://doi.org/10.1371/journal.pcbi.0030189
(Origin of the stiff vs. sloppy eigenspectrum framing used in the
identifiability tools.)

Vrugt, J. A. (2016). Markov chain Monte Carlo simulation using the DREAM
software package: Theory, concepts, and MATLAB implementation.
*Environmental Modelling & Software*, 75, 273--316.
https://doi.org/10.1016/j.envsoft.2015.08.013
(The DREAM sampler used for the Bayesian / uncertainty-quantification path.)

Software
~~~~~~~~

Harris, C. R., et al. (2020). Array programming with NumPy.
*Nature*, 585, 357--362.
https://doi.org/10.1038/s41586-020-2649-2

McKinney, W. (2010). Data structures for statistical computing in Python.
In *Proceedings of the 9th Python in Science Conference* (Vol. 445, pp. 56--61).

Houska, T., Kraft, P., Chamorro-Chavez, A., & Breuer, L. (2015). SPOTting model
parameters using a ready-made Python package.
*PLoS ONE*, 10(12), e0145180.
https://doi.org/10.1371/journal.pone.0145180
(The SPOTPY package driving the in-process SCE-UA and DREAM calibration runners.)

External Resources
~~~~~~~~~~~~~~~~~~

* `CSDMS: Community Surface Dynamics Modeling System <https://csdms.colorado.edu/>`_
* `USGS Water Resources <https://www.usgs.gov/mission-areas/water-resources>`_
