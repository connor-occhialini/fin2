# Comparison of the Convergence of Monte Carlo Numerical Integrators
Here we compute the electrostatic potential energy between two cubes separated by a distance r each with the same, yet non-uniform charge density.  This is done by computing a six dimensional integral based upon each of the 3 local coordinates of each cube.  To compute this integral, we first use a high-quality GSL numerical integrator, VEGAS, which uses monte carlo methods along with chi^2 analysis which converges to the true value faster than naive applications of the monte carlo technique.  For comparison, we devise a simple monte carlo integrator which uses the mean value theorem for integrals to numerically compute the 6-dimensional integral in question.  To estimate the error in our "handmade" integrator, we compute the integral value 16 times and take the standard deviation of these values. Finally, we use the dipole approximation for each cube to estimate the potential energy between the cubes.

These three methods are shown separately in the following image, with a base-2 logarithmic x-axis and base-10 logarithmic y-axis.  The energies shown are the absolute value, with all energies actually being negative (indicating an attraction between the two cubes).

![](energy.png?)

Relative timing between the VEGAS integrator and the handmade integrator (with error estimation) yield run times of ~300 and ~200s respectively.  Furthermore, the uncertainty for the VEGAS integrator is approximately an order of magnitude smaller (10^{-1}).  Since naive monte carlo methods converge like 1/sqrt{N}, where N is the total number of monte carlo steps, we thus expect to increase the computations (and thus the runtime) by a factor of 10^2 on our handmade integrator to achieve the same accuracy as VEGAS.  This would translate into a runtime of approximately 20,000s, or about 6 hours.

Due to the long runtime of the main program, the output is supplied in the form of 'data.temp' file.

This code is in part inspired by code supplied by Dr. Rozman for Physics 2200 course.

This project is a collaboration between myself and Filip Bergabo (UConn).
