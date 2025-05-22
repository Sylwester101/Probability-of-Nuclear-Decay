# Probability-of-Nuclear-Decay
Calculates the probability for actinide atoms to undergo fission with respect to competing neutron and gamma decays.

## Dependencies
Uses the following python libraries:
- numpy
- matplotlib.pyplot
- scipy.interpolate
- math

## Features

The program uses the tanh-sinh quadrature algorithm in accordance with the paper: https://www.davidhbailey.com//dhbpapers/dhb-tanh-sinh.pdf

However, it remains unoptimized lacking the multithreading and calculation shortcuts that the paper suggests, therefore it may be slower than the tanh-sinh quadrature in the library 
python library mpmath.



The probabilities of the neutron and fission decay are based on the unintegrated decay width formulae derived in the paper: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.78.054604

The gamma decay width is only based on the equation presented in the GEF manual, and will need to be changed in the near future for a more theoretically driven equation rather than the experimental based equations GEF uses for their simulations. This equation is seen in section 8.1.2 in the [GEF manual](https://www.khschmidts-nuclear-web.eu/Preprints/db-doc2014-1.pdf).




The core of the program is the monte carlo simulation that an arbitrary number of samples undergoes. Based on the probabilites calculated from the decay widths, the atom selectes one of three paths: fission, neutron emission, or radiation.
If the atom undegoes fission the loop ends for that sample and the fission map adds a count to the respective atomic weight when the atom underwent fission.
If the atom emits a neutron, the atom's atomic weight is reduced by one, and the excited energy of this new atom is the excited energy of the original atom minus the neutron binding energy and some kinetic energy of the neutron that is selected through a monte carlo simulation going through a numerically inverted kinetic energy culumative distribution function. The atom then continues the loop, until it undergoes fission at some point. 
If the atom radiates, at the current moment it loses 2 MeV from the excited energy, however this will change in the near future to create a more accurate simulation. The atom also continues the loop until it undergoes fission.


## License

MIT © 2025 Sylwester Janowski
