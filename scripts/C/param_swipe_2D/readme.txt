A C script to calculate multiple observables for a two-parameter sweep. Requires gsl to solve the differential equation. The sweeped model parameters as well as the behavioral feedback function can be changed by slight changes in the code.

The script creates a folder including multiple files:
- maximaI.dat: listing all maxima of the I compartment per parameter set with timestamp
- minimaI.dat: listing all minima of the I compartment per parameter set with timestamp
- maximacases.dat: listing all maxima of the new cases (positive contribution to $\dot I$) per parameter set with timestamp
- minimacases.dat: listing all minima of the new cases (positive contribution to $\dot I$) per parameter set with timestamp
- summarydat.dat: listing summary observables per parameter set, i.e. 'winding' W (number of peaks per year), 'mean_I' (average I compartment), 'F_cost' (average mitigation cost given as $\langle 1-C_M(t)\rangle$)

By default the script creates a sweep in Mmax and tau, as in Fig. 4a of the article. To build the script type 'make sp' in the command line.

To change default parameters, adapt lines 64 onwards.

To change the model, i.e. the behavioral feedback function from the default softplus function to logistic or something else, change line 7 in param_swipe_2D.c from

#include <softplus.h>

e.g. to

#include <logistic.h>

, eventually adapt the necessary parameters (lines 64 onwards), and finally rebuild the script using 'make log'.