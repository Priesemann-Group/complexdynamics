#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <softplus.h>

/*------------------------  SOLVER PARAMETERS  -------------------------*/

const double T_max = 20000.0;			// End time of the simulation
										// T_max=720 for the measuring the dynamics of the transient regime as in the semitransparent lines in Fig.6c/d of the article
const double T_transient = 6000.0;		// Transient time being disregarded at the beginning of the simulation
										// T_transient=180 for the measuring the dynamics of the transient regime as in the semitransparent lines in Fig.6c/d of the article
const double delta_t = 0.5;				// Time step size of the solver

/*------------------------  SWEEPED PARAMETER RANGES  -------------------------*/
// Note: change other default model parameters below (lines 64 onwards)

const char paramA_name[] = "Mmax"; // Name is only used for the file output. To change which parameters are swiped change lines 64 onwards
const double paramA_min = 0.7;
const double paramA_max = 0.1;
const double paramA_res = 0.003;

const char paramB_name[] = "tau"; // Name is only used for the file output. To change which parameters are swiped change lines 64 onwards
const double paramB_min = 0.1;
const double paramB_max = 50.;
const double paramB_res = 0.499;

/*----------------------------  RUN  -----------------------------*/
int main(void){
	
	printf("Main starts\n");
	// Create output folder
	char dirpath[] = "../../data/simulations/param_swipe_2D/sp_Mmax_tau_s025_gamma01_nu100days_R5_g1";
	char tmp[120];
	sprintf(tmp, "mkdir -p %s", dirpath);
	system(tmp);

	// Create output files
	sprintf(tmp, "%s/maximaI.dat", dirpath);
	FILE *maximaIdat = fopen(tmp, "w");
	sprintf(tmp, "%s/minimaI.dat", dirpath);
	FILE *minimaIdat = fopen(tmp, "w");
	sprintf(tmp, "%s/maximacases.dat", dirpath);
	FILE *maximacasesdat = fopen(tmp, "w");
	sprintf(tmp, "%s/minimacases.dat", dirpath);
	FILE *minimacasesdat = fopen(tmp, "w");
	sprintf(tmp, "%s/summarydat.dat", dirpath);
	FILE *summarydat = fopen(tmp, "w");
	
	fprintf(maximaIdat, "%s,%s,%s,%s", paramA_name,paramB_name,"t","i\n"); //CHANGE HERE
	fprintf(minimaIdat, "%s,%s,%s,%s", paramA_name,paramB_name,"t","i\n"); //CHANGE HERE
	fprintf(maximacasesdat, "%s,%s,%s,%s", paramA_name,paramB_name,"t","i\n"); //CHANGE HERE
	fprintf(minimacasesdat, "%s,%s,%s,%s", paramA_name,paramB_name,"t","i\n"); //CHANGE HERE
	fprintf(summarydat, "%s,%s,%s,%s,%s", paramA_name,paramB_name,"W", "av_I", "av_cost_F\n"); //CHANGE HERE
	
	
	// Sweep parameters
	for (double paramA = paramA_min; paramA <paramA_max; paramA += paramA_res){
		printf("%s starts: %f\n",paramA_name, paramA); //CHANGE HERE
		
		for (double paramB = paramB_min; paramB <paramB_max; paramB += paramB_res){

			// Change default parameters
			struct Params my_params = default_params;
			my_params.s = 0.25; //CHANGE HERE
			my_params.tau = paramB;  //CHANGE HERE
			my_params.nu = 1./100.;  //CHANGE HERE
			my_params.beta = 0.5;    //CHANGE HERE
			//my_params.g = 1.;    //CHANGE HERE
			my_params.epsilon = 0.00025;    //CHANGE HERE
			my_params.Mmax = paramA;    //CHANGE HERE

			// Create the model
			gsl_odeiv2_system SYSTEM = {sirODE, NULL, dim, &my_params};

			gsl_odeiv2_step *STEPPER = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, dim);
			
			double *y = malloc(dim*sizeof(double));
			double *yerr = malloc(dim*sizeof(double));
			
			// Comment in for pre-defined initial conditions
			y[0] = 1.-0.01;
			y[1] = 0.01;
			y[2] = 0.0;
			y[3] = 0.005;
			y[4] = 0.005;
			// Comment in for random initial conditions
			/*float S0 = (float)rand()/(float)(RAND_MAX/0.5);
			float R0 = (float)rand()/(float)(RAND_MAX/0.5);
			float I0 = 1-S0-R0;
		
			
			y[0] = S0;
			y[1] = I0;
			y[2] = R0;
			y[3] = I0;
			y[4] = I0;*/

			// Prepare observable variables
			double t = 0;
			
			double array_I[2] = {0,0};
			double array_Idot[2] = {0,0};
			double array_cases[2] = {0,0};
            
			double mean_I = 0;
			double F_cost = 0;

			double npeaks = 0;
			double t_firstpeak = 0;
			double t_lastpeak = 0;
			double foi_val = 0;
			double cases = 0;
            
			// Solve the model, print minima and maxima when observed
			while (t<T_max){
				
				array_I[0] = array_I[1];
				array_Idot[0] = array_Idot[1];
				
				gsl_odeiv2_step_apply(STEPPER, t, delta_t, y, yerr, NULL, NULL, &SYSTEM);
				

				foi_val = foi(y[4],y[1],&my_params);
				array_I[1] = y[1];
				array_Idot[1] = y[0]*foi_val*seasonality(t,&my_params) - my_params.gamma*y[1];
                
                // Detect Min/Max in I compartment
				if (array_Idot [0] > 0 && array_Idot[1] < 0 && t>T_transient){
					fprintf(maximaIdat, "%g,%g,%g,%g\n", paramA, paramB, t, y[1]);
					if (npeaks == 0) t_firstpeak = t;
					t_lastpeak = t;
					
					npeaks = npeaks+1;
				}
				if (array_Idot [0] < 0 && array_Idot[1] > 0 && t>T_transient){
					fprintf(minimaIdat, "%g,%g,%g,%g\n", paramA, paramB, t, y[1]);
				}
				
				if (t>T_transient) {
				mean_I += y[1];
				F_cost += 1/(foi_val/(y[1]*my_params.beta));
                }


                // Detect Min/Max in daily new cases
				cases = array_Idot[1] + my_params.gamma*y[1];
				
				if (array_cases[0] < array_cases[1] && array_cases[1] >= cases && t>T_transient){
					fprintf(maximacasesdat, "%g,%g,%g,%g\n", paramA, paramB, t, array_cases[1]);
				}
				if (array_cases[0] > array_cases[1] && array_cases[1] <= cases  && t>T_transient){
					fprintf(minimacasesdat, "%g,%g,%g,%g\n", paramA, paramB, t, array_cases[1]);
				}
				
				array_cases[0] = array_cases[1];
				array_cases[1] = cases;

				t += delta_t;
			}
			
			double const timediff = t_lastpeak-t_firstpeak;
			double const winding = (npeaks-1)*360/timediff;
			//printf("%g\n", winding);
			//printf("Windingnumber %g\n", winding);
			mean_I = mean_I/(T_max/delta_t);
			F_cost = F_cost/(T_max/delta_t);
			
			// Print summary output per simulation
			fprintf(summarydat, "%g,%g,%g,%g,%g\n", paramA, paramB, winding, mean_I, F_cost);
            
			free(y);
			free(yerr);
		}
		
	}
	
	printf("Run finished\n");
	
	fclose(maximaIdat);
	fclose(minimaIdat);
	fclose(maximacasesdat);
	fclose(minimacasesdat);
    fclose(summarydat);
    return EXIT_SUCCESS;
	
}
