pwr: param_swipe_2D_pwr.o power_law.o
	gcc -Wall param_swipe_2D_pwr.o power_law.o -lgsl -lgslcblas -lm -o swipe2D_pwr

param_swipe_2D_pwr.o: param_swipe_2D.c ../models/power_law.h
	gcc -I../models/ -Wall -c param_swipe_2D.c -o param_swipe_2D_power.o
	
power_law.o: ../models/power_law.h ../models/power_law.c
	gcc -Wall -c ../models/power_law.c

	
exp: param_swipe_2D_exp.o exponential.o
	gcc -Wall param_swipe_2D_exp.o exponential.o -lgsl -lgslcblas -lm -o swipe2D_exp

param_swipe_2D_exp.o: param_swipe_2D.c ../models/exponential.h
	gcc -I../models/ -Wall -c param_swipe_2D.c -o param_swipe_2D_exp.o
	
exponential.o: ../models/exponential.h ../models/exponential.c
	gcc -Wall -c ../models/exponential.c

	
sp: param_swipe_2D_sp.o softplus.o
	gcc -Wall param_swipe_2D_sp.o softplus.o -lgsl -lgslcblas -lm -o swipe2D_sp

param_swipe_2D_sp.o: param_swipe_2D.c ../models/softplus.h
	gcc -I../models/ -Wall -c param_swipe_2D.c -o param_swipe_2D_sp.o
	
softplus.o: ../models/softplus.h ../models/softplus.c
	gcc -Wall -c ../models/softplus.c

	
spr: param_swipe_2D_spr.o softplus_reparam.o
	gcc -Wall param_swipe_2D_spr.o softplus_reparam.o -lgsl -lgslcblas -lm -o swipe2D_spr

param_swipe_2D_spr.o: param_swipe_2D.c ../models/softplus_reparam.h
	gcc -I../models/ -Wall -c param_swipe_2D.c -o param_swipe_2D_spr.o
	
softplus_reparam.o: ../models/softplus_reparam.h ../models/softplus_reparam.c
	gcc -Wall -c ../models/softplus_reparam.c


spd: param_swipe_2D_spd.o softplus_decay.o
	gcc -Wall param_swipe_2D_spdecay.o softplus_decay.o -lgsl -lgslcblas -lm -o swipe2D_spd

param_swipe_2D_spd.o: param_swipe_2D.c ../models/softplus_decay.h
	gcc -I../models/ -Wall -c param_swipe_2D.c -o param_swipe_2D_spdecay.o

softplus_decay.o: ../models/softplus_decay.h ../models/softplus_decay.c
	gcc -Wall -c ../models/softplus_decay.c

log: param_swipe_2D_log.o logistic.o
	gcc -Wall param_swipe_2D_log.o logistic.o -lgsl -lgslcblas -lm -o swipe2D_log

param_swipe_2D_log.o: param_swipe_2D.c ../models/logistic.h
	gcc -I../models/ -Wall -c param_swipe_2D.c -o param_swipe_2D_log.o

logistic.o: ../models/logistic.h ../models/logistic.c
	gcc -Wall -c ../models/logistic.c


clean:
	@rm -rf *.o
