#FLAGS = -O3 -Wall -I. -fnested-functions
FLAGS = -I. -O3 -Wall -Wno-misleading-indentation

fixtemp:main_fixtemp.c geometry.o energy.o obs.o utils.o misc.o \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o main main_fixtemp.c geometry.o energy.o obs.o \
	utils.o misc.o -lm 

simtemp:main_simtemp.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o main main_simtemp.c geometry.o energy.o obs.o \
	utils.o misc.o -lm 

misc.o: misc.c global.h defs.h sys.h main.h
	gcc $(FLAGS) -c misc.c -lm

geometry.o: geometry.c global.h defs.h sys.h main.h
	gcc $(FLAGS) -c geometry.c

energy.o: energy.c global.h defs.h sys.h main.h
	gcc $(FLAGS) -c energy.c

obs.o: obs.c global.h defs.h sys.h main.h
	gcc $(FLAGS) -c obs.c

utils.o: utils.c global.h defs.h sys.h main.h
	gcc $(FLAGS) -c utils.c

#####################

constants: tools/constants.c
	gcc $(FLAGS) -o constants tools/constants.c

conf2pdb:tools/conf2pdb.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -lm -o conf2pdb tools/conf2pdb.c geometry.o \
	energy.o obs.o utils.o misc.o -lm 

analys_conf_MBAR:tools/analys_conf_MBAR.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o analys_conf_MBAR tools/analys_conf_MBAR.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

analys_conf_MBAR_err:tools/analys_conf_MBAR_err.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o analys_conf_MBAR_err tools/analys_conf_MBAR_err.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

regen_checkpnt:tools/regen_checkpnt.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o regen_checkpnt tools/regen_checkpnt.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

conf2conf:tools/conf2conf.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o conf2conf tools/conf2conf.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

analys_conf:tools/analys_conf.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o analys_conf tools/analys_conf.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

analys_rt:tools/analys_rt.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o analys_rt tools/analys_rt.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

gfix:	tools/gfix.c defs.h
	gcc $(FLAGS) -o gfix tools/gfix.c -lm

gcalc_averages:	tools/gcalc_averages.c defs.h sys.h 
	gcc $(FLAGS) -o gcalc_averages tools/gcalc_averages.c -lm

rt_gcalc: tools/rt_gcalc.c defs.h 
	gcc $(FLAGS) -o rt_gcalc tools/rt_gcalc.c -lm

analys:	tools/analys.c 
	gcc $(FLAGS) -lm -o analys tools/analys.c 

analys_his:tools/analys_his.c 
	gcc $(FLAGS) -lm -o analys_his tools/analys_his.c 

analys_cv:tools/analys_cv.c
	gcc $(FLAGS) -o analys_cv tools/analys_cv.c -lm

analys_multihist:tools/analys_multihist.c 
	gcc $(FLAGS) -o analys_multihist tools/analys_multihist.c  -lm

analys_multihist_2d:tools/analys_multihist_2d.c 
	gcc $(FLAGS) -o analys_multihist_2d tools/analys_multihist_2d.c -lm

analys_his2d:tools/analys_his2d.c 
	gcc $(FLAGS) -lm -o analys_his2d tools/analys_his2d.c 

analys_singlehist:tools/analys_singlehist.c 
	gcc $(FLAGS) -o analys_singlehist tools/analys_singlehist.c -lm

superimpose_list_ca:tools/superimpose_list_ca.c
	gcc $(FLAGS) -o superimpose_list_ca tools/superimpose_list_ca.c -lm

clust_depth_first:tools/clust_depth_first.c geometry.o energy.o obs.o utils.o misc.o  \
	defs.h sys.h global.h main.h
	gcc $(FLAGS) -o clust_depth_first tools/clust_depth_first.c geometry.o \
	energy.o obs.o utils.o misc.o -lm

#####################

clean:
	rm -f *.o
	rm -f main 
	rm -f analys_conf

cleandata:
	rm -f _* tmp* rmin* emin* init_* snapshots/*

cleanall:
	rm -f *.o
	rm -f *~
	rm -f _* tmp* rmin* emin* init_*  
	rm -f main 
	rm -f constants 
