# Helmholtz and Laplace particle FMMs in R^2

PROJECT_LIB=hfmm2d




OBJSUF=o
MODSUF=mod
FC=h5fc -c 
FFLAGS=-O3 -fopenmp -fallow-argument-mismatch
FLINK=h5fc -o -fopenmp $(PROJECT)




.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c

# SOURCE FILE LIST
#
vpath %.f .

FMODS = 

FSRCS = d2tstrcr_omp.f d2mtreeplot.f  \
     cdjseval2d.f  

FSRCS += prini.f 

FSRCS += cfmm2dpart.f rfmm2dpart.f lfmm2drouts.f \
     l2dterms.f laprouts2d.f


ifeq ($(WITH_SECOND),1) 
FSRCS += second-r8.f
endif

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF)) 
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF)) 
#

lib$(PROJECT_LIB).a:   $(MODS)   $(OBJS)
	ar rvs lib$(PROJECT_LIB).a $(OBJS)

#
clean: 
	rm -f $(OBJS)
# 
distclean: 
	rm -f $(OBJS) int2 *.a fort.*
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
