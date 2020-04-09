##########################################################
##########################################################
 PROGRAM = 2ddg
##########################################################
# Suffix Rule for f90
##########################################################
.SUFFIXES : .o .f90
.f90.o:
	gfortran -O -c  -fcheck=all -fbacktrace $< 
##########################################################
SDIR = .
OBCTS =	$(SDIR)/dgflo2d.o\
        $(SDIR)/Postprocess.o\
        $(SDIR)/Preprocess.o\
        $(SDIR)/Getrhs.o\
        $(SDIR)/Getrhslag.o\
        $(SDIR)/Getrhslagmc.o\
        $(SDIR)/Reconstruction.o\
        $(SDIR)/Limiterlag.o\
        $(SDIR)/Getrhslagmctria.o\
        $(SDIR)/Getrhslaggauss.o\
        $(SDIR)/Lagfvm.o\
        $(SDIR)/Getrhslagrz.o\
        $(SDIR)/GetBClag.o\
        $(SDIR)/Getrhslaggd.o\
        $(SDIR)/LinearAlgebra.o\
        $(SDIR)/Getsolid.o\
        $(SDIR)/GetInitial.o\
        $(SDIR)/GetFE.o\
        $(SDIR)/Getrhslag_curv.o\
        $(SDIR)/Getrhslag_linear.o\
        $(SDIR)/WENO.o\
##########################################################
# Make executable "Hyperdg.o"
##########################################################
$(PROGRAM): $(OBCTS)
	gfortran -O -o $@ $(OBCTS) 
##########################################################
# Clean up: Type 'make clean'
##########################################################
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.M
##########################################################
# Clean up: Type 'make clean~'
##########################################################
clean~:
	rm -f *~
