# Makro-Definitions:

OBJS = hpair.o Cteq61Pdf.o
OBJS0 = hpair0.o Cteq61Pdf.o

#LIBS = /opt/lhapdf/lib/libLHAPDF.so

LIBS = -L/opt/lhapdf/lib/ -lLHAPDF -Wl,-rpath=/opt/lhapdf/lib/

#FFLAGS= -fno-emulate-complex -fno-automatic -ffixed-line-length-none -ffast-math -march=pentiumpro -malign-double -Wall -fno-silent

#FC=f77

#FFLAGS= -fno-fixed-form -fno-automatic -ffixed-line-length-none -ffast-math -march=pentiumpro -malign-double -Wall

#FFLAGS=

FC=gfortran

#FFLAGS= -pc 64 -g77libs

#FC=pgf77

# Commands:
.f.o:
	$(FC) -c $(FFLAGS) $*.f

hpair: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o run

hpair0: $(OBJS0)
	$(FC) $(FFLAGS) $(OBJS0) $(LIBS) -o run0

clean:
	rm -f $(OBJS)
