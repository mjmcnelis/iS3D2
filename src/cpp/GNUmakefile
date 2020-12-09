# ===========================================================================
#  Makefile iS3D
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile
##
##  Usage : 	(g)make	[all]		compile the whole project
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE
##			distclean	remove all objectsfiles and binaries
##

#choose a compiler.
#CC = pgc++ # pgi compiler for openacc support
CC = g++
#CC = icpc # icpc slower than g++ for sampler? (10/16/18)

#choose flags corresponding to compiler
#CFLAGS= -std=c++11 -O3 -fopenmp -lgsl -lgslcblas -lm # for g++
CFLAGS= -std=c++11 -O3 -lgsl -lgslcblas -lm # for g++
#CFLAGS= -std=c++11 -qopenmp -O3 -lgsl -lgslcblas # for icpc

RM = rm -f
O = .o
LDFLAGS = $(CFLAGS)
SYSTEMFILES = $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN = iS3D.e
endif

SRC = Main.cpp iS3D.cpp Arsenal.cpp emissionfunction.cpp MomentumSpectra.cpp SpacetimeDistribution.cpp ParticleSampler.cpp Polarization.cpp Table.cpp readindata.cpp ParameterReader.cpp DeltafData.cpp GaussThermal.cpp viscous_correction.cpp

INC = iS3D.h Arsenal.h emissionfunction.h Table.h readindata.h Stopwatch.h ParameterReader.h deltafReader.h GaussThermal.h viscous_correction.h Macros.h Particle.h


# -------------------------------------------------

OBJDIR = obj
SRCFILES = $(SRC) $(INC) GNUmakefile
OBJECTS = $(addprefix $(OBJDIR)/, $(addsuffix $O, \
		  $(basename $(SRC))))
TARGET = $(MAIN)
INSTPATH = ../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:	all mkobjdir clean distclean install

all:	mkobjdir $(TARGET)

help:
	@grep '^##' GNUmakefile

mkobjdir:
	-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

clean:
	-rm $(OBJECTS)

distclean:
	-rm $(TARGET)
	-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
# ./Main.cpp: iS3D.h
# ./iS3D.cpp: iS3D.h Table.h readindata.h emissionfunction.h arsenal.h ParameterReader.h
# ./arsenal.cpp: arsenal.h
# ./emissionfunction.cpp: iS3D.h readindata.h emissionfunction.h Stopwatch.h arsenal.h Table.h gaussThermal.h
# ./emissionfunction_momentum_spectra.cpp: iS3D.h readindata.h emissionfunction.h Stopwatch.h arsenal.h Table.h gaussThermal.h
# ./emissionfunction_spacetime_distribution.cpp: iS3D.h readindata.h emissionfunction.h Stopwatch.h arsenal.h Table.h gaussThermal.h
# ./emissionfunction_polarization.cpp: iS3D.h readindata.h emissionfunction.h Stopwatch.h arsenal.h Table.h gaussThermal.h
# ./emissionfunction_particle_sampler.cpp: iS3D.h readindata.h emissionfunction.h Stopwatch.h arsenal.h Table.h gaussThermal.h viscous_correction.h
# ./Table.cpp: Table.h arsenal.h
# ./readindata.cpp: iS3D.h readindata.h arsenal.h Table.h ParameterReader.h
# ./ParameterReader.cpp: arsenal.h ParameterReader.h
# ./gaussThermal.cpp: gaussThermal.h
##./emissionfunction_resonance_decays.cpp: iS3D.h readindata.h emissionfunction.h Stopwatch.h arsenal.h ParameterReader.h
