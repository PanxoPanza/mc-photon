mathlib_lib = ./mathlib
mathlib_dep = $(patsubst %,$(mathlib_lib)/%,localmathlib.h) $(patsubst %,$(mathlib_lib)/%,uni_localmathlib.h)

# 2.optical properties
optprop_lib = ./optical_properties
optprop_dep = $(patsubst %,$(optprop_lib)/%,materials.h)

DEPS = $(optprop_dep) $(mathlib_dep) MCRT_library.h

# -----------------------------------------------------------------
# objects directory
ODIR = obj

# compiling instructions
CC = g++ -std=c++11 -fopenmp -g
CFLAGS = -I$(mathlib_lib) -I$(optprop_lib) -I.

# object files
_OBJ = main.o ESource.o mcrtMiscellaneous.o Mie_Inclusion.o \
	ModelBuild.o Panel.o Photon.o Point3D.o preder_output.o \
	ReadGMSHFile.o Region.o Surface.o Vector3D.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# -----------------------------------------------------------------
# main target
mc_photon: $(OBJ) $(ODIR)/localmathlib.o $(ODIR)/uni_localmathlib.o $(ODIR)/materials.o
	$(CC) -o $@ $^ $(CFLAGS)

# object targets
# 1.compile dependencies

$(ODIR)/localmathlib.o: $(mathlib_lib)/localmathlib.cpp $(mathlib_dep)
	$(CC) -c -o $@ $< -I$(mathlib_lib)

$(ODIR)/uni_localmathlib.o: $(mathlib_lib)/uni_localmathlib.cpp $(mathlib_dep)
	$(CC) -c -o $@ $< -I$(mathlib_lib)

$(ODIR)/materials.o: $(optprop_lib)/materials.cpp $(optprop_dep) $(mathlib_dep)
	$(CC) -c -o $@ $< -I$(optprop_lib) -I$(mathlib_lib)

# 2. compile rest of the objects
$(ODIR)/%.o: %.cpp $(DEPS)
	$(shell [ ! -d $(@D) ] && mkdir -p $(@D))
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -rf $(ODIR) mc_photon

all:
	@echo $(DEPS)
	@echo $(OBJ)

MKDIR_P = mkdir -p

.PHONY: directories

all: directories program

directories: ${ODIR}

${ODIR}:
	${MKDIR_P} ${ODIR}
