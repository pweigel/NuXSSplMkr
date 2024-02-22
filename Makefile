#PATHS
ifeq ($(PREFIX),)
PREFIX=/usr/local/
endif

CURRENT_DIR 	= $(shell pwd)
# LHAPDF      	= $(PREFIX)
# BOOST       	= $(PREFIX)
# PHOTOSPLINE 	= $(PREFIX)

SOURCES 	= $(wildcard src/*.cpp)

OBJECTS 	= $(SOURCES:.cpp=.o)

INCLUDE_PATH 	= -I/usr/local/include -I./inc
INCLUDE_PATH  	+= -I$(CURRENT_DIR)/src
INCLUDE_PATH  	+= -I$(PREFIX)/include
INCLUDE_PATH    += -I$(SROOT)/include
INCLUDE_PATH    += -I$(PREFIX)/include/suitesparse  # TODO: Try to figure this out, needed for -

#Compiler
# CC 		= clang
# CXX 		= clang++
CC 		= gcc
CXX 		= g++

#Dynamic Library

#Flags
CXX_FLAGS       =  $(INCLUDE_PATH) -I. -O3 -fPIC -std=c++17 -g -Wall

# LD 		= clang++
LD 		= g++
LD_FLAGS 	= -L/usr/local/lib/ -L/usr/lib  -L/usr/local/lib64 -L/usr/lib64 -L$(PREFIX)/lib -L$(PREFIX)/lib64
LD_FLAGS  += -L$(SROOT)/lib
LD_FLAGS  += -L$(SROOT)/lib64

LD_FLAGS 	+= -lLHAPDF
LD_FLAGS 	+= -lgsl -lgslcblas
LD_FLAGS	+= -lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex 
# New features:
LD_FLAGS  += -lAPFEL -lapfelxx -lboost_program_options -lspglam -lphotospline


# Here are some things that need to be figured out...
# I have to do this to get photospline::ndsparse
INCLUDE_PATH += -I$(PREFIX)/include/suitesparse
LD_FLAGS += -lcfitsio -lspglam -lcholmod
CXX_FLAGS += -DPHOTOSPLINE_INCLUDES_SPGLAM

.PHONY: all clean

# CT_OBJ = $(CT)src/CT12Pdf.o $(CURRENT_DIR)src/ct10_xs.o

# all: bin/make_sf_splines bin/make_all_sf_splines bin/calculate_xs bin/calculate_all_xs bin/calculate_dsdy
all: bin/make_sf_splines bin/make_all_sf_splines bin/make_all_sf_splines_replicas bin/calculate_xs bin/calculate_all_xs bin/calculate_dsdy bin/construct_fonll bin/calculate_dsdxdy bin/calculate_dsdxdQ2 bin/calculate_all_dsdy bin/calculate_all_dsdy_replicas
# test: bin/test_CKMT bin/test_TMC
test: bin/test_grid bin/APFEL_check bin/test_dsdxdy #bin/apfelxx_test
replicas: bin/make_all_sf_splines_replicas bin/calculate_all_dsdy_replicas bin/calculate_all_xs_replicas

# bin/calculate_LO_xs: src/configuration.o src/structure_function.o src/physconst.o mains/calculate_LO_xs.o
# 	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
bin/calculate_dsdxdQ2: src/tools.o src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_dsdxdQ2.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_dsdxdy: src/tools.o src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_dsdxdy.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_dsdy: src/tools.o src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_dsdy.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_all_dsdy_replicas: src/tools.o src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_all_dsdy_replicas.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_all_dsdy: src/tools.o src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_all_dsdy.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_all_xs_replicas: src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_all_xs_replicas.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_xs: src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_xs.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/calculate_all_xs: src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o mains/calculate_all_xs.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

bin/make_all_sf_splines_replicas: src/configuration.o src/structure_function.o src/physconst.o mains/make_all_sf_splines_replicas.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/make_all_sf_splines: src/configuration.o src/structure_function.o src/physconst.o mains/make_all_sf_splines.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/make_sf_splines: src/configuration.o src/structure_function.o src/physconst.o mains/make_sf_splines.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/construct_fonll: src/configuration.o src/structure_function.o src/physconst.o mains/construct_fonll.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/test_CKMT: src/configuration.o src/structure_function.o src/physconst.o tests/test_CKMT.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
bin/test_TMC: src/configuration.o src/structure_function.o src/physconst.o tests/test_TMC.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
bin/test_grid: src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o src/tools.o tests/test_grid.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
bin/APFEL_check: tests/APFEL_check.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
bin/apfelxx_test: tests/apfelxx_test.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
bin/test_dsdxdy: src/configuration.o src/structure_function.o src/cross_section.o src/physconst.o tests/test_dsdxdy.o
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@
%.o:%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

clean:
	rm src/*.o bin/*.exe mains/*.o
