# HDF_INSTALL = /usr/local/Cellar/hdf5/1.10.4#/anaconda3#/usr/local/Cellar/hdf5/1.10.4#/Users/marchiano/src/hdf5#/usr/local/Cellar/hdf5/1.10.4#/usr/lib/x86_64-linux-gnu/hdf5/serial#
# EIGEN_INSTALL = /Users/marchiano/src/eigen

# Detect Linux of MacOS
OS := $(shell uname)
PC := $(shell hostname)
ifeq ($(OS),Darwin)
	ifeq ($(PC), Manuels-MacBook.local)
		HDF_INSTALL = /usr/local/Cellar/hdf5/1.12.0_3
		EIGEN_INSTALL = /usr/local/Cellar/eigen/3.4/include/eigen3
	endif
else
	ifeq ($(PC), B17-postdoc2020)
		HDF_INSTALL = /mnt/c/Users/mdiazesc/phdf5
		EIGEN_INSTALL = /mnt/c/Users/mdiazesc/eigen3/include/eigen3
	endif
endif

# Make Directory command
MKDIR_P := mkdir -p

# Program subfolders
APP_DIR = src
EXE_DIR = build
INC_DIR = include
#OBJ=$(SRC:.cpp=.o)
dep = $(obj:.o=.d)

# Edit the following variables as needed
EXTLIB  = -L$(HDF_INSTALL)/lib 
# CC      = g++
CC      = mpicxx
CFLAGS  = -O3 -Wall -Wextra  -std=c++11 
LDFLAGS =
LIB     = -lz -lm -ldl 
INCLUDE = -I$(INC_DIR) -I$(HDF_INSTALL)/include -I$(EIGEN_INSTALL)
#LIBSHDF = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a
#LIBSFFTW= $(EXTLIB) $(FFTW_INSTALL)/lib/libfftw3.a $(FFTW_INSTALL)/lib/libfftw3-mpi.a
LIBSHDF = $(EXTLIB) -lhdf5

# APPs to be compiled
EXE_MSH1=$(EXE_DIR)/test_GMSHparserV2_Global.run
EXE_MSH2=$(EXE_DIR)/test_GMSHparserV4_Global.run

# Compilation process for:
all: tests
tests: $(EXE_DIR) $(EXE_MSH1) $(EXE_MSH2) 

$(EXE_DIR):
	${MKDIR_P} $(EXE_DIR)

$(EXE_MSH1): $(APP_DIR)/test_GMSHparserV2_Global.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(INCLUDE) $(LIBSHDF) $(LIB)

$(EXE_MSH2): $(APP_DIR)/test_GMSHparserV4_Global.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(INCLUDE) $(LIBSHDF) $(LIB)

$(APP_DIR)/%.o: $(APP_DIR)/%.cpp
	$(CC) -o $@ -c -MMD $< $(CFLAGS) $(INCLUDE) 

-include $(dep) 

%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

clean:
	rm -rf $(EXE_DIR)/*.run $(APP_DIR)/*.o $(APP_DIR)/*.d $(EXE_DIR)/*.h5 $(EXE_DIR)/*.xmf

.SUFFIXES: .o .cpp

# Dependencies
