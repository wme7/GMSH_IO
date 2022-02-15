# Detect Linux of MacOS
OS := $(shell uname)
PC := $(shell hostname)

# Make Directory command
MKDIR_P := mkdir -p

# Program subfolders
APP_DIR = src
EXE_DIR = build
INC_DIR = include
dep = $(obj:.o=.d)

# Edit the following variables as needed
EXTLIB  = -L$(HDF_INSTALL)/lib 
CC      = g++
CFLAGS  = -O3 -Wall -Wextra  -std=c++17 
LDFLAGS =
LIB     = -lz -lm -ldl 
INCLUDE = -I$(INC_DIR) 

# APPs to be compiled
EXE_MSH1=$(EXE_DIR)/TestParsers.run

# Compilation process for:
all: tests
tests: $(EXE_DIR) $(EXE_MSH1)

$(EXE_DIR):
	${MKDIR_P} $(EXE_DIR)

$(EXE_MSH1): $(APP_DIR)/TestParsers.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(INCLUDE) $(LIBSHDF) $(LIB)

$(APP_DIR)/%.o: $(APP_DIR)/%.cpp
	$(CC) -o $@ -c -MMD $< $(CFLAGS) $(INCLUDE) 

-include $(dep) 

%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

clean:
	rm -rf $(EXE_DIR)/*.run $(APP_DIR)/*.o $(APP_DIR)/*.d $(EXE_DIR)/*.h5 $(EXE_DIR)/*.xmf

.SUFFIXES: .o .cpp
