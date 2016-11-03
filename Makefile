# msa_astar and msa_pastar Makefile

# Choose -std=c++11 or -std=c++0x
CXXVERSION = $(shell $(CXX) -dumpversion | cut -b 1-3)
#ONLY FOR OPENMPI
#MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
#MPI_LINK_FLAGS = $(shell mpic++ --showme:link)
#CXX = $(shell mpic++ --showme:command)

#FOR MPICH
MPI_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich 
MPI_LINK_FLAGS = -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich

#CXX = clang
ifneq "$(filter g++, $(CXX))" ""
ifeq "$(CXXVERSION)" "4.6"
CPPSTD = -std=c++0x
endif
ifeq "$(VERSION)" "4.4"
$(error Bad $(CXX) version $(CXXVERSION). Atomic operations are required)
endif
endif

ifeq "$(CPPSTD)" ""
CPPSTD = -std=c++11
endif

BIN_DIR     = ./bin
PASTAR_BIN  = $(BIN_DIR)/pastar

TARGET      = $(PASTAR_BIN)

SRC_DIR     = ./pastar
INC_DIR     = ./pastar/include
OBJ_DIR     = ./obj
CPPFLAGS   += -W $(MPI_COMPILE_FLAGS) $(CPPSTD)
LDFLAGS    += -pthread -lstdc++ -lm -lboost_program_options -lboost_system -lboost_filesystem \-lboost_serialization -llz4 $(MPI_LINK_FLAGS) $(CPPSTD)


ifdef THREADS
    CPPFLAGS += -Wall -DTHREADS_NUM=$(THREADS)
endif

ifdef HASH_SHIFT
    CPPFLAGS += -DHASH_SHIFT=$(HASH_SHIFT)
endif

ifdef NO_LIB_BOOST
    CPPFLAGS += -DNO_LIB_BOOST
endif

#ifndef DEBUG
   CPPFLAGS += -O3 -msse3
   LDFLAGS += -s
#else
    #CPPFLAGS += -g
    #LDFLAGS += -g
#endif

ifdef PROFILE_GENERATE
    CPPFLAGS += -fprofile-generate
    LDFLAGS  += -fprofile-generate
endif

ifdef PROFILE_USE
    CPPFLAGS += -fprofile-use
    LDFLAGS += -fprofile-use
endif

ifdef PROFILE_INFORMATION
    CPPFLAGS += -pg
    LDFLAGS += -pg
endif

COMMON_CPP_SRCS += \
    $(SRC_DIR)/backtrace.cpp \
    $(SRC_DIR)/Coord.cpp \
    $(SRC_DIR)/CoordHash.cpp \
    $(SRC_DIR)/Cost.cpp \
    $(SRC_DIR)/HeuristicHPair.cpp \
    $(SRC_DIR)/msa_options.cpp \
    $(SRC_DIR)/Node.cpp \
    $(SRC_DIR)/read_fasta.cpp \
    $(SRC_DIR)/PairAlign.cpp \
    $(SRC_DIR)/Sequences.cpp \
    $(SRC_DIR)/TimeCounter.cpp \
    $(SRC_DIR)/lz4sup.cpp \

PASTAR_SRCS = \
    $(SRC_DIR)/msa_pastar_main.cpp \
    $(SRC_DIR)/pastar.cpp \

INC_PATH += \
    -I$(INC_DIR) \
    -I/usr/include \

CPPFLAGS += \
    $(INC_PATH) \

COMMON_OBJS = $(COMMON_CPP_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
PASTAR_OBJS = $(PASTAR_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

all:	$(TARGET)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<
$(COMMON_OBJS):	| $(OBJ_DIR)
$(PASTAR_OBJS):	| $(OBJ_DIR)

$(PASTAR_BIN):	$(COMMON_OBJS) $(PASTAR_OBJS) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET) $(COMMON_OBJS) $(PASTAR_OBJS) 
