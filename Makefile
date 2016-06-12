# msa_astar and msa_pastar Makefile

# Choose -std=c++11 or -std=c++0x
#CXXVERSION = $(shell $(CXX) -dumpversion | cut -b 1-3)
BOOST_COMPILE_FLAGS = -I"C:\tools\boost_1_61_0\"
BOOST_LINK_FLAGS = /LIBPATH:"C:\tools\boost_1_61_0\lib64_msvc_14.0" -boost_program_options -lboost_system -lboost_filesystem \-lboost_mpi -lboost_serialization
CXX = /c/Program\ Files\ \(x86\)/IntelSWTools/mpi/5.1.3.207/intel64/bin/mpicxx.bat
#ifneq "$(filter g++, $(CXX))" ""
#ifeq "$(CXXVERSION)" "4.6"
#CPPSTD = -std=c++0x
#endif
#ifeq "$(VERSION)" "4.4"
#$(error Bad $(CXX) version $(CXXVERSION). Atomic operations are required)
#endif
#endif

#ifeq "$(CPPSTD)" ""
#CPPSTD = -std=c++11
#endif

BIN_DIR     = ./bin
PASTAR_BIN  = $(BIN_DIR)/msa_pastar

TARGET      = $(PASTAR_BIN)

SRC_DIR     = ./pastar
INC_DIR     = ./pastar
OBJ_DIR     = ./obj
CPPFLAGS   += -W -g $(BOOST_COMPILE_FLAGS) -pthread $(CPPSTD)
LDFLAGS    += -g -pthread -lstdc++ -lm  $(BOOST_LINK_FLAGS) $(CPPSTD)


ifdef THREADS
    CPPFLAGS += -Wall -DTHREADS_NUM=$(THREADS)
endif

ifdef HASH_SHIFT
    CPPFLAGS += -DHASH_SHIFT=$(HASH_SHIFT)
endif

ifdef NO_LIB_BOOST
    CPPFLAGS += -DNO_LIB_BOOST
endif

ifndef DEBUG
    OPTIMIZE = yes
    LDFLAGS += -s
else
    CPPFLAGS += -g
endif

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

PASTAR_SRCS = \
    $(SRC_DIR)/msa_pastar_main.cpp \
    $(SRC_DIR)/PAStar.cpp \

INC_PATH += \
    -I$(INC_DIR) \

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
