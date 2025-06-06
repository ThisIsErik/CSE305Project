CXX = g++
NVCC = /usr/local/cuda/bin/nvcc

# Compiler flags
CXXFLAGS = -std=c++20 -Wall -Wextra -Iexternal/thread-pool/include -Isrc $(DEFINES)
CUDA_HOME ?= /usr/local/cuda
ARCH = -arch=sm_60
NVCCFLAGS = -std=c++20 -O2 -Iexternal/thread-pool/include -I$(CUDA_HOME)/include -Isrc $(ARCH) $(DEFINES)

# Sources
SRC := $(wildcard src/*.cpp src/parallel/*.cpp src/sequential/*.cpp src/concurrent/*.cpp src/utils/*.cpp src/parallel/mm/*.cpp)
CU_SRC := $(wildcard src/gpu/*.cu)

# Object files
OBJ_CPP := $(SRC:.cpp=.o)
OBJ_CU := $(CU_SRC:.cu=.o)
OBJ_ALL := $(OBJ_CPP) $(OBJ_CU)

# Executables
OUT = run
CPU_OUT = run_cpu

# ====================
# Targets
# ====================

# Full GPU + CPU build
all: DEFINES = -DUSE_CUDA
all: $(OUT)

# CPU-only build (safe anywhere)
cpu: $(OBJ_CPP)
	$(CXX) $(CXXFLAGS) $(OBJ_CPP) -o $(CPU_OUT)

# Link everything with nvcc in the full build
$(OUT): $(OBJ_ALL)
	$(NVCC) $(ARCH) $(OBJ_ALL) -o $(OUT)

# Compile C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA files
src/gpu/%.o: src/gpu/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OUT) $(CPU_OUT) $(OBJ_ALL)