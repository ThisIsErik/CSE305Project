CXX = g++
NVCC = nvcc

# Flags
CXXFLAGS = -std=c++20 -Wall -Wextra -Iexternal/thread-pool/include $(DEFINES)
NVCCFLAGS = -std=c++20 -O2 -Iexternal/thread-pool/include

# Sources
SRC := $(wildcard src/*.cpp src/parallel/*.cpp src/sequential/*.cpp src/concurrent/*.cpp src/utils/*.cpp)
CU_SRC := $(wildcard src/gpu/*.cu)

# Objects
OBJ_CPP := $(SRC:.cpp=.o)
OBJ_CU := $(CU_SRC:.cu=.cu.o)
OBJ_ALL := $(OBJ_CPP) $(OBJ_CU)

# Output executables
OUT = run
CPU_OUT = run_cpu

# ====================
# Targets
# ====================

# Default GPU + CPU build
all: DEFINES = -DUSE_CUDA
all: $(OUT)

# CPU-only build (no CUDA)
cpu: $(OBJ_CPP)
	$(CXX) $(CXXFLAGS) $(OBJ_CPP) -o $(CPU_OUT)

# Link full build (with GPU)
$(OUT): $(OBJ_ALL)
	$(CXX) $(CXXFLAGS) $(OBJ_ALL) -o $(OUT)

# Compile C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA files
%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Clean everything
clean:
	rm -f $(OUT) $(CPU_OUT) $(OBJ_ALL)
