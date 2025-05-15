CXX = g++
NVCC = nvcc

CXXFLAGS = -std=c++20 -Wall -Wextra -Iexternal/thread-pool/include
NVCCFLAGS = -std=c++20 -O2 -Iexternal/thread-pool/include

SRC := $(wildcard src/*.cpp src/parallel/*.cpp src/sequential/*.cpp src/concurrent/*.cpp src/utils/*.cpp)
CU_SRC := $(wildcard src/gpu/*.cu)

OBJ_CPP := $(SRC:.cpp=.o)
OBJ_CU := $(CU_SRC:.cu=.o)
OBJ_ALL := $(OBJ_CPP) $(OBJ_CU)

OUT = run

all: $(OUT)

$(OUT): $(OBJ_ALL)
	$(NVCC) $(OBJ_ALL) -o $(OUT)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(OUT) $(OBJ_ALL)
