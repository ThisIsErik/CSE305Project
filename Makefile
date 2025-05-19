CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -Iexternal/thread-pool/include

SRC := $(wildcard src/*.cpp src/parallel/*.cpp src/sequential/*.cpp src/concurrent/*.cpp src/utils/*.cpp)

OUT = run

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(SRC)

clean:
	rm -f $(OUT)