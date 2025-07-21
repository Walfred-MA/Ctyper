# === Compiler and base flags ===
CXX := g++
CXXFLAGS := -std=c++17 -O3 -march=native -flto -Wall -Wextra
SAN :=
# SAN := -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all

# === Include and library paths ===
INCLUDES := -I$(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include
LIBS := -L/usr/local/lib -L$(CONDA_PREFIX)/lib -lhts -lbz2 -lpthread -lz -lm -lstdc++fs

# === Auto-detect all .cpp files ===
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

# === Output target ===
TARGET := ctyper

# === Phony targets ===
.PHONY: all clean ctyper_pgo_generate ctyper_pgo_use

# === Default build ===
all: $(TARGET)

# === Linking step ===
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SAN) -o $@ $^ $(INCLUDES) $(LIBS)

# === Compile each .cpp to .o ===
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(SAN) -c $< -o $@ $(INCLUDES)

# === PGO generation ===
ctyper_pgo_generate:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fprofile-generate" $(TARGET)

# === PGO usage ===
ctyper_pgo_use:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fprofile-use" $(TARGET)

# === Clean all intermediate files ===
clean:
	@echo "Cleaning..."
	rm -f *.o *.gch $(TARGET)

