CXX := g++
SAN := # Default empty. Uncomment and set as needed.
# SAN := -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment

BASE_CXXFLAGS := -std=c++17 -O3 -march=haswell -mtune=generic -Wall -Wextra
PGO_GENERATE_FLAGS := -fprofile-generate
PGO_USE_FLAGS := -fprofile-use
LTO_FLAG := -flto # Separated LTO

INCLUDES := -I$(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include
LIBS := -L/usr/local/lib -L$(CONDA_PREFIX)/lib -lhts -lbz2 -lpthread -lz -lm -lstdc++fs

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

TARGET := ctyper
PGO_GEN_EXE := $(TARGET)-pgo-gen

.PHONY: all clean ctyper_pgo_generate ctyper_pgo_use normal clean_objs

all: ctyper_pgo_use

%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(LTO_FLAG) $(INCLUDES) $(COMPILER_FLAGS) -c $< -o $@

normal: clean_objs $(TARGET)
	@echo "--------------------------------------------------------"
	@echo "Normal Build Complete."
	@echo "Executable: ./$(TARGET)"
	@echo "--------------------------------------------------------"

$(TARGET): COMPILER_FLAGS :=
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) (normal build)..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(LTO_FLAG) $(INCLUDES) $^ -o $@ $(LIBS)

# --- PGO Generation Phase ---
ctyper_pgo_generate: clean_objs $(PGO_GEN_EXE)
	@echo "--------------------------------------------------------"
	@echo "PGO Generation Build Complete."
	@echo "Now run the instrumented executable to generate profile data:"
	@echo "  GCOV_PREFIX='$(CURDIR)' GCOV_PREFIX_STRIP=100 ./$(PGO_GEN_EXE) <args>"
	@echo "--------------------------------------------------------"

$(PGO_GEN_EXE): COMPILER_FLAGS := $(PGO_GENERATE_FLAGS)
$(PGO_GEN_EXE): $(OBJS)
	@echo "Linking $(PGO_GEN_EXE) (PGO generate)..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(PGO_GENERATE_FLAGS) $(LTO_FLAG) $(INCLUDES) $^ -o $@ $(LIBS)

# --- PGO Use Phase ---
ctyper_pgo_use: clean_objs $(TARGET)
	@echo "--------------------------------------------------------"
	@echo "PGO Optimization Build Complete."
	@echo "Executable: ./$(TARGET)"
	@echo "--------------------------------------------------------"

$(TARGET): COMPILER_FLAGS := $(PGO_USE_FLAGS)
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) (PGO use)..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(PGO_USE_FLAGS) $(LTO_FLAG) $(INCLUDES) $^ -o $@ $(LIBS)

# --- Cleanup Rules ---
clean:
	@echo "Cleaning all..."
	rm -f $(OBJS) *.gch $(TARGET) $(PGO_GEN_EXE) *.profraw *.gcno
	@echo "Removed generated files."

clean_objs:
	@echo "Cleaning all object files for recompilation..."
	rm -f $(OBJS) *.gch

clean_all:
	@echo "Cleaning everything including PGO data..."
	rm -f $(OBJS) *.gch $(TARGET) $(PGO_GEN_EXE) *.profraw *.gcno *.gcda
