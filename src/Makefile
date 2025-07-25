CXX := g++
SAN := # Default empty. Uncomment and set as needed.
# SAN := -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment

BASE_CXXFLAGS := -std=c++17 -O3 -march=native -Wall -Wextra
PGO_GENERATE_FLAGS := -fprofile-generate
PGO_USE_FLAGS := -fprofile-use
LTO_FLAG := -flto # Separated LTO

INCLUDES := -I$(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include
LIBS := -L/usr/local/lib -L$(CONDA_PREFIX)/lib -lhts -lbz2 -lpthread -lz -lm -lstdc++fs

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o) # ALL .o files are in the current directory

TARGET := ctyper
PGO_GEN_EXE := $(TARGET)_pgo_generate # The executable for generating profile data

.PHONY: all clean ctyper_pgo_generate ctyper_pgo_use normal clean_objs

# Default target: builds normal executable
all: ctyper_pgo_use

# --- Object File Compilation Rule ---
# This is the ONE pattern rule for how to build any *.o file.
# The `$(COMPILER_FLAGS)` variable will be set via target-specific variables below.
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(LTO_FLAG) $(INCLUDES) $(COMPILER_FLAGS) -c $< -o $@

# --- Normal Build ---
normal: clean_objs $(TARGET)
	@echo "--------------------------------------------------------"
	@echo "Normal Build Complete."
	@echo "Executable: ./$(TARGET)"
	@echo "--------------------------------------------------------"

# Linking rule for the normal 'ctyper' executable
# Set COMPILER_FLAGS to empty for prerequisites of this target (i.e., normal compilation)
$(TARGET): COMPILER_FLAGS :=
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) (normal build)..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(LTO_FLAG) $(INCLUDES) $^ -o $@ $(LIBS)


# --- PGO Generation Phase ---
ctyper_pgo_generate: clean_objs $(PGO_GEN_EXE)
	@echo "--------------------------------------------------------"
	@echo "PGO Generation Build Complete."
	@echo "Now, run the instrumented executable to generate profile data:"
	@echo "  ./$(PGO_GEN_EXE) <your_typical_workload_arguments>"
	@echo "This will create *.gcda files."
	@echo "--------------------------------------------------------"

# Apply -fprofile-generate to object files when building PGO_GEN_EXE.
# This variable affects the %.o: %.cpp rule when its prerequisites are built.
$(PGO_GEN_EXE): COMPILER_FLAGS := $(PGO_GENERATE_FLAGS)

# Rule to link the PGO generation executable
$(PGO_GEN_EXE): $(OBJS)
	@echo "Linking $(PGO_GEN_EXE) (PGO generate)..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(PGO_GENERATE_FLAGS) $(LTO_FLAG) $(INCLUDES) $^ -o $@ $(LIBS)


# --- PGO Use Phase ---
ctyper_pgo_use: clean_objs $(TARGET)
	@echo "--------------------------------------------------------"
	@echo "PGO Optimization Build Complete."
	@echo "Executable: ./$(TARGET)"
	@echo "--------------------------------------------------------"

# Apply -fprofile-use to object files when building the PGO-optimized TARGET.
# This variable affects the %.o: %.cpp rule when its prerequisites are built.
$(TARGET): COMPILER_FLAGS := $(PGO_USE_FLAGS)

# Rule to link the PGO-optimized 'ctyper' executable
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) (PGO use)..."
	$(CXX) $(SAN) $(BASE_CXXFLAGS) $(PGO_USE_FLAGS) $(LTO_FLAG) $(INCLUDES) $^ -o $@ $(LIBS)


# --- Cleanup Rules ---
clean:
	@echo "Cleaning all..."
	rm -f $(OBJS) *.gch $(TARGET) $(PGO_GEN_EXE) *.profraw *.gcno
	@echo "Removed generated files (preserving *.gcda profile data)."

clean_objs:
	@echo "Cleaning all object files for recompilation..."
	rm -f $(OBJS) *.gch
