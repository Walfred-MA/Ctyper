SAN=
#-fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment

CPPFLAGS=-std=c++17 -flto -O3 -march=native 
PGO_GENERATE_FLAGS=-fprofile-generate
PGO_USE_FLAGS=-fprofile-use


ctyper: main.cpp CramReader.o  FastaReader.o FastqReader.o KmerHash.o KmerCounter.o KmerMatrix.o KmerWindow.o KtableReader.o PriorData.o Regression.o TreeRound.o Processor.o 
	$(CXX) $(SAN) $(CPPFLAGS) -o $@ $^ -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include/ -L/usr/local/lib -L$(CONDA_PREFIX)/lib/ -lhts -lbz2 -lpthread -lz -lm -lstdc++fs  


ctyper_pgo_generate: CPPFLAGS += $(PGO_GENERATE_FLAGS)
ctyper_pgo_generate: ctyper

ctyper_pgo_use: CPPFLAGS += $(PGO_USE_FLAGS)
ctyper_pgo_use: ctyper

%.o: %.cpp %.hpp FileReader.hpp
	$(CXX) $(SAN) $(CPPFLAGS) -c $^  -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include/ -L/usr/local/lib -L$(CONDA_PREFIX)/lib/ -lhts -lbz2 -lpthread -lz -lm -lstdc++fs 

clean:
	rm *.o *.gch  ctyper






