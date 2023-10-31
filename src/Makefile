SAN=
#-fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment

CPPFLAGS=-std=c++17 -flto -O3 

ctyper: main.cpp CramReader.o  FastaReader.o FastqReader.o KmerHash.o KmerCounter.o KmerMatrix.o KmerWindow.o KtableReader.o PriorData.o Regression.o TreeRound.o Processor.o 
	g++ $(SAN) $(CPPFLAGS) -o $@ $^ -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include/ -L/usr/local/lib -L$(CONDA_PREFIX)/lib/ -lhts -lbz2 -lpthread -lz -lm -lstdc++fs  


%.o: %.cpp %.hpp FileReader.hpp
	g++ $(SAN) $(CPPFLAGS) -c $^  -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include/ -L/usr/local/lib -L$(CONDA_PREFIX)/lib/ -lhts -lbz2 -lpthread -lz -lm -lstdc++fs 

clean:
	rm *.o *.gch  ctyper






