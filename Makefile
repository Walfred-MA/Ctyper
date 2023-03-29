SAN=
#-fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment


ctyper: main.cpp CramReader.o  FastaReader.o FastqReader.o KmerCounter.o KmerMatrix.o KtableReader.o PriorData.o Regression.o TreeRound.o Processor.o
	g++ $(SAN) -std=c++17 -g -O2 -o $@ $^ -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include/ -L$(CONDA_PREFIX)/lib/ -lhts -ldeflate -lbz2 -lpthread -lz -lm 


%.o: %.cpp %.hpp
	g++ $(SAN) -std=c++17 -g -O2 -c $^  -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include/ -L$(CONDA_PREFIX)/lib/ -lhts -ldeflate -lbz2 -lpthread -lz -lm 

clean:
	rm *.o *.gch  ctyper






