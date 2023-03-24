SAN=
#-fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment


ctyper: main.cpp CramReader.o FastaReader.o FastqReader.o KmerCounter.o KmerMatrix.o KtableReader.o PriorData.o Regression.o TreeRound.o Processor.o
	g++ $(SAN) -std=c++11 -g -o $@ $^ -I $(EIGEN_ROOT)/include/eigen3 -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib -lpthread -lz -L $(CONDA_PREFIX)/lib -lhts 


%.o: %.cpp %.hpp
	g++ $(SAN) -std=c++11 -g -c $^ -I$(CONDA_PREFIX)/include/python3.11/  -I$(CONDA_PREFIX)/include -I $(EIGEN_ROOT)/include/eigen3

clean:
	rm *.o ctyper






