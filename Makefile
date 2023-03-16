Ctype: main.cpp CramReader.o FastaReader.o FastqReader.o KmerCounter.o KmerMatrix.o KtableReader.o PriorData.o Regression.o TreeRound.o
	gcc -o $@ $^ -I$(CONDA_PREFIX)/lib/python3.11/site-packages/numpy/core/include -I$(CONDA_PREFIX)/include/python3.11/ -L$(CONDA_PREFIX) -l python3


%.o: %.cpp %.hpp
	gcc -c $^ -I$(CONDA_PREFIX)/include/python3.11/ -I$(CONDA_PREFIX)/lib/python3.11/site-packages/numpy/core/include






