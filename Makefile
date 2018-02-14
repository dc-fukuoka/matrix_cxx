all:
	icpc -std=c++11 -Wall -g -O3 -mavx -fopenmp test.cxx -mkl
#	g++ -I$(MKLROOT)/include -std=c++11 -Wall -g -O3 -mavx -fopenmp test.cxx -L$(MKLROOT)/lib/intel64 -lmkl_gnu_thread -lmkl_rt -lmkl_core -lmkl_intel_lp64

clean:
	rm -f a.out *~
