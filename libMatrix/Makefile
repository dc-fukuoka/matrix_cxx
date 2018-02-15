CXX       = $(shell printenv CXX)
CXXFLAGS  = -std=c++11 -g -O3 -mavx -fpic -MMD -MP -Wall
OPENMP    = -fopenmp
ifeq ($(CXX), icpc)
CPPFLAGS  = -I.
LIBS      = -mkl
LDFLAGS   = -shared
else ifeq ($(CXX), g++)
CPPFLAGS  = -I. -I$(MKLROOT)/include
LIBS      = -lmkl_gnu_thread -lmkl_rt -lmkl_core -lmkl_intel_lp64
LDFLAGS   = -shared -L$(MKLROOT)/lib/intel64
endif
SRCS      = Matrix.cxx
OBJS      = $(SRCS:%.cxx=%.o)
DEPS      = $(SRCS:%.cxx=%.d)
BIN       = libMatrix.so

.SUFFIXES: .cxx .o
.PHONY: clean

.cxx.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP) -c $<

$(BIN): $(OBJS)
	$(CXX) $(OPENMP) $(LDFLAGS) $(LIBS) $^ -o $@

ALL: $(BIN)

-include $(DEPS)

clean:
	rm -f *.o *.d *~ $(BIN)
