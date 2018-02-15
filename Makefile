CXX       = $(shell printenv CXX)
CXXFLAGS  = -std=c++11 -g -O3 -mavx -MMD -MP -Wall
OPENMP    = -fopenmp
CPPFLAGS  = -I./libMatrix
LIBS      = -lMatrix
LDFLAGS   = -L./libMatrix -Wl,-rpath=$(PWD)/libMatrix
SRCS      = main.cxx
OBJS      = $(SRCS:%.cxx=%.o)
DEPS      = $(SRCS:%.cxx=%.d)
BIN       = a.out

.SUFFIXES: .cxx .o
.PHONY: clean

.cxx.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

$(BIN): $(OBJS) $(LIBS)
	$(CXX) $(LDFLAGS) $(LIBS) $^ -o $@

$(LIBS):
	make -C libMatrix

ALL: $(BIN)

-include $(DEPS)

clean:
	rm -f *.o *.d *~ $(BIN)
	make -C libMatrix clean
