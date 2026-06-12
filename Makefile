CXX ?= g++
CPPFLAGS = -I. -I/usr/local/include

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CXX := clang++
	CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pthread
	LDFLAGS  = -g -L/usr/local/lib
	LDLIBS   = -lm -pthread
else
	CXX := $(CXX)
	CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -march=native -mtune=native -DNDEBUG
	LDFLAGS  = -g -L/usr/local/lib
	LDLIBS   = -lm
endif

DEPS = graph.h utility.h algorithms.h
OBJ  = graph.o utility.o algorithms.o main.o \
	kcored/noise.o kcored/lds.o kcored/kcore_ldp.o

.PHONY: all clean distclean

all: counting

%.o: %.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

counting: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f counting *.o kcored/*.o

distclean: clean
	rm -f results/facebook_run.log results/facebook_metrics.csv figures/facebook_metrics.png figures/facebook_metrics.pdf
