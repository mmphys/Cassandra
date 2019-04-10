# CERN root library config utility
APFEL     ?= $(HOME)/bin/apfel_gnu
RC        := root-config
CXX       ?= g++
#CXX       ?= clang++

CXXFLAGS  := $(shell $(RC) --cflags)
CXXFLAGS  := $(filter -I%, $(CXXFLAGS))
CXXFLAGS  := -std=c++17 $(CXXFLAGS) -I/usr/local/include -I$(APFEL)/include

LDFLAGS    := $(shell $(RC) --libs)
LDFLAGS    := $(filter -L%, $(LDFLAGS))
LDFLAGS   := $(LDFLAGS) $(addprefix -L,$(sort $(shell lhapdf-config --libdir) $(APFEL)/lib $(ROOTSYS)))

LDLIBS     := $(filter -l%, $(shell $(RC) --libs))
LDLIBS    := $(addprefix -l,APFEL LHAPDF) $(LDLIBS)

srcfiles  := $(shell find . -type f -iname "*.cpp")
objects   := $(patsubst %.cpp, %.o, $(srcfiles))

all: Cassandra Twiggy

Cassandra: Cassandra.o CassMain.o
	$(CXX) $(CXXFLAGS) -o Cassandra Cassandra.o CassMain.o $(LDFLAGS) $(LDLIBS)

Twiggy: Cassandra.o Twiggy.o
	$(CXX) $(CXXFLAGS) -o Twiggy Cassandra.o Twiggy.o $(LDFLAGS) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)
	rm -f *~ .depend

include .depend
