# CERN root library config utility
APFEL     ?= $(HOME)/.localMSc
RC        := root-config
CXX       ?= g++
#CXX       ?= clang++

CXXFLAGS  := $(shell $(RC) --cflags)
CXXFLAGS  := $(filter -I%, $(CXXFLAGS))
CXXFLAGS  := -std=c++17 $(CXXFLAGS) -I/usr/local/include -I$(shell lhapdf-config --incdir)

RPATH     := -Xlinker -rpath -Xlinker $(GridPkg)/lib/gcc13
LDFLAGS    := $(shell $(RC) --libs)
#LDFLAGS    := $(filter -L%, $(LDFLAGS))
LDFLAGS   := $(addprefix -L,$(sort $(shell lhapdf-config --libdir))) $(LDFLAGS) $(RPATH)

#LDLIBS     := $(filter -l%, $(shell $(RC) --libs))
LDLIBS    := $(addprefix -l,APFEL LHAPDF) $(LDLIBS)

srcfiles  := $(shell find . -type f -iname "*.cpp")
objects   := $(patsubst %.cpp, %.o, $(srcfiles))

.PHONY: all
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
