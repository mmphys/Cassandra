# CERN root library config utility
APFEL     := $(HOME)/bin/apfel-3.0.2
RC        := root-config
#CXX      := g++
#CXX       := g++-mp-8
CXX       := clang++

#CXXFLAGS := -std=c++17 -I/usr/local/include -I$(ROOTSYS)/include
CXXFLAGS  := $(shell $(RC) --cflags)
#CXXFLAGS  := $(subst -stdlib=libc++,,$(CXXFLAGS))
#CXXFLAGS  := $(subst -m64,,$(CXXFLAGS))
#CXXFLAGS  := $(subst -pthread,,$(CXXFLAGS))
#CXXFLAGS  := $(subst -std=c++1z,-std=c++11,$(CXXFLAGS))
CXXFLAGS  := $(CXXFLAGS) -Wall -Wextra -Werror
CXXFLAGS  := $(CXXFLAGS) -I/usr/local/include -I$(APFEL)/include

LDDIRS    := $(shell lhapdf-config --libdir) $(APFEL)/lib $(ROOTSYS)
LDFLAGS   := $(addprefix -L,$(sort $(LDDIRS)))

LDLIBS    := $(shell $(RC) --libs)
#LDLIBS    := $(subst -stdlib=libc++,,$(LDLIBS))
#LDLIBS    := $(subst -m64,,$(LDLIBS))
#LDLIBS    := $(subst -pthread,,$(LDLIBS))
#LDLIBS    := $(subst -std=c++1z,-std=c++17,$(LDLIBS))
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
