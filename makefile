#Makefile for twoPointSpin.cpp

# This include sets include_flags and libs for
# different machines.
Sourcecode = /Users/catherine/Documents/School/SourceCode/
include $(Sourcecode)/Makefile.machine

#CXX = g++
CXXFLAGS = $(warn_flags) $(include_flags) \
	-I$(Sourcecode)/diagon/include \
	-I$(Sourcecode)/Utility


sources := $(wildcard *.cpp)
objects := $(sources:.cpp=.o)
depends := $(sources:.cpp=.d)
diagonlib = $(Sourcecode)/diagon/lib/

twoPointSpin: $(objects) $(diagonlib)/libdiagon.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(objects) -L$(diagonlib) -ldiagon $(threadlib)

check: twoPointSpin
	./twoPointSpin ../input/twoPointUnitTest.in > ../stdout/Test.stdout 2>&1 &

ifneq ($(MAKECMDGOALS), clean)
-include $(depends)
endif

# Automatically generate prerequisites.  (See "info make" for details)
%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) $(CXXFLAGS) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	rm -f $(objects) $(depends) twoPointSpin Test.out Test.stdout a.out *.d.* *~
