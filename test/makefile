# makefile for twoPointTest
# 
# Catherine J. Stevenson
# December 10, 2007
# 
Sourcecode = /Users/catherine/Documents/School/SourceCode/
include $(Sourcecode)/Makefile.machine

#CXX = g++
CXXFLAGS = $(warn_flags) $(include_flags) \
	-I../ \
	-I$(Sourcecode)/diagon/include \
	-I$(Sourcecode)/Utility


sources := $(wildcard *.cpp)
objects := $(sources:.cpp=.o)
depends := $(sources:.cpp=.d)
diagonlib = $(Sourcecode)/diagon/lib/

otherObjects = ../twoPoint.o

twoPointTest: $(objects) $(diagonlib)/libdiagon.a $(otherObjects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(objects) $(otherObjects) \
	-L$(diagonlib) -ldiagon $(testlib)

check: twoPointTest
	./twoPointTest < twoPointUnitTest.in

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
	rm -f $(objects) $(depends) twoPointTest a.out *~
