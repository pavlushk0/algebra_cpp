
CPP      := g++
CC		 := gcc
CPPFLAGS := -Wall -std=c++17 
CFLAGS	 := -Wall -std=c11
CPPLIBS  := -lstdc++
CLIBS    :=	-lc -lm
RM       := rm -f

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

$(info COMMON MAKEFILE)

dummy:
	$(info no make target)

cpp_version:
	$(info CPP_VERSION file compile)
	$(CPP) $(CPPFLAGS) test.cpp vec3.cpp mtrx3.cpp mtrx4.cpp -o main $(CPPLIBS)

clean:
	$(RM) *.o
	$(RM) main
