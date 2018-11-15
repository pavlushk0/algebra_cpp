
CPP      := g++
CC		 := gcc
CPPFLAGS := -Wall -std=c++17 
CFLAGS	 := -Wall -std=c11
CPPLIBS  := -lstdc++ -lSDL2 -lGL -lGLU
CLIBS    :=	-lc -lm
RM       := rm -f

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

$(info COMMON MAKEFILE)

dummy:
	$(info no make target)

cpp_version:
	$(info CPP_VERSION file compile)
	$(CPP) $(CPPFLAGS) test.cpp vec3.cpp mtrx3.cpp -o main $(CPPLIBS)

c_version:
	$(info C_VERSION file compile)
	$(CC) $(CFLAGS) test.c vec3.c mtrx3.c mtrx4.c mtrxC.c qtnn.c -o main $(CLIBS)

clean:
	$(RM) *.o
	$(RM) main
