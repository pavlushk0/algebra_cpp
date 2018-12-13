
CC       := g++
FLAGS    := -Wall -std=c++17
LIBS 	 := -lstdc++
RM       := rm -f

SRCS := test.cpp mtrx234.cpp vec234.cpp qtnn.cpp plane.cpp
OBJS := $(SRCS:.cpp=.o)

$(info COMMON MAKEFILE)

all: compile link 

compile:
	$(CC) $(FLAGS) -c $(SRCS)

link:
	$(CC) $(OBJS) -o main $(LIBS)

clean:
	$(RM) *.o
	$(RM) main
