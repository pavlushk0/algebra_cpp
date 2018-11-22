
CC       := g++
FLAGS    := -Wall -std=c++17 
LIBS 	 := -lstdc++
RM       := rm -f

SRCS := mtrx34.cpp test.cpp vec3.cpp
OBJS := $(SRCS:.cpp=.o)

$(info COMMON MAKEFILE)

all: compile link 

compile:
	$(CC) $(FLAGS) -c mtrx34.cpp test.cpp vec3.cpp

link:
	$(CC) $(OBJS) -o main $(LIBS)

clean:
	$(RM) *.o
	$(RM) main
