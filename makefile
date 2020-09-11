
FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: cluster
clean:
	rm -rf *.o cluster

cluster: graph.o group.o partition.o computations.o algorithms.o main.o
	gcc $(FLAGS) graph.o group.o partition.o computations.o algorithms.o main.o -o cluster $(LIBS)

graph.o: graph.c graph.h
	gcc $(FLAGS) -c graph.c

group.o: group.c group.h graph.h
	gcc $(FLAGS) -c group.c

partition.o: partition.c partition.h group.h
	gcc $(FLAGS) -c partition.c

computations.o: computations.c computations.h partition.h
	gcc $(FLAGS) -c computations.c

algorithms.o: algorithms.c algorithms.h computations.h
	gcc $(FLAGS) -c algorithms.c

main.o: main.c algorithms.h
	gcc $(FLAGS) -c main.c