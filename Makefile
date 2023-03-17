#make

.PHONY:all
all: main

main: main.o Matrix.o
	g++ main.o Matrix.o -o main
main.o: Matrix.hpp main.cpp
	g++ -Wall -c main.cpp
Matrix.o: Matrix.hpp Matrix.cpp
	g++ -Wall -c Matrix.cpp

.PHONY:clean
clean:
	rm main *.o

.PHONY:tar
tar:
	tar zcvf Matrix.tar.gz main.cpp Matrix.cpp Matrix.hpp Makefile
