CXXFLAG = -std=c++17
CXX = g++

all: EE

EE : main.o EllipticEquation.o iters.o
	g++ EllipticEquation.o main.o iters.o  -o EE

main.o: main.cpp iters.hpp EllipticEquation.hpp
	g++ -c main.cpp 

EllipticEquation.o: EllipticEquation.cpp EllipticEquation.hpp iters.hpp
	g++ -c EllipticEquation.cpp 	

iters.o : iters.cpp
	g++ -c iters.cpp


plot_ans: EE
	./EE > ans.txt
	gnuplot plot_ans.gpi
	rm ans.txt

clean:
	rm *.o EE

