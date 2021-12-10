CXXFLAG = -std=c++20
CXX = g++
BIN = P2D

all: ${BIN}

${BIN} : main.o ${BIN}.o Progonka.o
	${CXX} ${BIN}.o main.o Progonka.o  -o ${BIN}

main.o: main.cpp ${BIN}.hpp Progonka.hpp
	${CXX} -c main.cpp -o main.o

${BIN}.o : ${BIN}.hpp ${BIN}.cpp Progonka.hpp
	${CXX} -c ${BIN}.cpp -o	${BIN}.o

Progonka.o : Progonka.hpp Progonka.cpp
	${CXX} -c Progonka.cpp -o	Progonka.o



plot_ans: ${BIN}
	./${BIN} > ans.txt
	gnuplot plot_ans.gpi
	rm ans.txt

clean:
	rm *.o ${BIN}

