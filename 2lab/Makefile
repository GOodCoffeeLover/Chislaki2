all: H1D

H1D : main.o H1D.o 
	g++ H1D.o main.o  -o H1D

main.o: main.cpp Progonka.hpp H1D.h
	g++ -c main.cpp 

H1D.o: H1D.cpp H1D.h Progonka.hpp
	g++ -c H1D.cpp 	

plot_ans: H1D
	./H1D > ans.txt
	./plot_ans.gpi
	rm ans.txt

clean:
	rm H1D.o main.o

