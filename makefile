all: mstruct.cpp mstruct.h
	g++ -o mstruct mstruct.cpp

clean:
	rm mstruct mstruct.o
