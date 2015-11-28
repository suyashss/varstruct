all: varstruct.cpp varstruct.h
	g++ -o varstruct varstruct.cpp

clean:
	rm -f varstruct varstruct.o
