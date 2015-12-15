bin/varstruct: varstruct.cpp varstruct.h
	mkdir -p bin
	g++ -o varstruct varstruct.cpp
	mv varstruct bin

clean:
	rm -f bin/varstruct varstruct.o

test:
	python tests/cram-0.4/cram.py -q tests/tiny.t
