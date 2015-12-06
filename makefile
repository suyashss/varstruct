varstruct: varstruct.cpp varstruct.h
	g++ -o varstruct varstruct.cpp

clean:
	rm -f varstruct varstruct.o

test:
	python tests/cram-0.4/cram.py -q tests/tiny.t
