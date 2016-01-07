bin/varstruct: src/varstruct.cpp src/varstruct.h
	mkdir -p bin
	g++ -o bin/varstruct src/varstruct.cpp

clean:
	rm -f bin/varstruct

test:
	python tests/cram-0.4/cram.py -q tests/tiny.t
