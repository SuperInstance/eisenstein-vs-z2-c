CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c99 -O2 -Iinclude
SRCS = src/eisenstein_vs_z2.c
TESTS = tests/test_eisenstein.c

.PHONY: all lib test clean

all: test

lib: $(SRCS)
	$(CC) $(CFLAGS) -c src/eisenstein_vs_z2.c -o build/eisenstein_vs_z2.o
	ar rcs build/libeisenstein_vs_z2.a build/eisenstein_vs_z2.o

test: $(SRCS) $(TESTS)
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRCS) $(TESTS) -lm -o build/test_eisenstein
	./build/test_eisenstein

clean:
	rm -rf build/
