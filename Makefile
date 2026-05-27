CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude
SRC = src/eisenstein_vs_z2.c
LIB = libeisenstein.a

.PHONY: all lib test clean

all: lib

build:
	mkdir -p build

lib: $(LIB)

$(LIB): $(SRC) include/eisenstein_vs_z2.h | build
	$(CC) $(CFLAGS) -c $(SRC) -o build/eisenstein_vs_z2.o -lm
	ar rcs $@ build/eisenstein_vs_z2.o

test: test_runner
	./test_runner

test_runner: tests/test_eisenstein.c $(LIB) include/eisenstein_vs_z2.h
	$(CC) $(CFLAGS) tests/test_eisenstein.c -L. -leisenstein -lm -o test_runner

clean:
	rm -f $(LIB) test_runner build/*.o
