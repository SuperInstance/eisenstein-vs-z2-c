CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2
AR = ar

SRC = src/eisenstein_vs_z2.c
OBJ = $(SRC:.c=.o)
LIB = libeisenstein_vs_z2.a
HDR = include/eisenstein_vs_z2.h

.PHONY: all lib test clean

all: lib

lib: $(LIB)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

src/%.o: src/%.c $(HDR)
	$(CC) $(CFLAGS) -Iinclude -c -o $@ $<

test: $(LIB) tests/test_eisenstein.c $(HDR)
	$(CC) $(CFLAGS) -Iinclude -o test_eisenstein tests/test_eisenstein.c $(LIB) -lm
	./test_eisenstein

clean:
	rm -f src/*.o $(LIB) test_eisenstein
