TARGET = prog.out
PROGRAM_NAME = prog.out
LIBS = -lm
CC = gcc
CFLAGS = -g -Wall

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) -o mmio.o -c mmio.c
	ar rc libmmio.a mmio.o
	ranlib libmmio.a
	mpicc mult_vector.c -o prog.out -L. -lmmio

run:
	mpirun ${PROGRAM_NAME} ${ARGS}

clean:
	-rm -f *.o
	-rm -f $(TARGET)

