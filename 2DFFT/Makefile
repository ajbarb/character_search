CC=gcc

all: final clean

final: 2dfft.o fft.o
	$(CC) fft.o 2dfft.c -lm -o 2dfft

clean:
	rm -rf *.o
