CC=gcc

all: final clean

final: fft.o 2dfft.o img_write.o rotate.o
	$(CC) fft.o 2dfft.o main.c img_write.o rotate.o -lm -o find_char

clean:
	rm -rf *.o
