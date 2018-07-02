CC= g++
CFLAGS= -I./third/gsl-2.4/


test: main.o
	$(CC) -g -o test main.o

main.o: main.cpp wsbt.cpp
	$(CC) -c -g main.cpp $(CFLAGS)

clean:
	rm main.o
