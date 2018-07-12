CC= g++ -Wall
CFLAGS= -I/usr/gsl/include/
LIBS= -L/third/gsl/lib/



test: main.o 
	$(CC) $(LFLAGS) -g -o test main.o -lgsl -lgslcblas -lm

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp wsbt.cpp 

gtest: gtest.o
	$(CC) -g -o gtest gtest.o -lgtest -lpthread -lgsl -lgslcblas

gtest.o: gtest.cpp
	$(CC) -I/usr/local/include -c gtest.cpp




clean:
	#rm main.o
	#rm wsbt.o
	#rm gtest.o
	rm *.o
	rm out.*
	rm test
	rm gtest
