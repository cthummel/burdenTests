CC= g++ -Wall -std=c++11
CFLAGS= -I/usr/gsl/include/
LIBS= -L/third/gsl/lib/



test: main.o 
	$(CC) $(LFLAGS) -g -o test main.o -lgsl -lgslcblas -lm

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp wsbt.cpp genericBurdenTest.cpp

gtest: gtest.o
	$(CC) -g -o gtest gtest.o -lgtest -lpthread -lgsl -lgslcblas

gtest.o: gtest.cpp
	$(CC) -I/usr/local/include -c gtest.cpp




clean:
	rm *.o out.* *.stats *.txt anno.* test gtest
