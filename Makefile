CC=g++
GPP=g++

DISABLEWARN = -Wno-unused-but-set-variable -Wno-unused-result -Wno-unused-variable
#ARCHFLAGS =  -march=x86-64  -mno-avx
#CXXFLAGS +=  -O3 $(ARCHFLAGS)  -Wall -static -I.  -std=c++11

OIGRAPH = -I/usr/local/include/igraph/ -L/usr/local/lib -ligraph -DUSE_IGRAPH

CXXFLAGS +=  -O3 $(ARCHFLAGS)  -Wall -I.  -std=c++11 

LIBS = -lpthread

CFLAGS = -O3

DEPS = *.cpp *.h
OBJS=argtable3.o 
.PHONY:	all clean

PROGS= gclu gclu_ig

all: gclu

#Argtable should support compiling with g++, but there was an error message.
argtable3.o:
	gcc -c $(CFLAGS) contrib/argtable3.c

options.o:
	$(CC) -c $(CXXFLAGS) options.c

gclu: $(DEPS) $(OBJS)
	$(CC) $(CXXFLAGS) $(DISABLEWARN) graphclu.cpp $(LIBS) $(OBJS) -o gclu  -static

gclu_gdb: $(DEPS) $(OBJS)
	$(CC) -g -O1 $(ARCHFLAGS)  -Wall -I.  -std=c++11 $(DISABLEWARN) graphclu.cpp $(LIBS) $(OBJS) -o gclu_gdb  -static


gclu_ig: $(DEPS) $(OBJS)
	$(CC) $(CXXFLAGS) $(DISABLEWARN) graphclu.cpp $(LIBS) $(OBJS) -o gclu  -static $(OIGRAPH)

clean:
	rm -f $(PROGS) *.o

