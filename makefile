CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic
BIN=test_sage
OBJS = test_sage.o sage.o isi_sage.o

test_sage:$(OBJS) 
	$(CXX) -o $(BIN) $(OBJS) $(CXXFLAGS)

isi_sage.o:isi_sage.h
sage.o:sage.h

.PHONY:clean
clean:
	-rm $(BIN) $(OBJS)
