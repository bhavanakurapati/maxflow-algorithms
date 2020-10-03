all: MaxFlow

MaxFlow: MaxFlow.cpp
	g++ -o $@ $^ 

clean:
	rm -f *.o MaxFlow
