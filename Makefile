all:
	g++ -O3 utilities/src/*.cpp main.cpp -std=c++11 -o main -pthread
clean:
	rm main
