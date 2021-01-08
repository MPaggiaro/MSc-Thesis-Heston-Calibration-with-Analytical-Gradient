# OLD version:

# all : main

# main : main.cpp
# 	g++ -Wall main.cpp -o main

# clean:
# 	rm main

objects = main.o EuropeanOption.o

all : main

main : $(objects)
	g++ $(objects) -o main

main.o : main.cpp
	g++ -Wall -c main.cpp

EuropeanOption.o: EuropeanOption.cpp
	g++ -Wall -c EuropeanOption.cpp

.PHONY : clean
clean :
	rm main $(objects)