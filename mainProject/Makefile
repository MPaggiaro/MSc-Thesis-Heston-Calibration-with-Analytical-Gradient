# OLD version:

# all : main

# main : main.cpp
# 	g++ -Wall main.cpp -o main

# clean:
# 	rm main

objects = main.o Calibration.o EuropeanOption.o Market.o

all : main

main : $(objects)
	g++ $(objects) -o main -llevmar -llapack -lblas -lm

#main.o : main.cpp
#	g++ -Wall -c main.cpp
#
#EuropeanOption.o: EuropeanOption.cpp
#	g++ -Wall -c EuropeanOption.cpp

%.o: %.cpp
	g++ -Wall -c $<

.PHONY : clean
clean :
	rm main $(objects)