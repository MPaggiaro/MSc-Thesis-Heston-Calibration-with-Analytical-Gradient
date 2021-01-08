all : hc



hc : HestonCalibrator.cpp
	g++ -Wall HestonCalibrator.cpp -o hc -llevmar -llapack -lblas -lm

clean:
	rm hc