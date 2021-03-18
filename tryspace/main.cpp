#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <boost/math/quadrature/gauss.hpp>
using boost::math::quadrature::gauss;
#define xGauss gauss<double,64>::abscissa()
#define wGauss gauss<double,64>::weights()

bool AreSame(double a, double b);


int main() {

    int a = 5;
    // try with square roots:
    double square = sqrt(a);
    double squareDouble = sqrt((double) a);
    auto squareAuto = sqrt(a);

    std::cout << "End of program" << std::endl;
}

bool AreSame(double a, double b)
{
    return fabs(a - b) < std::numeric_limits<double>::epsilon();
}