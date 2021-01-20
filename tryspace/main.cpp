#include <iostream>
#include <boost/math/quadrature/gauss.hpp>
using boost::math::quadrature::gauss;
#define xGauss gauss<double,64>::abscissa()
#define wGauss gauss<double,64>::weights()


int main() {

    // right weights and abscissa.
    auto x = xGauss;
    auto w = wGauss;

    // Let's see how they are transformed in the code of Cui.
    const double lb = 0.0, // our zero
    ub = 200, // our N
    Q = 0.5*(ub - lb), // N/2
    P = 0.5*(ub + lb); // N/2

    std::vector<double> hey(10);
    // OK: default-initialized with zeros!

    std::cout << "End of program" << std::endl;
}
