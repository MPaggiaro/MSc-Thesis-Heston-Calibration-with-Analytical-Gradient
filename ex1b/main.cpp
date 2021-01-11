#include <iostream>
#include <ctime>
#include <vector>
#include "EuropeanOption.h"
#include "Market.h"

template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (auto ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}



int main() 
{
    Market market(0.02,0,3.0,0.1,0.25,-0.8,0.08,1);
    EuropeanOption::market = market;
    EuropeanOption myOption(1.1,1,"C");
    double start = std::clock();

    std::cout.precision(12);
    std::cout << "Option (" << myOption.cfType << ") price = " << myOption.Price() << std::endl;
    double stop = std::clock();

    std::cout << "Jacobian: " << myOption.ComputeJacobian() << std::endl;

    std::cout << "Elapsed time: " << double(stop - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
    return 0;
}