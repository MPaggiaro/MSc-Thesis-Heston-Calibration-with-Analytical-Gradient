//
// Created by Marco Paggiaro on 04/02/2021.
//

#include "doublecomparison.h"
#include <limits>
#include <cmath>

bool AreSame(double a, double b)
{
    return fabs(a - b) < std::numeric_limits<double>::epsilon();
}
