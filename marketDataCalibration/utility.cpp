//
// Created by Marco Paggiaro on 04/02/2021.
//

#include "utility.h"
#include <limits>
#include <cmath>

bool AreSame(double a, double b)
{
    return fabs(a - b) < std::numeric_limits<double>::epsilon();
}

std::vector<double> interpolate (const std::vector<double> &x,
                                 const std::vector<double> &y, const std::vector<double> &xH)
{
    std::vector<double> yH(xH.size());
    for (int i = 0; i < xH.size(); ++i)
    {
        for (int j = 1; j < x.size(); ++j) {
            if  (xH[i] <= x[j])
            {
                double slope = ( y[j] - y[j-1]) / ( x[j] - x[j-1] );
                yH[i] = y[j-1] + slope * (xH[i] - x[j-1]);
                break;
            }
        }
    }
    return yH;
}