/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#include <iostream>
//#include <fstream>
#include <Eigen/Dense>

#include "Functionnal/functionnal.h"

namespace Test_Polynomial{
typedef double Scalar;
enum{Dim=50}; // we test on 50 samples at once
enum{Degree=3};

using Curve = Functionnal::CWisePolynomial<Scalar, Degree, Dim>;

}


int main(int /*argc*/, char */*argv*/[])
{
    using namespace Test_Polynomial;
    using namespace std;

    // c1(x) = 4.5x**3 + 3.5x**2 + 2.5x + 1.5
    //Curve c1 ( {  4.5, 3.5, 2.5, 1.5  });
    Curve c1 /*( { 0, 0, 0, 1  })*/;

    auto input  = Curve::EmbeddedVectorType::Random();
    auto output = c1.eval( input );

    for (int i = 0; i < Dim; ++i)
        std::cout << input(i) << " " << output(i) << "\n";
    std::cout << std::flush;

    return EXIT_SUCCESS;
}
