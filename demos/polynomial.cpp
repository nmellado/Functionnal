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

using Quadric = functionnal::Quadric<Scalar,  2>;
using Linear  = functionnal::Linear <Scalar,  10>;
}


int main(int /*argc*/, char */*argv*/[])
{
    using namespace Test_Polynomial;
    using namespace std;


    std::cout << "Testing 10D linear function f(x1, ...,x10) = \\sum_{i=1}^10(x_i)" << std::endl;
    Linear l1;

    // sample linear function in [-1:1], 100 samples
    for( int i = 0; i != 100; ++i) {
        typename Linear::InputVectorType v = Linear::InputVectorType::Random();

        std::cout << "f( " << v.transpose() << " ) = " << l1.eval (v) << "\n";
    }
    std::cout << std::flush;

    std::cout << "Testing 2D quadric f(x,y) = x^2 + y^2" << std::endl;
    Quadric q1;
    Quadric::EvalBase::getQMap( q1.coeffs ) = Quadric::EvalBase::QType::Identity();
    Quadric::EvalBase::getPMap( q1.coeffs ) = Quadric::EvalBase::PType::Zero();

    // sample quadric in [-1:1], 20 samples per dimension (not ordered)
    for( int i = 0; i != 10; ++i) {
        Scalar x = Scalar(i+1)/Scalar(10);
        for( int j = 0; j != 10; ++j) {
            Scalar y = Scalar(j+1)/Scalar(10);
            std::cout <<  x << " " <<  y << " " << q1.eval( {{ x,  y}} ) << "\n";
            std::cout <<  x << " " << -y << " " << q1.eval( {{ x, -y}} ) << "\n";
            std::cout << -x << " " <<  y << " " << q1.eval( {{-x,  y}} ) << "\n";
            std::cout << -x << " " << -y << " " << q1.eval( {{-x, -y}} ) << "\n";
        }
    }
    std::cout << std::flush;

    return EXIT_SUCCESS;
}
