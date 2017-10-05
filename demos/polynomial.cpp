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
enum{NbElemt=50}; // we test on 50 samples at once
enum{Dim = 2};
enum{Degree=3};

using Curve = functionnal::CWisePolynomial<Scalar, Degree, NbElemt>;
using Quadric1D = functionnal::Quadric<Scalar,  Dim>;

}


int main(int /*argc*/, char */*argv*/[])
{
    using namespace Test_Polynomial;
    using namespace std;

    Curve c1 ( { 0, 0, 0, 1 } );             // x3
    Curve c2 ( { 0.5, 0, -2, 0.7 } );        // 0.5 - 2x^2 + 0.7x^3
    Curve c3 ( Curve::CoeffType::Random() ); //

    typename Curve::InputVectorType  input  = Curve::InputVectorType::Random();
    typename Curve::OutputVectorType o1 = c1.eval( input );
    typename Curve::OutputVectorType o2 = c2.eval( input );
    typename Curve::OutputVectorType o3 = c3.eval( input );

    std::cout << "# \"f(x) = " << c1 << "\""
              << "  \"g(x) = " << c2 << "\""
              << "  \"h(x) = " << c3 << "\"" << std::endl;

    for (int i = 0; i < NbElemt; ++i)
        std::cout << input(i) << " "
                  << o1(i) << " "
                  << o2(i) << " "
                  << o3(i) << "\n";
    std::cout << std::flush;

    std::cout << "Testing 2D quadric f(x,y) = x^2 + y^2" << std::endl;
    Quadric1D q1;
    Quadric1D::EvalBase::getQMap( q1.coeffs ) = Quadric1D::EvalBase::QType::Identity();
    Quadric1D::EvalBase::getPMap( q1.coeffs ) = Quadric1D::EvalBase::PType::Zero();

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
