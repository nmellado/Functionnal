/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#include "Functionnal/functionnal.h"
#include "Functionnal/constrainedBezier.h"

namespace Test_Bezier{
typedef double Scalar;
enum{Dim=2};
enum{Degree=3};

typedef functionnal::Bezier<Scalar, Degree, Dim> BezierCurve;
typedef functionnal::BezierMap<Scalar, Degree, Dim> BezierMap;
typedef functionnal::ConstrainedBezierMap<Scalar, Degree+2, Dim> ConstrBezierMap;

}


int main(int /*argc*/, char */*argv*/[])
{
    using namespace Test_Bezier;
    using namespace std;

    const int nbSample = 100;

    ofstream filebuf;

    // build a random bezier curve
    //BezierCurve bezier (BezierCurve::CoeffType::Random());
    BezierCurve bezier ( {
                             0.20, 1.0 ,
                             0.35, 0.25,
                             0.45, 0.25,
                             0.80, 1.0
                         } );

    // create a new bezier curve from the first one, sharing the same memory,
    // with two extra end-points in (0,0) and (1,1)
    ConstrBezierMap cBezier ( bezier.coeffs.data() );
    cBezier.startPoint << 0.0, 0.0;
    cBezier.endPoint   << 1.0, 0.0;


    filebuf.open ("control0.txt");
    filebuf << "# Random control points \n"
            << "# x y " << std::endl;
    filebuf << bezier << std::endl;
    filebuf.close();

    filebuf.open ("control1.txt");
    filebuf << "# Constrained + Random control points \n"
            << "# x y " << std::endl;
    filebuf << cBezier << std::endl;
    filebuf.close();



    filebuf.open ("curves.txt");
    {
        filebuf<<"#x0 y0 x1 y1" << std::endl;
        ConstrBezierMap::EmbeddedVectorType input;

        for (int t = 0; t != nbSample; ++t){
            input << Scalar(t)/Scalar(nbSample-1);

            filebuf << bezier.eval( input ).transpose() << " "
                    << cBezier.eval( input ).transpose() << std::endl;
        }
    }
    filebuf.close();

    // test map to object conversion
    BezierMap bezierMap ( bezier.coeffs.data() );       // changes in bezier will affect bezierMap
    BezierCurve cBezierCopy = bezierMap.asFunctionnal(); // duplicate memory: objects are independants
    std::cout << "Autonomous copy of Constrained control points"<<std::endl;
    std::cout << cBezierCopy << std::endl;

    // compute bezier hodograph
    BezierCurve::Derivative bezierDer = bezier.derivative();
    filebuf.open ("hcontrol0.txt");
    filebuf << "# Hodograph0"<<std::endl;
    filebuf << bezierDer << std::endl;
    filebuf.close();

    // compute constrained bezier hodograph
    ConstrBezierMap::Derivative cBezierDer = cBezier.derivative();
    filebuf.open ("hcontrol1.txt");
    filebuf << "# Hodograph1"<<std::endl;
    filebuf << cBezierDer << std::endl;
    filebuf.close();

    // compute stretch of cBezierCopy using its derivative
    filebuf.open ("hodographs.txt");
    {
        filebuf<<"# Bezier curve stretch" << std::endl;
        ConstrBezierMap::EmbeddedVectorType input;

        for (int t = 0; t != nbSample; ++t){
            input << Scalar(t)/Scalar(nbSample-1);

            filebuf << bezierDer.eval( input ).transpose() << " "
                    << cBezierDer.eval( input ).transpose() << std::endl;
        }
    }
    filebuf.close();

    return EXIT_SUCCESS;
}
