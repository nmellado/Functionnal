# Functionnal lib
Simple header-only set of classes for efficient Polynomials and Bezier curves evaluations, in arbitrary dimensions.
This library is made to stay simple and easy to maintain. 

## Compilation
Use Cmake to build the sources. For the root folder:

    mkdir build
    cd build
    cmake ..
    make

## Usage
Checkout the `demos/` folder for concrete examples.

### Types of curves
All the curves are implemented using the struct `template < typename EvalBase > struct Functionnal`, and provides the following API:

``` c++
OutputVectorType eval(const InputVectorType& x) const;  // evaluate the input function f(x)
Derivative derivative() const;                          // construct a new object f'(x)
operator << (Stream s) const;                           // print to input stream s: f(x) = 3x^2 + 1
```
This API is implemented by the `EvalBase` template class.

In addition to `Functionnal`, we also provides the classes `FunctionnalMap` and `ConstFunctionnalMap` to use the same API with objects storing the functionnal coefficients as `Eigen::Map`. This is very useful if you need to use the library with raw buffers for storage and synchronisation.
Both `FunctionnalMap` and `ConstFunctionnalMap` constructors take a pointer to the coefficient memory as input. They also provide with a `asFunctionnal()` method to create duplicate objects.
For instance:
```c++
BezierCurve bezier ( {
                        0.20, 1.0 ,
                        0.35, 0.25,
                        0.45, 0.25,
                        0.80, 1.0
                      } );
BezierMap bezierMap ( bezier.coeffs.data() );       // changes in bezier will affect bezierMap
BezierCurve cBezierCopy = bezierMap.asFunctionnal(); // duplicate memory: objects are independants
```

#### Bezier
Example with cubic bezier curves in 3D

``` c++
typedef double Scalar;
enum{Dim=2};
enum{Degree=3};
typedef functionnal::Bezier<Scalar, Degree, Dim> BezierCurve;
BezierCurve bezier ( {
                        0.20, 1.0 ,
                        0.35, 0.25,
                        0.45, 0.25,
                        0.80, 1.0
                      } );
const int nbSample = 100;
for (int t = 0; t != nbSample; ++t){
    BezierCurve::InputVectorType input;
    input << Scalar(t)/Scalar(nbSample-1);
    std::cout << bezier.eval( input ).transpose() << std::endl;
}

```

#### N-Dimensionnal Quadric
Defined as in https://en.wikipedia.org/wiki/Quadric. Example of use in 2D:

``` c++
enum{Dim = 2};
using Scalar  = float;
using Quadric = functionnal::Quadric<Scalar,  Dim>;

Quadric q;

// set f(x,y) = x^2 + y^2
Quadric::EvalBase::getQMap( q.coeffs ) = Quadric::EvalBase::QType::Identity();
Quadric::EvalBase::getPMap( q.coeffs ) = Quadric::EvalBase::PType::Zero();

// evaluate quadric:
Scalar res0 = q.eval( {{ 1., 1. }} );
Scalar res1 = q.eval( Quadric::InputVectorType::Random() );
```

#### Coefficient-wise polynomial of arbitrary degree
Evaluate a 1D polynomial function of arbitrary degree on a set of scalar values.
```c++
using Curve = functionnal::CWisePolynomial<Scalar, Degree, NbElemt>;
Curve c2 ( { 0.5, 0, -2, 0.7 } );        // 0.5 - 2x^2 + 0.7x^3

typename Curve::InputVectorType  input  = Curve::InputVectorType::Random();

// Process all the input samples at once
typename Curve::OutputVectorType o2 = c2.eval( input );

for (int i = 0; i < NbElemt; ++i)
        std::cout << input(i) << " " << o2(i) << "\n";
    std::cout << std::flush;

```
    
