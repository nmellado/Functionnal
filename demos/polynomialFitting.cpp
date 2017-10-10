/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>
#include <Eigen/StdVector>

#include "Functionnal/functionnal.h"


template <typename _Base>
struct Functor_w_df : public Eigen::NumericalDiff<_Base, Eigen::NumericalDiffMode::Central>
{
    using Base = _Base;
    using Scalar = typename Base::Scalar;

    inline Functor_w_df(const Base& fun, Scalar _epsfcn)
        : Eigen::NumericalDiff<_Base, Eigen::NumericalDiffMode::Central>(fun, _epsfcn) {}
};

template <typename _FunctionMap>
struct functor : Eigen::DenseFunctor<typename _FunctionMap::Scalar>
{
    using FunctionMap    = _FunctionMap;
    using Scalar         = typename FunctionMap::Scalar;
    using Base           = Eigen::DenseFunctor<Scalar>;
    using InputType      = typename Base::InputType;
    using ValueType      = typename Base::ValueType;
    using JacobianType   = typename Base::JacobianType;
    using QRSolver       = typename Base::QRSolver;


    template <typename Container >
    inline functor(const Container& samples_)
        : Base(FunctionMap::nbCoeff(), samples_.size()),
          samples(samples_)
    { }


    inline functor(const functor &other)
        : Base(FunctionMap::nbCoeff(), other.samples.size()),
          samples(other.samples)
    { }

    inline int operator() (const InputType &x, ValueType &fvec) const {
        FunctionMap q(x.data());

        const int NbSamples = Base::values();

        assert(samples.size() == samplesValues.size());
        assert(samples.size() == NbSamples);

        for (int i = 0; i < NbSamples; ++i) {
            // errors must not be squared here: sign is useful for
            // derivatives computations, and the squared sum is performed
            // by LM main class.
            fvec(i) = q(samples[i]) - samplesValues[i];
        }
        return 0;
    }

//    inline int df(const InputType& _x, JacobianType &jac) const
//    {
//        FunctionMap q(_x.data());
//        for(int d = 0; d < FunctionMap::dim(); ++d) {
//            typename FunctionMap::Derivative der = q.derivative(d);
//            der(_x, jac.row(d));
//        }
//    }

protected:

public:
    const std::vector<typename FunctionMap::InputVectorType,
                      Eigen::aligned_allocator<typename FunctionMap::InputVectorType> >& samples;
    std::vector<Scalar> samplesValues;
};


template <typename Function, typename FunctionMap, typename Container>
void runFit( const Container& samples,
             const Function& targetFunction = Function()) {

    using Point =  typename Function::InputVectorType;
    using Scalar = typename Function::Scalar;

    // Create Eigen data structures
    using FunctorType = functor<FunctionMap>;
    using FunctorTypeWithAutoDiffs = Functor_w_df<FunctorType>;
    FunctorTypeWithAutoDiffs fun(FunctorType (samples), 0.0001);
//    FunctorType fun( samples );

    fun.samplesValues.clear();
    fun.samplesValues.reserve(samples.size());

    // Generate Input values and samples:
    // - random sampling in 10 dimension
    // - value computed by targetFunction
    for (int i = 0; i < samples.size(); ++i)
    {
        fun.samplesValues.push_back( targetFunction(samples[i]) );
    }

    // prepare solver
    Eigen::LevenbergMarquardt<FunctorTypeWithAutoDiffs> lm(fun);
//    Eigen::LevenbergMarquardt<FunctorType> lm(fun);

    // init problem
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> coeffs(Function::nbCoeff(), 1);
    coeffs.setRandom();
    coeffs *= 0.1;
    coeffs += targetFunction.coeffs.transpose();

    //std::cout << "input coeffs: " << coeffs.transpose() << std::endl;
    std::cout << "input error: "
        << std::sqrt((coeffs.transpose() - targetFunction.coeffs).squaredNorm())
        << std::endl;

    // configure solver parameters
    lm.setXtol(1.0e-6); // change between two consecutive function value
    lm.setGtol(1.0e-6); // change between two consecutive gradient norm
    lm.setMaxfev(1000); // max iterations

                        // solve
    lm.minimize(coeffs);

    std::cout << "nb iterations: " << lm.iterations() << std::endl;
    std::cout << "residuals: " << lm.fnorm() << std::endl;

    //std::cout << "computed coeffs: " << coeffs.transpose() << std::endl;
    //std::cout << "targetFunction coeffs: " << targetFunction.coeffs << std::endl;

    std::cout << "fitting error (coeff-wise): "
        << std::sqrt((coeffs.transpose() - targetFunction.coeffs).squaredNorm())
        << std::endl;

}

template <int Dim>
void run(int nbSamples) {
    typedef double Scalar;

    using Quadric = functionnal::Quadric < Scalar, Dim >;
    using QuadricMap = functionnal::ConstQuadricMap<Scalar, Dim>;

    using Linear = functionnal::Linear<Scalar, Dim>;
    using LinearMap = functionnal::ConstLinearMap<Scalar, Dim>;

    using QuadraticForm = functionnal::QuadraticForm < Scalar, Dim>;
    using QuadraticFormMap = functionnal::ConstQuadraticFormMap < Scalar, Dim>;

    std::vector<typename Quadric::InputVectorType,
                Eigen::aligned_allocator<typename Quadric::InputVectorType> > samples;
    samples.resize(nbSamples);

    // Generate Input values and samples:
    // - random sampling in 10 dimension
    // - value computed by targetFunction
    for (int i = 0; i < nbSamples; ++i) {
        samples[i] = Quadric::InputVectorType::Random();
    }

    std::cout << "## Linear fit" << std::endl;
    runFit <Linear, LinearMap>(samples);

    std::cout << "## Quadric fit" << std::endl;
    runFit <Quadric, QuadricMap>(samples);

    std::cout << "## Quadratic Form fit" << std::endl;
    runFit <QuadraticForm, QuadraticFormMap>(samples);
}

/// In this example we fit a quadric function on real-valued points embedded in a
/// 10-dimensionnal space, by minimizing the orthogonal distance between the
/// function and the values v_i at position x_i in the least-squares sense,
/// e.g. for $n$ samples, we find the quadric $argmin_q \sum_i^n (q(x_i) - v_i)^2$.
///
/// The optimization is performed using the Eigen::LevenbergMarquardt solver, and
/// numerical differentiation (for demonstration purpose).
int main(int, char **)
{
    using namespace std;

    std::cout << "# 1D" << std::endl;
    run<2>(1000);
    std::cout << "# 3D" << std::endl;
    run<3>(10000);
    std::cout << "# 7D" << std::endl;
    run<7>(10000);

    return EXIT_SUCCESS;
}
