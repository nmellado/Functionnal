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

#include "Functionnal/functionnal.h"


template <typename _Base>
struct Functor_w_df : public Eigen::NumericalDiff<_Base, Eigen::NumericalDiffMode::Central>
{
    using Base = _Base;
    using Scalar = typename Base::Scalar;

    inline Functor_w_df(Scalar _epsfcn)
        : Eigen::NumericalDiff<_Base, Eigen::NumericalDiffMode::Central>(_epsfcn) {}
};

template <int NbSamples, typename _FunctionMap>
struct functor : Eigen::DenseFunctor<typename _FunctionMap::Scalar>
{
    using FunctionMap    = _FunctionMap;
    using Scalar         = typename FunctionMap::Scalar;
    using Base           = Eigen::DenseFunctor<Scalar>;
    using InputType      = typename Base::InputType;
    using ValueType      = typename Base::ValueType;
    using JacobianSolver = typename Base::JacobianType;
    using QRSolver       = typename Base::QRSolver;


    inline functor() : Base(FunctionMap::nbCoeff(), NbSamples) { }

    inline int operator() (const InputType &x, ValueType &fvec) const {
        FunctionMap q(x.data());

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

public:
    std::vector<typename FunctionMap::InputVectorType> samples;
    std::vector<Scalar> samplesValues;
};


template <int nbSamples, typename Function, typename FunctionMap>
void runFit( const std::vector<typename Function::InputVectorType>& samples,
          const Function& targetFunction = Function()) {

    using Point =  typename Function::InputVectorType;
    using Scalar = typename Function::Scalar;

    // Create Eigen data structures
    using FunctorType = Functor_w_df<functor<nbSamples, FunctionMap>>;
    FunctorType fun(std::numeric_limits<Scalar>::epsilon());

    fun.samples = samples;

    // Generate Input values and samples:
    // - random sampling in 10 dimension
    // - value computed by targetFunction
    for (const auto& p : fun.samples)
    {
        fun.samplesValues.push_back( targetFunction(p) );
    }

    // prepare solver
    Eigen::LevenbergMarquardt<FunctorType> lm(fun);

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

template <int Dim, int nbSamples>
void run() {
    typedef double Scalar;

    using Quadric = functionnal::Quadric < Scalar, Dim >;
    using QuadricMap = functionnal::ConstQuadricMap<Scalar, Dim>;

    using Linear = functionnal::Linear<Scalar, Dim>;
    using LinearMap = functionnal::ConstLinearMap<Scalar, Dim>;

    using QuadraticForm = functionnal::QuadraticForm < Scalar, Dim>;
    using QuadraticFormMap = functionnal::ConstQuadraticFormMap < Scalar, Dim>;

    std::vector<typename Quadric::InputVectorType> samples;
    samples.resize(nbSamples);

    // Generate Input values and samples:
    // - random sampling in 10 dimension
    // - value computed by targetFunction
    for (int i = 0; i < nbSamples; ++i) {
        samples[i] = Quadric::InputVectorType::Random();
    }

    std::cout << "## Linear fit" << std::endl;
    runFit <nbSamples, Linear, LinearMap>(samples);

    std::cout << "## Quadric fit" << std::endl;
    runFit <nbSamples, Quadric, QuadricMap>(samples);

    std::cout << "## Quadratic Form fit" << std::endl;
    runFit <nbSamples, QuadraticForm, QuadraticFormMap>(samples);
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

    std::cout << "# 2D" << std::endl;
    run<2, 10000>();
    std::cout << "# 3D" << std::endl;
    run<3, 1000>();
    std::cout << "# 10D" << std::endl;
    run<10, 1000>();

    return EXIT_SUCCESS;
}
