/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include <Eigen/Core>

namespace functionnal{

/*!
 *\brief Helper class defining a Functionnal object storing its own memory,
 *
 */
template < typename _EvalBase >
struct Functionnal {

public:
    typedef _EvalBase EvalBase;
    typedef typename EvalBase::Scalar Scalar;

    typedef typename EvalBase::InputVectorType InputVectorType;
    typedef typename EvalBase::OutputVectorType OutputVectorType;
    typedef typename EvalBase::CoeffType CoeffType;
    typedef Functionnal<typename EvalBase::Derivative> Derivative;

protected:
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = EvalBase::NbCoeff
    };

public:
    static constexpr int nbCoeff() { return NbCoeff; }
    static constexpr int dim() { return Dim; }


    //! Actual storage of the coefficients
    CoeffType coeffs;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline
    Functionnal(){ EvalBase::initCoeffs( coeffs); }


    template <typename Derived>
    inline
    Functionnal(const Eigen::MatrixBase<Derived>& otherCoeffs)
          : coeffs(otherCoeffs){
    }   // update with input param


    /*!
     * \brief Functionnal
     * \param coeffArray curve constructors
     *
     * Can be called using brace-initializer lists:
     * \code
     * using MyFonctionnal = Functionnal<SomeEvalFunction>;
     * MyFonctionnal f ( {1.0 , 2., } ); // MyFonctionnal::NbCoeff = 2
     * \endcode
     */
    inline
    Functionnal(const std::array<Scalar, NbCoeff>& coeffArray)
          : coeffs(coeffArray.data()){
    }   // update with input param


    //! \brief Init with default coefficients
    inline
    void initCoeffs()
    {
        EvalBase::initCoeffs(coeffs);
    }

    inline
    OutputVectorType eval(const InputVectorType& x) const{
        return EvalBase::staticEval(x, coeffs);
    }

    /// Convenience function, allowing direct use of brace-initializer lists:
    /// \code
    /// f.eval( {{ x,  y}} );
    /// \endcode
    inline
    OutputVectorType eval(const std::array<Scalar, InputVectorType::SizeAtCompileTime> & x) const{
        return EvalBase::staticEval(Eigen::Map<const InputVectorType>(x.data()), coeffs);
    }

    /*!
     * \brief Build and return de derivative of the Functionnal
     */
    inline
    Derivative derivative() const {
        return Derivative(EvalBase::staticDerivative(coeffs));
    }

public:
    template <class StreamT>
    inline
    void print(StreamT& stream) const{
        return EvalBase::staticPrint(stream, coeffs);
    }
};


namespace internal{
/*!
 *\brief Helper class defining a Functionnal object using external memory
 */
template <typename _EvalBase,
          template <typename> class MapType >
struct FunctionnalMapBase {

public:
    typedef _EvalBase EvalBase;
    typedef typename EvalBase::Scalar Scalar;

    typedef typename EvalBase::InputVectorType InputVectorType;
    typedef typename EvalBase::OutputVectorType OutputVectorType;
    typedef typename EvalBase::CoeffType CoeffType;
    typedef Functionnal<typename EvalBase::Derivative> Derivative;

protected:
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = EvalBase::NbCoeff
    };

public:
    static constexpr int nbCoeff() { return NbCoeff; }
    static constexpr int dim() { return Dim; }

    //! Map to the actual coefficients
    MapType<CoeffType> coeffs;


    inline FunctionnalMapBase(Scalar* data)
        :coeffs(data){}
    inline FunctionnalMapBase(const Scalar* data)
        :coeffs(data){}

    inline Functionnal<EvalBase> asFunctionnal() const{
        return Functionnal<EvalBase>(coeffs);
    }

    //! \brief Init with default coefficients
    inline
    void initCoeffs()
    { EvalBase::initCoeffs(coeffs); }

    inline
    OutputVectorType eval(const InputVectorType& x) const{
        return EvalBase::staticEval(x, coeffs);
    }

    //! Build and return de derivative of the FunctionnalMap as Functionnal
    inline
    Derivative derivative() const {
        return EvalBase::derivative(coeffs);
    }

public:
    template <class StreamT>
    inline
    void print(StreamT& stream) const{
        return EvalBase::staticPrint(stream, coeffs);
    }
};

template<class T> using Map      = Eigen::Map<T>;
template<class T> using ConstMap = const Eigen::Map<const T>;

} // namespace internal
} // namespace Functionnal

////////////////////////////////////////////////////////////////////////////////
// Predefined types

#include "Functionnal/functionnalEvalFunc.h"


namespace functionnal{

template<typename _EvalBase>
using FunctionnalMap =
internal::FunctionnalMapBase<_EvalBase, internal::Map >;

template<typename _EvalBase>
using ConstFunctionnalMap =
internal::FunctionnalMapBase<_EvalBase, internal::ConstMap >;



// Polynomial typedefs
template<typename _Scalar, int _Dim>
using Quadric =
Functionnal< GeneralizedQuadricEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using QuadricMap =
FunctionnalMap< GeneralizedQuadricEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using ConstQuadricMap =
ConstFunctionnalMap< GeneralizedQuadricEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using Linear =
Functionnal< LinearEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using LinearMap =
FunctionnalMap< LinearEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using ConstLinearMap =
ConstFunctionnalMap< LinearEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using QuadraticForm =
Functionnal< QuadraticFormEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using QuadraticFormMap =
FunctionnalMap< QuadraticFormEvalFunc< _Scalar, _Dim > >;

template<typename _Scalar, int _Dim>
using ConstQuadraticFormMap =
ConstFunctionnalMap< QuadraticFormEvalFunc< _Scalar, _Dim > >;



// Polynomial typedefs
template<typename _Scalar, int _Degree, int _Dim>
using CWisePolynomial =
Functionnal< CWisePolynomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using CWisePolynomialMap =
FunctionnalMap< CWisePolynomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using CWiseConstPolynomialMap =
ConstFunctionnalMap< CWisePolynomialEvalFunc< _Scalar, _Degree, _Dim > >;


// Bezier typedefs
template<typename _Scalar, int _Degree, int _Dim>
using Bezier =
Functionnal< BezierBinomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using BezierMap =
FunctionnalMap< BezierBinomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using ConstBezierMap =
ConstFunctionnalMap< BezierBinomialEvalFunc< _Scalar, _Degree, _Dim > >;

} // namespace Functionnal

#include "Functionnal/functionnal.hpp"

