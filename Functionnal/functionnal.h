/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include <Eigen/Core>

namespace Functionnal{

/*!
 *\brief Helper class defining a Functionnal object storing its own memory
 */
template < typename _EvalBase >
struct Functionnal {
    typedef _EvalBase EvalBase;
    typedef typename EvalBase::Scalar Scalar;
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = EvalBase::NbCoeff
    };
    typedef typename EvalBase::EmbeddedVectorType EmbeddedVectorType;
    typedef typename EvalBase::EmbeddingVectorType EmbeddingVectorType;
    typedef typename EvalBase::CoeffType CoeffType;
    typedef Functionnal<typename EvalBase::Derivative> Derivative;

    //! Actual storage of the coefficients
    CoeffType coeffs;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline
    Functionnal(){ }


    template <typename Derived>
    inline
    Functionnal(const Eigen::MatrixBase<Derived>& otherCoeffs)
          : coeffs(otherCoeffs){
    }   // update with input param


    inline
    Functionnal(const std::array<Scalar, Dim*NbCoeff>& coeffArray)
          : coeffs(coeffArray.data()){
    }   // update with input param


    //! \brief Init with default coefficients
    inline
    void initCoeffs()
    {
        EvalBase::initCoeffs(coeffs);
    }

    inline
    EmbeddingVectorType eval(const EmbeddedVectorType& x) const{
        return EvalBase::staticEval(x, coeffs);
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

    typedef _EvalBase EvalBase;
    typedef typename EvalBase::Scalar Scalar;
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = EvalBase::NbCoeff
    };
    typedef typename EvalBase::EmbeddedVectorType EmbeddedVectorType;
    typedef typename EvalBase::EmbeddingVectorType EmbeddingVectorType;
    typedef typename EvalBase::CoeffType CoeffType;
    typedef Functionnal<typename EvalBase::Derivative> Derivative;

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
    EmbeddingVectorType eval(const EmbeddedVectorType& x) const{
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


namespace Functionnal{

template<typename _EvalBase>
using FunctionnalMap =
internal::FunctionnalMapBase<_EvalBase, internal::Map >;

template<typename _EvalBase>
using ConstFunctionnalMap =
internal::FunctionnalMapBase<_EvalBase, internal::ConstMap >;



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

