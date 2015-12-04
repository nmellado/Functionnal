/*
    This Source Code Form is subject to the terms of the 
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed 
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
    
    \author Nicolas Mellado nmellado0@gmail.com
*/

#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include <Eigen/Core>

namespace functional{

/*!
 *\brief Helper class defining a Functional object storing its own memory
 */
template < typename _EvalBase >
struct Functional {
    typedef _EvalBase EvalBase;
    typedef typename EvalBase::Scalar Scalar;
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = EvalBase::NbCoeff
    };
    typedef typename EvalBase::EmbeddedVectorType EmbeddedVectorType;
    typedef typename EvalBase::EmbeddingVectorType EmbeddingVectorType;
    typedef typename EvalBase::CoeffType CoeffType;
    typedef Functional<typename EvalBase::Derivative> Derivative;

    //! Actual storage of the coefficients
    CoeffType coeffs;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline
    Functional(){ }


    template <typename Derived>
    inline
    Functional(const Eigen::MatrixBase<Derived>& otherCoeffs)
          : coeffs(otherCoeffs){
    }   // update with input param


    inline
    Functional(const std::array<Scalar, Dim*NbCoeff>& coeffArray)
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
     * \brief Build and return de derivative of the Functional
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
 *\brief Helper class defining a Functional object using external memory
 */
template <typename _EvalBase,
          template <typename> class MapType >
struct FunctionalMapBase {

    typedef _EvalBase EvalBase;
    typedef typename EvalBase::Scalar Scalar;
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = EvalBase::NbCoeff
    };
    typedef typename EvalBase::EmbeddedVectorType EmbeddedVectorType;
    typedef typename EvalBase::EmbeddingVectorType EmbeddingVectorType;
    typedef typename EvalBase::CoeffType CoeffType;
    typedef Functional<typename EvalBase::Derivative> Derivative;

    //! Map to the actual coefficients
    MapType<CoeffType> coeffs;


    inline FunctionalMapBase(Scalar* data)
        :coeffs(data){}
    inline FunctionalMapBase(const Scalar* data)
        :coeffs(data){}

    inline Functional<EvalBase> asFunctional() const{
        return Functional<EvalBase>(coeffs);
    }

    //! \brief Init with default coefficients
    inline
    void initCoeffs()
    { EvalBase::initCoeffs(coeffs); }

    inline
    EmbeddingVectorType eval(const EmbeddedVectorType& x) const{
        return EvalBase::staticEval(x, coeffs);
    }

    //! Build and return de derivative of the FunctionalMap as Functional
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
} // namespace Functional

////////////////////////////////////////////////////////////////////////////////
// Predefined types

#include "functionalEvalFunc.h"


namespace functional{

template<typename _EvalBase>
using FunctionalMap =
internal::FunctionalMapBase<_EvalBase, internal::Map >;

template<typename _EvalBase>
using ConstFunctionalMap =
internal::FunctionalMapBase<_EvalBase, internal::ConstMap >;



// Polynomial typedefs
template<typename _Scalar, int _Degree, int _Dim>
using Polynomial =
Functional< PolynomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using PolynomialMap =
FunctionalMap< PolynomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using ConstPolynomialMap =
ConstFunctionalMap< PolynomialEvalFunc< _Scalar, _Degree, _Dim > >;


// Bezier typedefs
template<typename _Scalar, int _Degree, int _Dim>
using Bezier =
Functional< BezierBinomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using BezierMap =
FunctionalMap< BezierBinomialEvalFunc< _Scalar, _Degree, _Dim > >;

template<typename _Scalar, int _Degree, int _Dim>
using ConstBezierMap =
ConstFunctionalMap< BezierBinomialEvalFunc< _Scalar, _Degree, _Dim > >;

} // namespace Functional

#include "functional.hpp"

#endif // FUNCTIONAL_H
