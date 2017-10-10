/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include <Eigen/Core>
#include "Functionnal/common.h"
#include "Functionnal/evalFunc/polynomial.hpp"
#include <iostream>

namespace functionnal{


/*!
 * \brief Base class used to evaluate a bezier curve of arbitrary degree
 * This class defines only the evaluation part as static and member functions.
 *
 * The Dimension parameter defines the control points, while the parametric space
 * is always defined in \f$ [0 : 1]\f$
 *
 */
template <typename _Scalar, int _Degree, int _Dim >
struct BezierEvalFunc{

    typedef _Scalar Scalar;
    enum{
        Degree  = _Degree,
        Dim     = _Dim,
        NbCoeff = _Dim*(_Degree+1) //need the different control points
    };
    //! Vector type in the parametric domain, defines a position in the parametric domain
    typedef Eigen::Matrix<Scalar, 1, 1> InputVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> OutputVectorType;

    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    //! Use ternaire operator to avoid infite template recursion
    typedef BezierEvalFunc<Scalar, ( Degree>0 ? Degree-1 : 0), Dim> Derivative;

    inline
    BezierEvalFunc(){}

    /*!
     * \brief Update the control point as linear interpolations of the end points
     * \param coeffs
     */
    template <typename CType>
    static inline
    void
    initCoeffs(CType& coeffs){
        auto start = coeffs.template head<Dim>();
        // do not use .tail to support larger coeffs vectors
        auto end   = coeffs.template block<1, Dim>(0, Degree*Dim);

        for( int i = 0; i < Degree; ++i ){
            coeffs.template block<1, Dim>(0, i*Dim) =
                    Scalar(i) / Scalar(Degree) * (end - start) + start;
        }
    }

    /*!
     * \brief Evaluate the Bezier curve for the variable \f$x \in [0 : 1]\f$
     */
    static inline
    OutputVectorType staticEval(
            const InputVectorType& x,
            const CoeffType& coeffs)
    {
        return recursiveCall(_Degree, 0, x(0), coeffs);
    }

    static inline
    typename Derivative::CoeffType staticDerivative( const CoeffType& coeffs,
                                                     int DId){
        assert(DId == 0);
        // coefficients are P_{i+1} - P{i}
        return    coeffs.template block<1, Derivative::NbCoeff>(0, 0).array()
                - coeffs.template block<1, Derivative::NbCoeff>(0, Dim).array();
    }

    template <class StreamT>
    static inline
    void staticPrint(StreamT& stream,
                     const CoeffType& coeffs,
                     int /*n*/ = NbCoeff-1)
    {

        Eigen::IOFormat myFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
        for (int d = 0; d != Degree+1; ++d){
            stream << coeffs.template segment<_Dim>(d*Dim).format(myFormat) << std::endl;
        }
    }

private:
    static constexpr
    OutputVectorType recursiveCall(int r, int i,
                                      const Scalar& t,
                                      const CoeffType& coeffsVec){
        return r == 0 ?
                    OutputVectorType(coeffsVec.template segment<Dim>(Dim*i))
                  :
                    (Scalar(1)-t) * recursiveCall(r-1, i,   t, coeffsVec).array() +
                               t  * recursiveCall(r-1, i+1, t, coeffsVec).array();
    }
};

/*!
 * \brief Base class used to evaluate a bezier curve of arbitrary degree
 * using binomial coefficients.
 *
 * The Dimension parameter defines the control points, while the parametric space
 * is always defined in \f$ [0 : 1]\f$
 *
 */
template <typename _Scalar, int _Degree, int _Dim >
struct BezierBinomialEvalFunc{

    typedef _Scalar Scalar;
    enum{
        Degree  = _Degree,
        Dim     = _Dim,
        NbCoeff = _Dim*(_Degree+1) //need the different control points
    };
    //! Vector type in the parametric domain, defines a position in the parametric domain
    typedef Eigen::Matrix<Scalar, 1, 1> InputVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> OutputVectorType;

    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    typedef BezierBinomialEvalFunc<Scalar, ( Degree>0 ? Degree-1 : 0), Dim> Derivative;

    inline
    BezierBinomialEvalFunc(){}

    /*!
     * \brief Update the control point as linear interpolations of the end points
     * \param coeffs
     */
    template <typename CType>
    static inline
    void
    initCoeffs(CType& coeffs){
        BezierEvalFunc<Scalar, Degree, Dim>::initCoeffs(coeffs);
    }

    /*!
     * \brief Evaluate the Bezier curve for the variable \f$x \in [0 : 1]\f$
     */
    static inline
    OutputVectorType staticEval(const InputVectorType& x,
                                   const Eigen::Ref<const CoeffType>& coeffs)
    {
        return recursiveCall(x(0), 0, coeffs);
    }
    /*!
     * \brief Build and return the partial derivative of the Functionnal in
     * dimension Did
     */
    static inline
    typename Derivative::CoeffType staticDerivative( const CoeffType& coeffs,
                                                     int DId = 0){
        // delegate computation to BezierEvalFunc
        return BezierEvalFunc<Scalar, Degree, Dim>::staticDerivative(coeffs, DId);
    }

    template <class StreamT>
    static inline
    void staticPrint(StreamT& stream,
                     const CoeffType& coeffs,
                     int /*n*/ = NbCoeff-1)
    {

        Eigen::IOFormat myFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
        for (int d = 0; d != Degree+1; ++d){
            stream << coeffs.template segment<_Dim>(d*Dim).format(myFormat) << std::endl;
        }
    }

private:
    static constexpr
    Scalar Bi(Scalar t, int i){
        return Utils::PascalTriangle::Row<Degree>::get()[i] * Utils::POW(Scalar(1)-t, Degree-i) * Utils::POW(t, i);
    }

    static constexpr
    OutputVectorType recursiveCall(const Scalar& t,
                                      const int i,
                                      const CoeffType& coeffsVec){
        return i == Degree ?
                      OutputVectorType(
                        coeffsVec.template tail<Dim>()
                        * Bi(t, Degree))
                  :
                      OutputVectorType(
                        coeffsVec.template segment<Dim>(Dim*i)
                        * Bi(t, i))
                    + recursiveCall(t, i+1, coeffsVec);
    }
};

/*!
 * \brief Base class used to evaluate a bezier curve of arbitrary degree
 * using binomial coefficients.
 *
 * The Dimension parameter defines the control points, while the parametric space
 * is always defined in \f$ [0 : 1]\f$
 *
 */
template <typename _Scalar, int _Degree, int _Dim >
struct RationnalBezierEvalFunc{

    typedef _Scalar Scalar;
    enum{
        Degree  = _Degree,
        Dim     = _Dim,
        NbCoeff = (_Dim+1)*(_Degree+1) // need the different weighted control points
                                       // the first _Dim*(_Degree+1) coeffs describe
                                       // the points, and the last (_Degree+1)
                                       // coeffs describe the weights.
    };
    //! Vector type in the parametric domain, defines a position in the parametric domain
    typedef Eigen::Matrix<Scalar, 1, 1> InputVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> OutputVectorType;

    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    static const std::array<int, Degree+1> BinomialCoeffs;

    inline
    RationnalBezierEvalFunc(){}

    /*!
     * \brief Update the control point as linear interpolations of the end points
     * \param coeffs
     */
    template <typename CType>
    static inline
    void
    initCoeffs(CType& coeffs){
        coeffs.template tail<Degree+1>().fill(Scalar(1));
        BezierEvalFunc<Scalar, Degree, Dim>::initCoeffs(coeffs);
    }

    /*!
     * \brief Evaluate the Bezier curve for the variable \f$x \in [0 : 1]\f$
     */
    static inline
    OutputVectorType staticEval(const InputVectorType& x,
                                   const Eigen::Ref<const CoeffType>& coeffs)
    {
        return recursiveCall(x(0), 0,
                             coeffs,
                             // compute the weighted sum of the binomial coefficients
                             sumWBi(x(0), coeffs.template tail<_Degree+1>()));
    }

    template <class StreamT>
    static inline
    void staticPrint(StreamT& stream,
                     const CoeffType& coeffs,
                     int /*n*/ = NbCoeff-1)
    {

        Eigen::IOFormat myFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
        for (int d = 0; d != Degree+1; ++d){
            stream << coeffs.template segment<_Dim>(d*Dim).format(myFormat) << std::endl;
        }
    }

private:
    //! Binomial functor
    static constexpr
    Scalar Bi(Scalar t, int i){
        return Utils::PascalTriangle::Row<_Degree>::get()[i]
                * Utils::POW(Scalar(1)-t, Degree-i)
                * Utils::POW(t, i);
    }

    template<typename Derived>
    static
    Scalar sumWBi(Scalar t, const Eigen::MatrixBase<Derived>& /*w*/){
        Eigen::Map<Eigen::Matrix<int, 1, Degree+1> > BinomialMap
                (Utils::PascalTriangle::Row<Degree>::get().data());

        Scalar sumW = 0;
        for(int i = 0; i != Degree+1; ++i){
            sumW += BinomialMap(i)
                    * Utils::POW(Scalar(1)-t, Degree-i)
                    * Utils::POW(t, i);
        }

        return sumW;
    }

    static constexpr
    OutputVectorType recursiveCall(const Scalar& t,
                                      const int i,
                                      const CoeffType& coeffsVec,
                                      const Scalar& wSum){
        return i == Degree ?
                      OutputVectorType(
                        coeffsVec.template segment<Dim>(Dim*i)
                        * Bi(t, Degree)   // binomial coefficient
                        * coeffsVec.template tail<Degree+1>()(i) // element weight
                        / wSum)           // normalization
                  :
                      OutputVectorType(
                        coeffsVec.template segment<Dim>(Dim*i)
                        * Bi(t, i)        // binomial coefficient
                        * coeffsVec.template tail<Degree+1>()(i) // element weight
                        / wSum)           // normalization

                    + recursiveCall(t, i+1, coeffsVec, wSum);
    }
};

}
