/*
    This Source Code Form is subject to the terms of the 
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed 
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
    
    \author Nicolas Mellado nmellado0@gmail.com
*/

#ifndef FUNCTIONNAL_EVAL_FUNC_H
#define FUNCTIONNAL_EVAL_FUNC_H

#include <Eigen/Core>
#include "common.h"
#include <iostream>

namespace Functionnal{
/*!
 * \brief Base class used to evaluate a polynomial function of arbitrary degree
 * in n dimension. Here embedding and embedded spaces share the same dimensions.
 *
 */
template <typename _Scalar, int _Degree, int _Dim >
struct PolynomialEvalFunc{

    typedef _Scalar Scalar;
    enum{
        Degree  = _Degree,
        Dim     = _Dim,
        NbCoeff = (_Degree+1) * Dim
    };
    //! Vector type in the parametric domain, defines a position in the parametric domain
    typedef Eigen::Matrix<Scalar, Dim, 1> EmbeddedVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> EmbeddingVectorType;

    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    typedef PolynomialEvalFunc<Scalar, Degree-1, Dim> Derivative;

    inline
    PolynomialEvalFunc(){}

    /*!
     * \brief Init the polynomial as \f$f(x) = x\f$
     * \param coeffs
     */
    template <typename CType>
    static inline
    void
    initCoeffs(CType& coeffs){
        coeffs = CoeffType::Zero();
        if (NbCoeff >= 2)
            coeffs.template block<1, _Dim>(0, _Dim) = EmbeddingVectorType::Ones();
    }

    /*!
     * \brief Evaluate the Polynomial for the variable x
     *
     * \todo Use Eigen::CwiseNullaryOp to initialize M and avoid the for loop
     *
     * \param coeffs Storing \f$coeffs = [a0, a1, ..., a_{Degree}]\f$     *
     *          Initialized as
     * \code
     * curve.coeffs << 0, 1, 2, 3; // x + 2x^2 + 3x^3
     * \endcode
     *
     * \note We need to use a map here in order to play with
     * coefficient modified from outside (e.g. by the Eigen solver)
     * without manual synchronisation
     */
    static inline
    EmbeddingVectorType staticEval(const EmbeddedVectorType& x,
                                   const Eigen::Ref<const CoeffType>& coeffs)
    {
        // generate evaluation matrix with:
        // X =
        // 1    1    1
        // x    y    z
        // x^2  y^2  z^2
        // ...  ...  ...
        // x^n  y^n  z^n
        // where n is the polynom degree
        Eigen::Matrix<Scalar, NbCoeff, Dim> M;
        M.row(0) = Eigen::Matrix<Scalar, 1, Dim>::Ones().array();
        for(int i = 1; i != NbCoeff; ++i)
            M.row(i) = M.row(i-1).array() * x.transpose().array();

        // eval
        return coeffs * M; //coeffs * M;
    }

    template <class StreamT>
    static inline
    void staticPrint(StreamT& stream,
                     const CoeffType& coeffs,
                     int n = NbCoeff-1)
    {
        if(n == 0){
            stream << "[" << coeffs.col(n).transpose() << "]" << std::flush;
        }else if (n == 1){
            if ( coeffs(n) != 0)
                stream << "["  << coeffs.col(n).transpose() << "] x + ";
            staticPrint(stream, coeffs,n-1);
        }else {
            if ( coeffs(n) != 0)
                stream << "["  << coeffs.col(n).transpose() << "] x^" << n << " + ";
            staticPrint(stream, coeffs, n-1);
        }
    }
};

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
    typedef Eigen::Matrix<Scalar, 1, 1> EmbeddedVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> EmbeddingVectorType;

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
    EmbeddingVectorType staticEval(
            const EmbeddedVectorType& x,
            const CoeffType& coeffs)
    {
        return recursiveCall(_Degree, 0, x(0), coeffs);
    }

    static inline
    typename Derivative::CoeffType staticDerivative( const CoeffType& coeffs){
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
    EmbeddingVectorType recursiveCall(int r, int i,
                                      const Scalar& t,
                                      const CoeffType& coeffsVec){
        return r == 0 ?
                    EmbeddingVectorType(coeffsVec.template segment<Dim>(Dim*i))
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
    typedef Eigen::Matrix<Scalar, 1, 1> EmbeddedVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> EmbeddingVectorType;

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
    EmbeddingVectorType staticEval(const EmbeddedVectorType& x,
                                   const Eigen::Ref<const CoeffType>& coeffs)
    {
        return recursiveCall(x(0), 0, coeffs);
    }

    static inline
    typename Derivative::CoeffType staticDerivative( const CoeffType& coeffs){
        // delegate computation to BezierEvalFunc
        return BezierEvalFunc<Scalar, Degree, Dim>::staticDerivative(coeffs);
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
    EmbeddingVectorType recursiveCall(const Scalar& t,
                                      const int i,
                                      const CoeffType& coeffsVec){
        return i == Degree ?
                      EmbeddingVectorType(
                        coeffsVec.template tail<Dim>()
                        * Bi(t, Degree))
                  :
                      EmbeddingVectorType(
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
    typedef Eigen::Matrix<Scalar, 1, 1> EmbeddedVectorType;

    //! Vector type in the embedding domain
    typedef Eigen::Matrix<Scalar, Dim, 1> EmbeddingVectorType;

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
    EmbeddingVectorType staticEval(const EmbeddedVectorType& x,
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
    EmbeddingVectorType recursiveCall(const Scalar& t,
                                      const int i,
                                      const CoeffType& coeffsVec,
                                      const Scalar& wSum){
        return i == Degree ?
                      EmbeddingVectorType(
                        coeffsVec.template segment<Dim>(Dim*i)
                        * Bi(t, Degree)   // binomial coefficient
                        * coeffsVec.template tail<Degree+1>()(i) // element weight
                        / wSum)           // normalization
                  :
                      EmbeddingVectorType(
                        coeffsVec.template segment<Dim>(Dim*i)
                        * Bi(t, i)        // binomial coefficient
                        * coeffsVec.template tail<Degree+1>()(i) // element weight
                        / wSum)           // normalization

                    + recursiveCall(t, i+1, coeffsVec, wSum);
    }
};

}
#endif // FUNCTIONNAL_EVAL_FUNC_H
