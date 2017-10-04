/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include <Eigen/Core>
#include "Functionnal/common.h"
#include <iostream>

namespace Functionnal{
/*!
 * \brief Base class used to evaluate a polynomial function of arbitrary degree
 * in n dimension. Here embedding and embedded spaces share the same dimensions.
 * For multi-dimensionnal vectors, each dimension is treated a scalar value and
 * evaluated with the same polynomial. For instance, when dimension = 3,
 * polynomials with degree 3 are evaluated as:
 * \code
 * c =  a0   a1   a2   ....  an
 * X =
 *      1    1    1
 *      x    y    z
 *      x^2  y^2  z^2
 *      ...  ...  ...
 *      x^n  y^n  z^n
 *
 * f(x) = c * X;
 * \endcode
 */
template <typename _Scalar, int _Degree, int _Dim >
struct CWisePolynomialEvalFunc{

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

    typedef CWisePolynomialEvalFunc<Scalar, Degree-1, Dim> Derivative;

    inline
    CWisePolynomialEvalFunc(){}

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


}
