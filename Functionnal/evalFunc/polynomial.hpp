/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include <Eigen/Core>
#include "Functionnal/common.h"

namespace functionnal{

template <typename _Scalar, int _Dim >
struct ConstantEvalFunc {
    using Scalar = _Scalar;
    constexpr static const int Dim      = _Dim;
    constexpr static const int Degree   = 0;
    constexpr static const int NbCoeff  = 1;

    typedef Eigen::Matrix<Scalar, Dim, 1> InputVectorType;
    typedef Scalar OutputVectorType;
    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    inline ConstantEvalFunc(){}

    static inline
    OutputVectorType staticEval(const InputVectorType& /*x*/,
                                const Eigen::Ref<const CoeffType>& coeffs)
    {
        return coeffs(0);
    }
};

/*!
 * \brief Base class used to evaluate a n-dimensionnal linear function defined as
 * $ P x^\mathrm{T} + R = 0$
 * with $X = {x_1, ..., x_{n+1}}$
 */
template <typename _Scalar, int _Dim >
struct LinearEvalFunc {

    using Scalar = _Scalar;
    constexpr static const int Degree   = 2;
    constexpr static const int Dim      = _Dim;

    constexpr static const int PSize    = (Dim);
    constexpr static const int RSize    = 1;

    constexpr static const int NbCoeff  = PSize + RSize;

    typedef Eigen::Matrix<Scalar, Dim, 1> InputVectorType;
    typedef Scalar OutputVectorType;
    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    typedef Eigen::Matrix<Scalar, 1,     (Dim)> PType;
    typedef Eigen::Matrix<Scalar, 1,       1  > RType;

    typedef void Derivative;


    inline
    LinearEvalFunc(){}

    static inline Eigen::Map<PType> getPMap(Eigen::Ref<CoeffType> coeffs) {
        return Eigen::Map<PType> (coeffs.data());
    }

    static inline Eigen::Ref<RType> getRMap(Eigen::Ref<CoeffType> coeffs) {
        return Eigen::Ref<RType> (coeffs.template tail<1>());
    }

    static inline Eigen::Map<const PType> getConstPMap(const Eigen::Ref<const CoeffType> coeffs) {
        return Eigen::Map<const PType> (coeffs.data());
    }

    static inline Eigen::Ref<const RType> getConstRMap(const Eigen::Ref<const CoeffType> coeffs) {
        return Eigen::Ref<const RType> (coeffs.template tail<1>());
    }

    /*!
     * \brief Init the quadric as \f$f(x) = x\f$
     * \param coeffs
     */
    static void initCoeffs (Eigen::Ref<CoeffType> coeffs,
                            Eigen::Ref<const PType> p = PType::Ones(),
                            const RType& r = RType(0.))
    {
        getPMap(coeffs) = p;
        getRMap(coeffs) = r;
    }


    static inline
    OutputVectorType staticEval(const InputVectorType& x,
                                const Eigen::Ref<const CoeffType>& coeffs)
    {
        return getConstPMap(coeffs).dot(x) + getConstRMap(coeffs)(0);
    }

};

/*!
 * \brief Base class used to evaluate a n-dimensionnal quadric, defined as
 * $x Q x^\mathrm{T} + P x^\mathrm{T} + R = 0$
 * with $X = {x_1, ..., x_{n+1}}$
 * \see https://en.wikipedia.org/wiki/Quadric
 */
template <typename _Scalar, int _Dim >
struct QuadricEvalFunc {

    using Scalar = _Scalar;
    constexpr static const int Degree   = 2;
    constexpr static const int Dim      = _Dim;

    constexpr static const int QSize    = (Dim)*(Dim);
    constexpr static const int PSize    = (Dim);
    constexpr static const int RSize    = 1;

    constexpr static const int NbCoeff  = QSize + PSize + RSize;

    typedef Eigen::Matrix<Scalar, Dim, 1> InputVectorType;
    typedef Scalar OutputVectorType;
    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    typedef Eigen::Matrix<Scalar, (Dim), (Dim)> QType;
    typedef Eigen::Matrix<Scalar, 1,     (Dim)> PType;
    typedef Eigen::Matrix<Scalar, 1,       1  > RType;

    typedef LinearEvalFunc< Scalar, Dim> Derivative;


    inline
    QuadricEvalFunc(){}

    static inline Eigen::Map<QType> getQMap(Eigen::Ref<CoeffType> coeffs) {
        return Eigen::Map<QType> (coeffs.data());
    }

    static inline Eigen::Ref<PType> getPMap(Eigen::Ref<CoeffType> coeffs) {
        return Eigen::Ref<PType> (coeffs.template block<1, PType::ColsAtCompileTime>(0, QSize));
    }

    static inline Eigen::Ref<RType> getRMap(Eigen::Ref<CoeffType> coeffs) {
        return Eigen::Ref<RType> (coeffs.template tail<1>());
    }

    static inline Eigen::Map<const QType> getConstQMap(const Eigen::Ref<const CoeffType> coeffs) {
        return Eigen::Map<const QType> (coeffs.data());
    }

    static inline Eigen::Ref<const PType> getConstPMap(const Eigen::Ref<const CoeffType> coeffs) {
        return Eigen::Ref<const PType> (coeffs.template block<1, PType::ColsAtCompileTime>(0, QSize));
    }

    static inline Eigen::Ref<const RType> getConstRMap(const Eigen::Ref<const CoeffType> coeffs) {
        return Eigen::Ref<const RType> (coeffs.template tail<1>());
    }

    /*!
     * \brief Init the quadric as \f$f(x) = x\f$
     * \param coeffs
     */
    static void initCoeffs (Eigen::Ref<CoeffType> coeffs,
                            Eigen::Ref<const QType> q = QType::Zero(),
                            Eigen::Ref<const PType> p = PType::Ones(),
                            const RType& r = RType(0.))
    {
        getQMap(coeffs) = q;
        getPMap(coeffs) = p;
        getRMap(coeffs) = r;
    }


    static inline
    OutputVectorType staticEval(const InputVectorType& x,
                                const Eigen::Ref<const CoeffType>& coeffs)
    {
        auto q = getConstQMap(coeffs);
        auto p = getConstPMap(coeffs);
        auto r = getConstRMap(coeffs);

        return (x.transpose()*q).dot(x) + p.dot(x) + r(0);
    }

};

/*!
 * \brief Base class used to evaluate a n-dimensionnal quadratic form, defined
 * as $x Q x^\mathrm{T} = 0$ with $X = {x_1, ..., x_{n+1}}$
 * \see https://en.wikipedia.org/wiki/Quadratic_form
 */
template <typename _Scalar, int _Dim >
struct QuadraticFormEvalFunc {

    using Scalar = _Scalar;
    constexpr static const int Degree   = 2;
    constexpr static const int Dim      = _Dim;

    constexpr static const int QSize    = (Dim)*(Dim);

    constexpr static const int NbCoeff  = QSize;

    typedef Eigen::Matrix<Scalar, Dim, 1> InputVectorType;
    typedef Scalar OutputVectorType;
    typedef Eigen::Matrix<Scalar, 1, NbCoeff> CoeffType;

    typedef Eigen::Matrix<Scalar, (Dim), (Dim)> QType;

    typedef LinearEvalFunc< Scalar, Dim> Derivative;


    inline QuadraticFormEvalFunc(){}

    static inline Eigen::Map<QType> getQMap(Eigen::Ref<CoeffType> coeffs) {
        return Eigen::Map<QType> (coeffs.data());
    }

    static inline Eigen::Map<const QType> getConstQMap(const Eigen::Ref<const CoeffType> coeffs) {
        return Eigen::Map<const QType> (coeffs.data());
    }

    /*!
     * \brief Init the quadric as \f$f(x) = x\f$
     * \param coeffs
     */
    static void initCoeffs (Eigen::Ref<CoeffType> coeffs,
                            Eigen::Ref<const QType> q = QType::Identity())
    {
        getQMap(coeffs) = q;
    }


    static inline
    OutputVectorType staticEval(const InputVectorType& x,
                                const Eigen::Ref<const CoeffType>& coeffs)
    {
        return (x.transpose()*getConstQMap(coeffs)).dot(x);
    }

};

}
