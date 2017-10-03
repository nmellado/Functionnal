/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include "Functionnal/functionnal.h"

#include <iostream>


namespace Functionnal{

namespace internal{

/*!
 * \brief Bezier curve with constrained endpoints.
 *
 * The purpose of this class is to hide the end-points in the coefficient map.
 *
 * \warning End-points are *not* stored in the mapped memory, but in the
 * instance. So don't destroy it...
 *
 * ConstrainedBezier can be copied at construction time in order to copy the
 * end-points. However, you still need to give the memory area you want to be
 * mapped to in order to allow a correct evaluation (Eigen::Map cannot be
 * duplicated).
 *
 */
template <typename _Scalar, int _Degree, int _Dim,
          template <typename, int, int> class  _EvalFuncT,
          template <typename> class MapType >
struct ConstrainedBezierMapBase {
private:
    //! This type defines the coefficient types used for evaluation
    typedef _EvalFuncT< _Scalar, _Degree, _Dim> EvalBase;

    //! This type defines the coefficient types that will be optimized
    typedef _EvalFuncT< _Scalar, _Degree-2, _Dim> ConstrainedEvalBase;

public:
    typedef typename EvalBase::Scalar Scalar;
    enum{
        Dim     = EvalBase::Dim,
        NbCoeff = ConstrainedEvalBase::NbCoeff,
        Degree  = EvalBase::Degree
    };

    typedef typename EvalBase::EmbeddedVectorType EmbeddedVectorType;
    typedef typename EvalBase::EmbeddingVectorType EmbeddingVectorType;
    typedef typename ConstrainedEvalBase::CoeffType CoeffType;
    typedef Functionnal<_EvalFuncT < Scalar,
                                     EvalBase::Derivative::Degree,
                                     EvalBase::Derivative::Dim> > Derivative;

    //! Map to the unconstrained coefficients
    MapType<CoeffType> coeffs;

    EmbeddingVectorType startPoint, //! <\brief Constrained starting point
                        endPoint;   //! <\brief Constrained end point


    inline
    ConstrainedBezierMapBase(Scalar* data,
                             const EmbeddingVectorType& start = EmbeddingVectorType::Zero(),
                             const EmbeddingVectorType& end   = EmbeddingVectorType::Zero())
        : coeffs(data),
          startPoint(start),
          endPoint(end) {}

    inline
    ConstrainedBezierMapBase(const Scalar* data,
                             const EmbeddingVectorType& start = EmbeddingVectorType::Zero(),
                             const EmbeddingVectorType& end   = EmbeddingVectorType::Zero())
        : coeffs(data),
          startPoint(start),
          endPoint(end) {}
    inline
    ConstrainedBezierMapBase(const Scalar* data,
                             const ConstrainedBezierMapBase& other)
        : coeffs(data),
          startPoint(other.startPoint),
          endPoint(other.endPoint) {}

    inline
    void initCoeffs()
    {
        typename EvalBase::CoeffType ccoeffs = concatCoeffs();
        EvalBase::initCoeffs(ccoeffs);
        ConcatCoeff<EvalBase>::unroll(ccoeffs, coeffs);
//        // copy vertex coordinates
//        constexpr int Deg = (Degree-1);
//        coeffs.template head<Dim*Deg>() = ccoeffs.template block<1, Dim*Deg>(0, Dim);
//        // copy vertex weights
//        coeffs.template block<1,Deg>(0,Dim*Deg) = ccoeffs.template block<1, Deg>(0, Dim*(Deg+2));
    }

    inline
    EmbeddingVectorType eval(const EmbeddedVectorType& x) const{
        return EvalBase::staticEval(x, concatCoeffs());
    }

    /*!
     * \brief Build and return de derivative of the Functionnal
     */
    inline
    Derivative derivative() const {
        return Derivative(EvalBase::staticDerivative(concatCoeffs()));
    }

private:
    /*!
     * \brief default version do nothing and crashes at compile time
     * \see specialized versions below
     */
    inline
    typename EvalBase::CoeffType concatCoeffs() const{
        return ConcatCoeff<EvalBase>::roll(startPoint, endPoint, coeffs);
    }

    template <typename Base, typename DUMMY = void>
    struct ConcatCoeff {
        static inline void unroll(
                typename EvalBase::CoeffType& input,
                MapType<CoeffType>&output){
            // copy vertex coordinates
            constexpr int Deg = (Degree-1);
            output.template head<Dim*Deg>() = input.template block<1, Dim*Deg>(0, Dim);
        }

        static inline typename EvalBase::CoeffType roll(
                const EmbeddingVectorType& /*start*/,
                const EmbeddingVectorType& /*end*/,
                const CoeffType& /*c*/) {
            //static_assert(false, "Unsupported Evaluation class");
            return EvalBase::CoeffType();
        }
    };

    template <typename DUMMY>
    struct ConcatCoeff<RationnalBezierEvalFunc< _Scalar, _Degree, _Dim>, DUMMY>{
        static inline void unroll(
                typename EvalBase::CoeffType& input,
                MapType<CoeffType>&output){
            // copy vertex coordinates
            constexpr int Deg = (Degree-1);
            output.template head<Dim*Deg>() = input.template block<1, Dim*Deg>(0, Dim);
            // copy vertex weights
            output.template block<1,Deg>(0,Dim*Deg) = input.template block<1, Deg>(0, Dim*(Deg+2));
        }

        static typename EvalBase::CoeffType roll(
                const EmbeddingVectorType& start,
                const EmbeddingVectorType& end,
                const CoeffType& c) {
            typename EvalBase::CoeffType cCoeffs;

            const int nBControlPoint = Degree-1;

            int offset = 0;
            // set starting point, middle coeffs, end point
            cCoeffs.template head<Dim>()
                    = start;
            cCoeffs.template segment<Dim*nBControlPoint>
                    (offset += Dim)
                    = c.template head<Dim*nBControlPoint>();
            cCoeffs.template segment<Dim>
                    (offset += Dim*nBControlPoint)
                    = end;

            // set weights
            cCoeffs.template segment<1>(offset += Dim)
                    << 1;
            cCoeffs.template segment<nBControlPoint>(offset += 1)
                  = c.template tail<nBControlPoint>();
            cCoeffs.template segment<1>(offset += nBControlPoint) << 1;

            return cCoeffs;
        }
    };

    template <typename DUMMY>
    struct ConcatCoeff<BezierEvalFunc< _Scalar, _Degree, _Dim>, DUMMY>{
        static inline void unroll(
                typename EvalBase::CoeffType& input,
                MapType<CoeffType>&output){
            // copy vertex coordinates
            constexpr int Deg = (Degree-1);
            output.template head<Dim*Deg>() = input.template block<1, Dim*Deg>(0, Dim);
        }

        static typename EvalBase::CoeffType roll(
                const EmbeddingVectorType& start,
                const EmbeddingVectorType& end,
                const CoeffType& c) {
            typename EvalBase::CoeffType cCoeffs;
            cCoeffs.template head<EmbeddingVectorType::RowsAtCompileTime>()
                    = start;
            cCoeffs.template tail<EmbeddingVectorType::RowsAtCompileTime>()
                    = end;
            cCoeffs.template segment<NbCoeff>(EmbeddingVectorType::RowsAtCompileTime)
                    = c;

            return cCoeffs;
        }
    };

    template <typename DUMMY>
    struct ConcatCoeff<BezierBinomialEvalFunc< _Scalar, _Degree, _Dim>, DUMMY>{
        static inline void unroll(
                typename EvalBase::CoeffType& input,
                MapType<CoeffType>&output){
            // copy vertex coordinates
            constexpr int Deg = (Degree-1);
            output.template head<Dim*Deg>() = input.template block<1, Dim*Deg>(0, Dim);
        }

        static typename EvalBase::CoeffType roll(
                const EmbeddingVectorType& start,
                const EmbeddingVectorType& end,
                const CoeffType& c) {
            typename EvalBase::CoeffType cCoeffs;
            cCoeffs.template head<EmbeddingVectorType::RowsAtCompileTime>()
                    = start;
            cCoeffs.template tail<EmbeddingVectorType::RowsAtCompileTime>()
                    = end;
            cCoeffs.template segment<NbCoeff>(EmbeddingVectorType::RowsAtCompileTime)
                    = c;

            return cCoeffs;
        }
    };

public:
    template <class StreamT>
    inline
    void print(StreamT& stream) const{
        return EvalBase::staticPrint(stream, concatCoeffs());
    }
};

} // namespace internal


//  ConstrainedBezier typedefs
template<typename _Scalar, int _Degree, int _Dim>
using ConstrainedBezierMap =
internal::ConstrainedBezierMapBase< _Scalar, _Degree, _Dim, BezierBinomialEvalFunc, internal::Map >;

template<typename _Scalar, int _Degree, int _Dim>
using ConstConstrainedBezierMap =
internal::ConstrainedBezierMapBase< _Scalar, _Degree, _Dim, BezierBinomialEvalFunc, internal::ConstMap >;

//  ConstrainedBezier typedefs
template<typename _Scalar, int _Degree, int _Dim>
using ConstrainedRationnalBezierMap =
internal::ConstrainedBezierMapBase< _Scalar, _Degree, _Dim, RationnalBezierEvalFunc, internal::Map >;

template<typename _Scalar, int _Degree, int _Dim>
using ConstConstrainedRationnalBezierMap =
internal::ConstrainedBezierMapBase< _Scalar, _Degree, _Dim, RationnalBezierEvalFunc, internal::ConstMap >;

} // namespace Functionnal


template <class StreamT,
          typename _S, int _De, int _Di,
          template <typename, int, int> class _EvalFuncT,
          template <typename> class _MapT>
inline
StreamT& operator<< (StreamT& stream,
                     const Functionnal::internal::
                     ConstrainedBezierMapBase<_S,_De,_Di,_EvalFuncT,_MapT> &p){
    p.print(stream);
    return stream;
}


