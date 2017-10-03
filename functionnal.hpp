/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include "functionnal.h"

namespace functionnal{
namespace internal{
} // namespace internal
} // namespace Functionnal

template <class StreamT,
          typename _EvalBase,
          template <typename> class _MapType>
inline
StreamT& operator<< (StreamT& stream,
                     const functionnal::internal::
                     FunctionnalMapBase<_EvalBase, _MapType> &p){
    p.print(stream);
    return stream;
}

template <class StreamT,
          typename _EvalBase>
inline
StreamT& operator<< (StreamT& stream,
                     const functionnal::
                     Functionnal<_EvalBase> &p){
    p.print(stream);
    return stream;
}

