/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include "Functionnal/functionnal.h"

namespace functionnal{
namespace internal{
} // namespace internal
} // namespace Functionnal

template <class StreamT,
          typename EvalBase, typename CoeffStorageType >
inline
StreamT& operator<< (StreamT& stream,
                     functionnal::FunctionnalBase<EvalBase, CoeffStorageType> &p){
    p.print(stream);
    return stream;
}

