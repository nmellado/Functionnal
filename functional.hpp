/*
    This Source Code Form is subject to the terms of the 
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed 
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
    
    \author Nicolas Mellado nmellado0@gmail.com
*/

#ifndef FUNCTIONNAL_HPP
#define FUNCTIONNAL_HPP

// The following directive is used only to get code completion and doesn't pass
// the pre-compilation stage.
#ifndef FUNCTIONNAL_H
#include "functionnal.h"
#endif

namespace Functionnal{
namespace internal{
} // namespace internal
} // namespace Functionnal

template <class StreamT,
          typename _EvalBase,
          template <typename> class _MapType>
inline
StreamT& operator<< (StreamT& stream,
                     const Functionnal::internal::
                     FunctionnalMapBase<_EvalBase, _MapType> &p){
    p.print(stream);
    return stream;
}

template <class StreamT,
          typename _EvalBase>
inline
StreamT& operator<< (StreamT& stream,
                     const Functionnal::
                     Functionnal<_EvalBase> &p){
    p.print(stream);
    return stream;
}


#endif // FUNCTIONNAL_HPP
