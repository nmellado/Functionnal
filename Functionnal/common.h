/*
    This Source Code Form is subject to the terms of the
    Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

    \author Nicolas Mellado nmellado0@gmail.com
*/

#pragma once

#include <array>
#include <cmath>


namespace Functionnal {
namespace Utils{

//! \brief Compile time pow
template<typename baseT, typename expoT>
constexpr baseT POW(baseT base, expoT expo)
{
    return (expo != 0 )? base * POW(base, expo -1) : 1;
}


namespace PascalTriangle{

namespace internal{

/*!
 * \todo Add static assert to check that p <= r
 * \param row in the pascal triangle
 * \param column in the pascal triangle
 *
 *
 */
template <typename IntT>
constexpr
IntT value(IntT r, IntT p) {
    return  (p == 0 || p == r) ? 1
                               : value(r-1, p-1) + value(r-1, p);
}

template<int Size, int... i>
constexpr std::array<int, Size+1> fs() {
    return  std::array<int, Size+1> {{ value(Size, i)... }};
}
template<int Size, int... i> struct S;

template<int Size, int... i> struct S<Size, 0,i...>
{ static constexpr std::array<int, Size+1> get() { return fs<Size, 0,i...>(); } };

template<int Size, int i, int... j> struct S<Size, i,j...>
{ static constexpr std::array<int, Size+1> get() { return S<Size, i-1,i,j...>::get(); } };

} // namespace internal

/*!
 * \todo Add way to use other types than int.
 *
 * This could be done using a main class templated with the input type,
 * and nested classes implementing the get function.
 * Hence, the main class would be specialized, and the int type propagated to
 * nested ones without specialization
 *
 * \todo Check aforementionned instructions
 */
template<int Size> struct Row : public internal::S<Size, Size>{};

}
} // namespace utils
} // namespace functionnal
