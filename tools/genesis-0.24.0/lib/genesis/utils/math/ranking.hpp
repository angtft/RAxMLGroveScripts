#ifndef GENESIS_UTILS_MATH_RANKING_H_
#define GENESIS_UTILS_MATH_RANKING_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2020 Lucas Czech and HITS gGmbH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief
 *
 * @file
 * @ingroup utils
 */

#include "genesis/utils/core/algorithm.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

namespace genesis {
namespace utils {

// =================================================================================================
//     Ranking Standard
// =================================================================================================

/**
 * @brief Return the ranking of the values in the given range, using Standard competition ranking
 * ("1224" ranking).
 *
 * See https://en.wikipedia.org/wiki/Ranking for details.
 *
 * @see ranking_modified(), ranking_dense(), ranking_ordinal(), ranking_fractional() for other
 * ranking methods.
 */
template <class RandomAccessIterator>
std::vector<size_t> ranking_standard( RandomAccessIterator first, RandomAccessIterator last )
{
    // Prepare result, and get the sorting order of the vector.
    auto const size = static_cast<size_t>( std::distance( first, last ));
    auto result = std::vector<size_t>( size, 1 );
    auto const order = stable_sort_indices( first, last );

    // Shortcuts for better readability.
    auto ordered_value = [&]( size_t i ){
        return *( first + order[i] );
    };
    auto ordered_result = [&]( size_t i ) -> size_t& {
        return result[ order[i] ];
    };

    // Calculate ranks.
    for( size_t i = 1; i < size; ++i ) {

        // Same values get the same rank. The next bigger one continues at the current i.
        if( ordered_value( i ) == ordered_value( i - 1 ) ) {
            ordered_result( i ) = ordered_result( i - 1 );
        } else {
            ordered_result( i ) = i + 1;
        }
    }

    return result;
}

/**
 * @copydoc ranking_standard( RandomAccessIterator first, RandomAccessIterator last )
 */
inline std::vector<size_t> ranking_standard( std::vector<double> const& vec )
{
    return ranking_standard( vec.begin(), vec.end() );
}

// =================================================================================================
//     Ranking Modified
// =================================================================================================

/**
 * @brief Return the ranking of the values in the given range, using Modified competition ranking
 * ("1334" ranking).
 *
 * See https://en.wikipedia.org/wiki/Ranking for details.
 *
 * @see ranking_standard(), ranking_dense(), ranking_ordinal(), ranking_fractional() for other
 * ranking methods.
 */
template <class RandomAccessIterator>
std::vector<size_t> ranking_modified( RandomAccessIterator first, RandomAccessIterator last )
{
    // Prepare result, and get the sorting order of the vector.
    auto const size = static_cast<size_t>( std::distance( first, last ));
    auto result = std::vector<size_t>( size, 1 );
    auto const order = stable_sort_indices( first, last );

    // Shortcuts for better readability.
    auto ordered_value = [&]( size_t i ){
        return *( first + order[i] );
    };
    auto ordered_result = [&]( size_t i ) -> size_t& {
        return result[ order[i] ];
    };

    // Calculate ranks. The loop variable is incremented at the end.
    for( size_t i = 0; i < size; ) {

        // Look ahead: How often does the value occur?
        size_t j = 1;
        while( i+j < size && ordered_value(i+j) == ordered_value(i) ) {
            ++j;
        }

        // Set the j-next entries.
        for( size_t k = 0; k < j; ++k ) {
            ordered_result( i + k ) = i + j;
        }

        // We can skip the j-next loop iterations, as we just set their values
        i += j;
    }

    return result;
}

/**
 * @copydoc ranking_modified( RandomAccessIterator first, RandomAccessIterator last )
 */
inline std::vector<size_t> ranking_modified( std::vector<double> const& vec )
{
    return ranking_modified( vec.begin(), vec.end() );
}

// =================================================================================================
//     Ranking Dense
// =================================================================================================

/**
 * @brief Return the ranking of the values in the given range, using Dense ranking ("1223" ranking).
 *
 * See https://en.wikipedia.org/wiki/Ranking for details.
 *
 * @see ranking_standard(), ranking_modified(), ranking_ordinal(), ranking_fractional() for other
 * ranking methods.
 */
template <class RandomAccessIterator>
std::vector<size_t> ranking_dense( RandomAccessIterator first, RandomAccessIterator last )
{
    // Prepare result, and get the sorting order of the vector.
    auto const size = static_cast<size_t>( std::distance( first, last ));
    auto result = std::vector<size_t>( size, 1 );
    auto const order = stable_sort_indices( first, last );

    // Shortcuts for better readability.
    auto ordered_value = [&]( size_t i ){
        return *( first + order[i] );
    };
    auto ordered_result = [&]( size_t i ) -> size_t& {
        return result[ order[i] ];
    };

    // Calculate ranks.
    for( size_t i = 1; i < size; ++i ) {

        // Same values get the same rank. The next bigger one continues by incrementing.
        if( ordered_value( i ) == ordered_value( i - 1 ) ) {
            ordered_result( i ) = ordered_result( i - 1 );
        } else {
            ordered_result( i ) = ordered_result( i - 1 ) + 1;
        }
    }

    return result;
}

/**
 * @copydoc ranking_dense( RandomAccessIterator first, RandomAccessIterator last )
 */
inline std::vector<size_t> ranking_dense( std::vector<double> const& vec )
{
    return ranking_dense( vec.begin(), vec.end() );
}

// =================================================================================================
//     Ranking Ordinal
// =================================================================================================

/**
 * @brief Return the ranking of the values in the given range, using Ordinal ranking ("1234" ranking).
 *
 * See https://en.wikipedia.org/wiki/Ranking for details.
 *
 * @see ranking_standard(), ranking_modified(), ranking_dense(), ranking_fractional() for other
 * ranking methods.
 */
template <class RandomAccessIterator>
std::vector<size_t> ranking_ordinal( RandomAccessIterator first, RandomAccessIterator last )
{
    // Prepare result, and get the sorting order of the vector.
    auto const size = static_cast<size_t>( std::distance( first, last ));
    auto result = std::vector<size_t>( size, 1 );
    auto const order = stable_sort_indices( first, last );

    // Shortcuts for better readability.
    auto ordered_result = [&]( size_t i ) -> size_t& {
        return result[ order[i] ];
    };

    // Calculate ranks. This is simply the order plus 1 (as ranks are 1-based).
    for( size_t i = 0; i < size; ++i ) {
        ordered_result( i ) = i + 1;
    }

    return result;
}

/**
 * @copydoc ranking_ordinal( RandomAccessIterator first, RandomAccessIterator last )
 */
inline std::vector<size_t> ranking_ordinal( std::vector<double> const& vec )
{
    return ranking_ordinal( vec.begin(), vec.end() );
}

// =================================================================================================
//     Ranking Fractional
// =================================================================================================

/**
 * @brief Return the ranking of the values in the given range, using Fractional ranking
 * ("1 2.5 2.5 4" ranking).
 *
 * See https://en.wikipedia.org/wiki/Ranking for details. This is the only raking method that
 * returns float values instead of integer values.
 *
 * @see ranking_standard(), ranking_modified(), ranking_dense(), ranking_ordinal() for other
 * ranking methods.
 */
template <class RandomAccessIterator>
std::vector<double> ranking_fractional( RandomAccessIterator first, RandomAccessIterator last )
{
    // Prepare result, and get the sorting order of the vector.
    auto const size = static_cast<size_t>( std::distance( first, last ));
    auto result = std::vector<double>( size, 1 );
    auto const order = stable_sort_indices( first, last );

    // Shortcuts for better readability.
    auto ordered_value = [&]( size_t i ){
        return *( first + order[i] );
    };
    auto ordered_result = [&]( size_t i ) -> double& {
        return result[ order[i] ];
    };

    // Calculate the average of the sum of numbers in the given inclusive range.
    auto sum_avg = []( size_t l, size_t r )
    {
        assert( l > 0 );
        assert( r > 0 );
        assert( l <= r );

        // Example:  l == 7, r == 9
        // We want:  (7 + 8 + 9) / 3 = 8.0
        // Upper:    1+2+3+4+5+6+7+8+9 = 45
        // Lower:    1+2+3+4+5+6       = 21
        // Diff:     45 - 21 = 24
        // Count:    9 - 7 + 1 = 3
        // Result:   24 / 3 = 8
        auto const upper = r * ( r + 1 ) / 2;
        auto const lower = ( l - 1 ) * l / 2;
        return static_cast<double>( upper - lower ) / static_cast<double>( r - l + 1 );
    };

    // Calculate ranks. The loop variable is incremented at the end.
    for( size_t i = 0; i < size; ) {

        // Look ahead: How often does the value occur?
        size_t j = 1;
        while( i+j < size && ordered_value(i+j) == ordered_value(i) ) {
            ++j;
        }

        // Set the j-next entries.
        auto entry = sum_avg( i + 1, i + j );
        for( size_t k = 0; k < j; ++k ) {
            ordered_result( i + k ) = entry;
        }

        // We can skip the j-next loop iterations, as we just set their values
        i += j;
    }

    return result;
}

/**
 * @copydoc ranking_fractional( RandomAccessIterator first, RandomAccessIterator last )
 */
inline std::vector<double> ranking_fractional( std::vector<double> const& vec )
{
    return ranking_fractional( vec.begin(), vec.end() );
}

} // namespace utils
} // namespace genesis

#endif // include guard
