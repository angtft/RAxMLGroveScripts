/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2019 Lucas Czech and HITS gGmbH

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
 * @ingroup sequence
 */

#include "genesis/sequence/functions/consensus.hpp"

#include "genesis/sequence/counts.hpp"
#include "genesis/sequence/sequence_set.hpp"
#include "genesis/sequence/sequence.hpp"
#include "genesis/sequence/functions/codes.hpp"
#include "genesis/sequence/functions/functions.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace genesis {
namespace sequence {

// =================================================================================================
//     Helper Template
// =================================================================================================

/**
 * @brief Helper struct to store characters sorted by their count.
 */
using CountPair = std::pair< size_t, char >;

/**
 * @brief Local helper function to sort a CountPair.
 *
 * The comparatorfirst sorts by count, and for equal counts, sorts alphanumercially by the char.
 * Thus, gaps are always first (their ascii code is smaller than all letters).
 */
static bool CountPairComparator ( CountPair const& lhs, CountPair const& rhs )
{
    return ( lhs.first > rhs.first ) || ( lhs.first == rhs.first && lhs.second < rhs.second );
}

/**
 * @brief Local helper function template that handles the common code for the nucleic acid
 * consensus sequence functions.
 */
template< typename CharSelector >
std::string consensus_sequence_template(
    SiteCounts const&     counts,
    bool const            allow_gaps,
    CharSelector const&   char_selector
) {
    // Check whether the counts object uses the needed character codes for this function.
    // The characters in the counts object are sorted, so we can directly compare then like this.
    if( counts.characters() != "ACGT" ) {
        throw std::runtime_error(
            "Computation of this consensus sequence function is only meant "
            "for nucleic acid codes (ACGT)."
        );
    }

    // Prepare result.
    std::string result;
    result.reserve( counts.length() );

    // Prepare some constants for simplicity.
    auto const chars     = counts.characters();
    auto const seq_count = counts.added_sequences_count();
    auto const num_chars = counts.characters().size();

    // Use a hard coded gap char here, as we have fixed character codes anyway.
    auto const gap_char  = '-';

    // We expect ACGT here.
    assert( num_chars == 4 );

    // Special case: empty counts object. In this case, return an all-gap sequence.
    if( seq_count == 0 ) {
        // We do not immediately return the string here, in order to facilitate copy elision.
        result = std::string( counts.length(), gap_char );
        return result;
    }

    // Process all sites of the sequence.
    for( size_t site_idx = 0; site_idx < counts.length(); ++site_idx ) {

        // Total sum of counts. Used for getting the number of gaps.
        size_t counts_sum = 0;

        // Map from counts to positions. We use this for sorting by count. It's a vector, because
        // it will only have 4 or 5 elements, thus this should be faster than complex containers.
        std::vector< CountPair > counts_map;
        counts_map.reserve( 4 );

        // Add all chars with their counts to the map.
        for( size_t char_idx = 0; char_idx < num_chars; ++char_idx ) {
            auto const char_count = counts.count_at( char_idx, site_idx );
            counts_sum += char_count;
            counts_map.push_back({ char_count, chars[ char_idx ] });
        }

        // We can never have a sum of counts higher then the number of sequences that were added
        // to the counts object.
        assert( counts_sum <= seq_count );

        // We already checked that the counts object is for nucleic acids.
        // Thus, we expect four values here.
        assert( counts_map.size() == 4 );

        // If we want to use gaps as a normal character, add their count to the map.
        // We want to compare the gap count with all other frequencies,
        // so it makes sense to just treat it as a normal character here.
        // In the char_selector function, some special care might need to be taken however.
        if( allow_gaps ) {
            auto const gap_count = seq_count - counts_sum;
            counts_sum += gap_count;
            counts_map.push_back({ gap_count, gap_char });

            // Now that we added gaps, the total sum matches the number of added sequences.
            assert( counts_sum == seq_count );
        }

        // Sort the counts so that the highest one is first.
        std::sort( counts_map.begin(), counts_map.end(), CountPairComparator );

        // Get the ambiguity code that represents the selected characters and add it to the seq.
        result += char_selector( counts_map, counts_sum );
    }

    return result;
}

// =================================================================================================
//     Majority
// =================================================================================================

std::string consensus_sequence_with_majorities(
    SiteCounts const&     counts,
    bool                  allow_gaps,
    char                  gap_char
) {
    std::string result;
    result.reserve( counts.length() );

    // Prepare some constants for simplicity.
    auto const chars     = counts.characters();
    auto const seq_count = counts.added_sequences_count();
    auto const num_chars = counts.characters().size();

    for( size_t site_idx = 0; site_idx < counts.length(); ++site_idx ) {

        size_t                        max_pos    = 0;
        SiteCounts::CountsIntType max_val    = 0;
        SiteCounts::CountsIntType counts_sum = 0;

        for( size_t char_idx = 0; char_idx < num_chars; ++char_idx ) {
            auto const char_count = counts.count_at( char_idx, site_idx );
            counts_sum += char_count;

            // We use a strict greater here, as this ensures to use the first character in cases
            // where many have the same count.
            if( char_count > max_val ) {
                max_pos = char_idx;
                max_val = char_count;
            }
        }

        // We can never have a max higher than the total sum of counts, and this again cannot be
        // higher then the number of sequences that were added to the counts object.
        assert( max_val    <= counts_sum );
        assert( counts_sum <= seq_count  );

        // We write a code char if it is the majority, that is, > 0 and > all other code counts.
        // In other cases, write a gap. That is, either no code has a count > 0, or, if we allow
        // gaps and gaps are more frquent than actual codes.
        auto const gap_count = seq_count - counts_sum;
        if(( max_val > 0 ) && (( ! allow_gaps ) || ( max_val > gap_count ))) {
            result += chars[ max_pos ];
        } else {
            result += gap_char;
        }
    }

    return result;
}

std::string consensus_sequence_with_majorities(
    SequenceSet const&    sequences,
    std::string const&    characters,
    bool                  allow_gaps,
    char                  gap_char
) {
    // Basic checks.
    if( sequences.size() == 0 ) {
        throw std::runtime_error( "Cannot calculate consensus sequence of empty SequenceSet." );
    }
    if( ! is_alignment( sequences ) ) {
        throw std::runtime_error(
            "Cannot calculate consensus sequence for SequenceSet that is not an alignment. "
            "That is, all Sequences need to have the same length."
        );
    }

    // Build counts object.
    auto counts = SiteCounts( characters, sequences[0].size() );
    counts.add_sequences( sequences );

    // Return consensus sequence.
    return consensus_sequence_with_majorities( counts, allow_gaps, gap_char );
}

std::string consensus_sequence_with_majorities(
    SequenceSet const&    sequences,
    bool                  allow_gaps
) {
    return consensus_sequence_with_majorities(
        sequences,
        nucleic_acid_codes_plain(),
        allow_gaps,
        '-'
    );
}

// =================================================================================================
//     Ambiguity
// =================================================================================================

std::string consensus_sequence_with_ambiguities(
    SiteCounts const&     counts,
    double                similarity_factor,
    bool                  allow_gaps
) {
    // Check the deviation range.
    if( similarity_factor < 0.0 || similarity_factor > 1.0 ) {
        throw std::invalid_argument(
            "Value of similarity_factor has to be in range [ 0.0, 1.0 ]."
        );
    }

    // Use a hard coded gap char here, as we have fixed character codes anyway.
    auto const gap_char  = '-';

    // Functor that selects the chars according to the consensus specification.
    auto selector = [&](
        std::vector<std::pair<size_t, char>> const& counts_map,
        size_t const                                counts_sum
    ){
        // So that everyone knows what we are dealing with.
        assert( std::is_sorted( counts_map.begin(), counts_map.end(), CountPairComparator ) );

        // Special case. All gap site, but allow_gaps == false.
        // In this case, return a gap instead of an 'N'.
        if( counts_sum == 0 ) {
            assert( ! allow_gaps );
            return gap_char;
        }

        // Special case. If gap is the most frequent char, we just return it.
        if( counts_map[0].second == gap_char ) {
            return gap_char;
        }

        // Prepare a string of characters codes for the ambiguities, init with the most frequent char.
        auto ambiguity_codes = std::string( 1, counts_map[0].second );

        // Every character that has at least this count is added to the ambiguity.
        auto const deviation_threshold = similarity_factor * static_cast<double>( counts_map[0].first );

        // Compare the less frequent codes to the most frequent one and
        // decide whether to add them to the ambiguities.
        for( size_t i = 1; i < counts_map.size(); ++i ) {
            auto const cur_count = static_cast<double>( counts_map[i].first );

            // If the count is below the threshold, we are done.
            // The map is sorted, so no other count will be high enough.
            // We also avoid zero counts, as this leads to wrong results with an
            // similarity_factor of 0.0. It would then just add all, ending up with all "N"s,
            // instead of just all codes that appear in the sequence.
            if( cur_count < deviation_threshold || cur_count == 0.0 ) {
                break;
            }

            // If it is a gap, we are done - the result is a gap, too.
            // If not, add it to the ambiguities.
            if( counts_map[i].second == gap_char ) {
                return gap_char;
            } else {
                ambiguity_codes += counts_map[i].second;
            }
        }

        // Return the ambiguity code that represents the selected characters.
        return nucleic_acid_ambiguity_code( ambiguity_codes );
    };

    return consensus_sequence_template( counts, allow_gaps, selector );
}

std::string consensus_sequence_with_ambiguities(
    SequenceSet const&    sequences,
    double                similarity_factor,
    bool                  allow_gaps
) {
    // Basic checks.
    if( sequences.size() == 0 ) {
        throw std::runtime_error( "Cannot calculate consensus sequence of empty SequenceSet." );
    }
    if( ! is_alignment( sequences ) ) {
        throw std::runtime_error(
            "Cannot calculate consensus sequence for SequenceSet that is not an alignment. "
            "That is, all Sequences need to have the same length."
        );
    }

    // Build counts object.
    auto counts = SiteCounts( nucleic_acid_codes_plain(), sequences[0].size() );
    counts.add_sequences( sequences );

    // Return consensus sequence.
    return consensus_sequence_with_ambiguities( counts, similarity_factor, allow_gaps );
}

// =================================================================================================
//     Threshold
// =================================================================================================

std::string consensus_sequence_with_threshold(
    SiteCounts const&     counts,
    double                frequency_threshold,
    bool                  allow_gaps,
    bool                  use_ambiguities
) {
    // Check the frequency threshold.
    if( frequency_threshold < 0.0 || frequency_threshold > 1.0 ) {
        throw std::invalid_argument(
            "Value of frequency_threshold has to be in range [ 0.0, 1.0 ]."
        );
    }

    // Use hard coded chars here, as we have fixed character codes anyway.
    auto const gap_char  = '-';
    auto const mask_char = 'X';

    // Functor that selects the chars according to the consensus specification.
    auto selector = [&](
        std::vector<std::pair<size_t, char>> const& counts_map,
        size_t const                                counts_sum
    ){
        // So that everyone knows what we are dealing with.
        assert( std::is_sorted( counts_map.begin(), counts_map.end(), CountPairComparator ) );

        // Special case. All gap site, but allow_gaps == false.
        // In this case, return a gap instead of an 'N'.
        if( counts_sum == 0 ) {
            assert( ! allow_gaps );
            return gap_char;
        }

        // Prepare a string of characters codes for the ambiguities.
        std::string ambiguity_codes;

        // Add up the counts and combine ambiguities until we reach the threshold.
        // If we still do not reach the threshold with all codes, we end up with an N.
        size_t accumulated_sum = 0;
        double fraction        = 0.0;
        for( size_t i = 0; i < counts_map.size(); ++i ) {

            // If there are no counts, we do not use it (and stop here, because in a sorted
            // counts order, all following counts will be zero anyway). This way, we only use
            // those codes for the ambiguity that actually appear at the site.
            if( counts_map[i].first == 0 ) {
                break;
            }

            // If it is a gap, we are done - the result is a gap, too.
            if( counts_map[i].second == gap_char ) {
                return gap_char;
            }

            // Use this char!
            accumulated_sum += counts_map[i].first;
            ambiguity_codes += counts_map[i].second;

            // Check if we already reached the threshold.
            // The division is okay, as we already checked that counts_sum > 0 before.
            fraction = static_cast<double>( accumulated_sum ) / static_cast<double>( counts_sum );
            if( fraction >= frequency_threshold ) {
                break;
            }
        }

        // We checked that counts_sum > 0 in the beginning. Thus, at least counts_map needs to
        // contain non-zero entries. Thus, we added at least one char to ambiguity_codes.
        assert( ambiguity_codes.size() > 0 );

        // Finally, return the needed code.
        if( ambiguity_codes.size() > 1 && ! use_ambiguities ) {
            return mask_char;
        } else {
            return nucleic_acid_ambiguity_code( ambiguity_codes );
        }
    };

    return consensus_sequence_template( counts, allow_gaps, selector );
}

std::string consensus_sequence_with_threshold(
    SequenceSet const&    sequences,
    double                frequency_threshold,
    bool                  allow_gaps,
    bool                  use_ambiguities
) {
    // Basic checks.
    if( sequences.size() == 0 ) {
        throw std::runtime_error( "Cannot calculate consensus sequence of empty SequenceSet." );
    }
    if( ! is_alignment( sequences ) ) {
        throw std::runtime_error(
            "Cannot calculate consensus sequence for SequenceSet that is not an alignment. "
            "That is, all Sequences need to have the same length."
        );
    }

    // Build counts object.
    auto counts = SiteCounts( nucleic_acid_codes_plain(), sequences[0].size() );
    counts.add_sequences( sequences );

    // Return consensus sequence.
    return consensus_sequence_with_threshold(
        counts,
        frequency_threshold,
        allow_gaps,
        use_ambiguities
     );
}

// =================================================================================================
//     Cavener
// =================================================================================================

std::string consensus_sequence_cavener(
    SiteCounts const&     counts,
    bool                  allow_gaps
) {
    // Functor that selects the chars according to the consensus specification.
    auto selector = [&](
        std::vector<std::pair<size_t, char>> const& counts_map,
        size_t const                                counts_sum
    ){
        // So that everyone knows what we are dealing with.
        assert( std::is_sorted( counts_map.begin(), counts_map.end(), CountPairComparator ) );

        // Special case. All gap site, but allow_gaps == false.
        // In this case, return a gap instead of an 'N'.
        if( counts_sum == 0 ) {
            assert( ! allow_gaps );
            return '-';
        }

        // Prepare
        std::string ambiguity_codes;

        // If the hightest freq is > 50% and > 2 * second highest freq, just use it.
        if(( 2 * counts_map[0].first > counts_sum ) && ( counts_map[0].first > 2 * counts_map[1].first )) {
            ambiguity_codes = std::string( 1, counts_map[0].second );

        // If the first two freqs > 75% (and both < 50%, which was checked above), use dual code.
        } else if( counts_map[0].first + counts_map[1].first > 3.0 * static_cast<double>( counts_sum ) / 4.0 ) {
            ambiguity_codes
                = std::string( 1, counts_map[0].second )
                + std::string( 1, counts_map[1].second )
            ;

        // If neither of the above, but one freq is 0, then use three codes.
        } else if( counts_map[3].first == 0 ) {
            ambiguity_codes
                = std::string( 1, counts_map[0].second )
                + std::string( 1, counts_map[1].second )
                + std::string( 1, counts_map[2].second )
            ;

        // Fall back case: Use 'N'.
        } else {
            ambiguity_codes = "ACGT";
        }

        // So far, we have treated gap chars as any other. As gaps are not mentioned in the
        // original method, this is the best we can do.
        // So now, if we have a gap in there, we return a gap as the end result,
        // as combining a gap with anything else results in a gap.
        std::sort( ambiguity_codes.begin(), ambiguity_codes.end() );
        if( ambiguity_codes.size() > 0 && ambiguity_codes[0] == '-' ) {
            return '-';
        }

        // Return the ambiguity code that represents the selected characters.
        return nucleic_acid_ambiguity_code( ambiguity_codes );
    };

    return consensus_sequence_template( counts, allow_gaps, selector );
}

std::string consensus_sequence_cavener(
    SequenceSet const&    sequences,
    bool                  allow_gaps
) {
    // Basic checks.
    if( sequences.size() == 0 ) {
        throw std::runtime_error( "Cannot calculate consensus sequence of empty SequenceSet." );
    }
    if( ! is_alignment( sequences ) ) {
        throw std::runtime_error(
            "Cannot calculate consensus sequence for SequenceSet that is not an alignment. "
            "That is, all Sequences need to have the same length."
        );
    }

    // Build counts object.
    auto counts = SiteCounts( nucleic_acid_codes_plain(), sequences[0].size() );
    counts.add_sequences( sequences );

    // Return consensus sequence.
    return consensus_sequence_cavener( counts, allow_gaps );
}

} // namespace sequence
} // namespace genesis
