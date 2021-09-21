#ifndef GENESIS_SEQUENCE_FUNCTIONS_LABELS_H_
#define GENESIS_SEQUENCE_FUNCTIONS_LABELS_H_

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

#include "genesis/utils/tools/hash/functions.hpp"

#include <string>
#include <utility>
#include <unordered_map>
#include <unordered_set>

namespace genesis {
namespace sequence {

// =================================================================================================
//     Forwad Declarations
// =================================================================================================

class Sequence;
class SequenceSet;

// =================================================================================================
//     Helper Structs
// =================================================================================================

struct LabelAttributes
{
    std::string label;
    std::unordered_map<std::string, std::string> attributes;
};

// =================================================================================================
//     General
// =================================================================================================

/**
 * @brief Return a pointer to a Sequence with a specific label, or `nullptr` iff not found.
 */
Sequence const* find_sequence( SequenceSet const& set, std::string const& label );

/**
 * @brief Return a set of all labels of the SequenceSet.
 */
std::unordered_set<std::string> labels( SequenceSet const& set );

/**
 * @brief Guess the abundance of a Sequence, using it's label.
 *
 * The function splits the label of a Sequence into two parts: the descriptive name of the
 * sequence, and an abundance value (weight or multiplicity of the sequence),
 * which are returned as a `std::pair`.
 *
 * The function accepts two patterns of reporting abundances via the
 * @link Sequence::label() label()@endlink of a Sequence:
 *
 *  * Appended via underscore: `name_123`. In this case, the number has to be the last in
 *    the label, that is, no other text may follow.
 *  * Using the format `;size=123;`. The semicoli are mandatory, except the second one if nothing
 *    else follows in the label. See label_attributes() for details.
 *
 * If neither of them is found, a default abundance of 1 is returned.
 */
std::pair<std::string, size_t> guess_sequence_abundance( Sequence const& sequence );

/**
 * @brief Guess the abundance of a Sequence, given it's label.
 *
 * This is the same as guess_sequence_abundance( Sequence const& ), but takes the label
 * as a string, instead of the Sequence object. See there for details.
 */
std::pair<std::string, size_t> guess_sequence_abundance( std::string const& label );

/**
 * @brief Get the attributes list (semicolons-separated) from a Sequence.
 *
 * It is common to store additional information in sequence headers, e.g., in the `fasta` format,
 * using a semicolon-separated list of attributes like this:
 *
 *     >some_name;size=123;thing=foo;
 *
 * This function disects this kind of information and returns it.
 * The returned struct contains the label (the part before the first semicolon),
 * as well as a map for the attributes. As this is not a multimap, later attributes with the same
 * key overwrite earlier ones.
 *
 * If the sequence label does not contain any information that is separated via a semicolon,
 * the attributes list is returned empty. However, if semicola are found in the label,
 * the correct format is expected (with the syntax `;key=value;`) for each attribute.
 * Otherwise, an exception is thrown. The last semicolon is optional; that is, the label
 * can simply end after the last value.
 */
LabelAttributes label_attributes( Sequence const& sequence );

/**
 * @brief Get the attributes list (semicolons-separated) from a Sequence, given it's label.
 *
 * This is the same as label_attributes( Sequence const& ), but takes the label
 * as a string, instead of the Sequence object. See there for details.
 */
LabelAttributes label_attributes( std::string const& label );

// =================================================================================================
//     Uniqueness
// =================================================================================================

/**
 * @brief Return true iff all labels of the Sequence%s in the SequenceSet are unique.
 *
 * The optional parameter `case_sensitive` controls how labels are compared. Default is `true`,
 * that is, Sequence%s are compared case-sensitively.
 */
bool has_unique_labels( SequenceSet const& set, bool case_sensitive = true );

/**
 * @brief Relabel the Sequence using the hash digest of its sites.
 *
 * See ::utils::HashingFunctions for the available hashing functions.
 */
void relabel_with_hash( Sequence& seq, utils::HashingFunctions hash_function );

/**
 * @brief Relabel all Sequence%s in the SequenceSet using the hash digest of the sites.
 *
 * See ::utils::HashingFunctions for the available hashing functions.
 *
 * If there are duplicate Sequence%s, this function will lead to multiple Sequence%s with the same
 * name, which might be an issue for downstream programs that expect unique labels.
 * See has_unique_labels() to check this.
 */
void relabel_with_hash( SequenceSet& set, utils::HashingFunctions hash_function );

// =================================================================================================
//     Validity
// =================================================================================================

/**
 * @brief Check whether a given string is a valid label for a Sequence.
 *
 * While we can work with any form of label (as long as it is a string), most file formats and
 * consequently most programs that read them restrict the set of valid characters for labels of
 * sequences. We thus provide this function, which uses the most common interpretation of valid
 * labels.
 *
 * A label is valid if its characters have a graphical representation (i.e., isgraph() is true) and
 * if non of these characters occurs:
 *
 *     :,();[]'
 *
 * Thus, all whitespaces, control characters, and the listed special characters are invalid.
 * See sanitize_label() for a function that replaces all invalid characters of the label by
 * underscores.
 */
bool is_valid_label(   std::string const& label );

/**
 * @brief Check whether a Sequence has a valid label.
 *
 * This might be important for printing the Sequence to a file that needs to be read by other
 * applications. See is_valid_label() for details on what is considered a valid label.
 * See sanitize_label() for a function that replaces all invalid characters of the label by
 * underscores.
 */
bool has_valid_label(  Sequence const&    seq );

/**
 * @brief Check whether all Sequence%s in a SequenceSet have valid labels.
 *
 * This might be important for printing the Sequences to a file that needs to be read by other
 * applications. See is_valid_label() for details on what is considered a valid label.
 * See sanitize_labels() for a function that replaces all invalid characters of the labels by
 * underscores.
 */
bool has_valid_labels( SequenceSet const& set );

/**
 * @brief Sanitize a label by replacing all invalid characters with underscores.
 *
 * This might be important for printing the Sequences to a file that needs to be read by other
 * applications. See is_valid_label() for details on what is considered a valid label.
 */
std::string sanitize_label( std::string const& label );

/**
 * @brief Sanitize a label by replacing all invalid characters with underscores.
 *
 * This might be important for printing the Sequences to a file that needs to be read by other
 * applications. See is_valid_label() for details on what is considered a valid label.
 */
void sanitize_label( Sequence&     seq );

/**
 * @brief Sanitize the labels of all Sequence%s in the SequenceSet by replacing all invalid
 * characters with underscores.
 *
 * This might be important for printing the Sequences to a file that needs to be read by other
 * applications. See is_valid_label() for details on what is considered a valid label.
 */
void sanitize_labels( SequenceSet& set );

// =================================================================================================
//     Modifiers
// =================================================================================================

/**
 * @brief Remove all those Sequence%s from a SequenceSet whose labels are in the given list.
 *
 * If `invert` is set to true, it does the same inverted: it removes all Sequence%s except those
 * whose label is in the list.
 */
void filter_by_label_list(
    SequenceSet&                           set,
    std::unordered_set<std::string> const& labels,
    bool                                   invert = false
);

} // namespace sequence
} // namespace genesis

#endif // include guard
