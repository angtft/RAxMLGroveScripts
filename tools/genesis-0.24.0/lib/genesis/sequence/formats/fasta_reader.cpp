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
 * @ingroup sequence
 */

#include "genesis/sequence/formats/fasta_reader.hpp"

#include "genesis/sequence/functions/labels.hpp"
#include "genesis/sequence/sequence_set.hpp"
#include "genesis/sequence/sequence.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/io/char.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/io/scanner.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace genesis {
namespace sequence {

// =================================================================================================
//     Constructor and Rule of Five
// =================================================================================================

FastaReader::FastaReader()
{
    lookup_.set_all( true );
}

// =================================================================================================
//     Reading
// =================================================================================================

SequenceSet FastaReader::read( std::shared_ptr< utils::BaseInputSource > source ) const
{
    SequenceSet result;
    utils::InputStream is( source );
    parse_document( is, result );
    return result;
}

void FastaReader::read(
    std::shared_ptr< utils::BaseInputSource > source,
    SequenceSet& sequence_set
) const {
    utils::InputStream is( source );
    parse_document( is, sequence_set );
}

// =================================================================================================
//     Parsing
// =================================================================================================

void FastaReader::parse_document(
    utils::InputStream& input_stream,
    SequenceSet&        sequence_set
) const {
    Sequence seq;

    if( parsing_method_ == ParsingMethod::kDefault ) {
        while( parse_sequence( input_stream, seq ) ) {
            sequence_set.add( seq );
        }

    } else if( parsing_method_ == ParsingMethod::kPedantic ) {
        while( parse_sequence_pedantic( input_stream, seq ) ) {
            sequence_set.add( seq );
        }

    } else {
        // There are no other methods currently implemented.
        assert( false );
    }
}

bool FastaReader::parse_sequence(
    utils::InputStream& input_stream,
    Sequence&           sequence
) const {
    // Init. Call clear in order to avoid not setting properties that might be added to
    // Sequence in the future. Should not noticeable affect speed, as the sequence string capacities
    // should not change when setting the strings to empty strings.
    auto& it = input_stream;
    sequence.clear();

    // Check for data.
    if( !it ) {
        return false;
    }

    // ---------------------------------------------
    //     Label
    // ---------------------------------------------

    // Scope to ensure that the label line is only used
    // while we actually are in that line.
    {

    // Check beginning of sequence.
    if( !it || *it != '>' ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting '>' at beginning of sequence at line " + std::to_string( it.line() ) + "."
        );
    }
    assert( it && *it == '>' );
    ++it;

    // Parse label.
    std::string const label = utils::read_while( it, isprint );
    if( label == "" ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting label after '>' in sequence at line " + std::to_string( it.line() ) + "."
        );
    }
    if( guess_abundances_ ) {
        auto const la = guess_sequence_abundance( label );
        sequence.label( la.first );
        sequence.abundance( la.second );
    } else {
        sequence.label( label );
    }

    // Check for unexpected end of file.
    if( !it || *it != '\n' ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Unexpected characters at the end of the label line in sequence at line "
            + std::to_string( it.line() ) + "."
        );
    }
    assert( it && *it == '\n' );
    ++it;

    } // End of line scope. We are done with the label line.

    // ---------------------------------------------
    //     Sites
    // ---------------------------------------------

    // Skip comments.
    while( it && *it == ';' ) {
        utils::skip_until( it, '\n' );
        assert( it && *it == '\n' );
        ++it;
    }

    // Check for unexpected end of file.
    if( !it ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting a sequence after the label line in sequence at line "
            + std::to_string( it.line() - 1 ) + "."
        );
    }
    assert( it );

    // Reserve some tmp memory. We will later copy the content, so that superfluous capacity
    // is stripped.
    // We could do a sites.reserve( ... ) here, but this yields only minor speedups.
    std::string sites;
    // sites.reserve( n );

    // Parse sequence. At every beginning of the loop, we are at a line start.
    while( it && *it != '>' ) {
        assert( it.column() == 1 );
        it.get_line( sites );
    }
    assert( !it || *it == '>' );

    if( sites.length() == 0 ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name() + ": Empty sequence at line "
            + std::to_string( it.line() - 1 ) + "."
        );
    }

    if( site_casing_ == SiteCasing::kToUpper ) {
        sequence.sites() = utils::to_upper_ascii( sites );
    } else if( site_casing_ == SiteCasing::kToLower ) {
        sequence.sites() = utils::to_lower_ascii( sites );
    } else {
        // We could do a move here instead, but this way, we save some memory, which might be
        // more reasonable for big sequences files than the small gain in speed.
        // sequence.sites() = std::move(sites);
        sequence.sites() = sites;
    }

    if( use_validation_ ) {
        for( auto const& c : sequence.sites() ) {
            if( !lookup_[c] ) {
                throw std::runtime_error(
                    "Malformed Fasta " + it.source_name() + ": Invalid sequence symbol "
                    + utils::char_to_hex( c )
                    + " in the sequence at/above line " + std::to_string( it.line() - 1 ) + "."
                );
            }
        }
    }

    return true;
}

bool FastaReader::parse_sequence_pedantic(
    utils::InputStream& input_stream,
    Sequence&           sequence
) const {
    // Init. Call clear in order to avoid not setting properties that might be added to
    // Sequence in the future. Should not noticeable affect speed, as the sequence string capacities
    // should not change when setting the strings to empty strings.
    auto& it = input_stream;
    sequence.clear();

    // Check for data.
    if( !it ) {
        return false;
    }

    // Check beginning of sequence.
    if( it.current() != '>' ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting '>' at beginning of sequence at " + it.at() + "."
        );
    }
    assert( it && *it == '>' );
    ++it;

    // Parse label.
    std::string label = utils::read_while( it, isprint );
    if( label == "" ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting label after '>' at " + it.at() + "."
        );
    }
    if( guess_abundances_ ) {
        auto const la = guess_sequence_abundance( label );
        sequence.label( la.first );
        sequence.abundance( la.second );
    } else {
        sequence.label( label );
    }

    // Check for unexpected end of stream.
    if( !it || ( *it != '\n' )) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting a sequence after the label line at " + it.at() + "."
        );
    }
    assert( it && (*it == '\n' ));

    // Check for unexpected end of file.
    if( !it || *it != '\n' ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting a sequence after the label line at " + it.at() + "."
        );
    }
    assert( it && *it == '\n' );

    // Skip comments.
    while( it && *it == ';' ) {
        utils::skip_while( it, isprint );
    }

    // Check for unexpected end of file.
    if( !it || *it != '\n' ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Expecting a sequence after the label line at " + it.at() + "."
        );
    }
    assert( it && *it == '\n' );
    ++it;

    // Parse sequence. At every beginning of the outer loop, we are at a line start.
    std::string sites;
    while( it && *it != '>' ) {
        assert( it.column() == 1 );

        size_t count = 0;
        while( it && *it != '\n' ) {

            // Weird C relicts need weird conversions...
            // See https://en.cppreference.com/w/cpp/string/byte/tolower
            char c = *it;
            if( site_casing_ == SiteCasing::kToUpper ) {
                c = static_cast<char>( std::toupper( static_cast<unsigned char>( c )));
            } else if( site_casing_ == SiteCasing::kToLower ) {
                c = static_cast<char>( std::tolower( static_cast<unsigned char>( c )));
            }
            if( use_validation_ && ! lookup_[c] ) {
                throw std::runtime_error(
                    "Malformed Fasta " + it.source_name() + ": Invalid sequence symbol "
                    + utils::char_to_hex( c ) + " in sequence at " + it.at() + "."
                );
            }

            sites += c;
            ++it;
            ++count;
        }

        if( count == 0 ) {
            throw std::runtime_error(
                "Malformed Fasta " + it.source_name()
                + ": Empty sequence line at " + it.at() + "."
            );
        }

        if( !it ) {
            throw std::runtime_error(
                "Malformed Fasta " + it.source_name()
                + ": Sequence line does not end with '\\n' at " + it.at() + "."
            );
        }
        assert( it && *it == '\n' );
        ++it;
    }
    assert( !it || *it == '>' );

    if( sites.length() == 0 ) {
        throw std::runtime_error(
            "Malformed Fasta " + it.source_name()
            + ": Empty sequence at " + it.at() + "."
        );
    }

    // Copy the sequence. We do not use move here, as we can save some memory this way.
    sequence.sites() = sites;

    return true;
}

// =================================================================================================
//     Properties
// =================================================================================================

FastaReader& FastaReader::parsing_method( FastaReader::ParsingMethod value )
{
    parsing_method_ = value;
    return *this;
}

FastaReader::ParsingMethod FastaReader::parsing_method() const
{
    return parsing_method_;
}

FastaReader& FastaReader::site_casing( SiteCasing value )
{
    site_casing_ = value;
    return *this;
}

FastaReader::SiteCasing FastaReader::site_casing() const
{
    return site_casing_;
}

FastaReader& FastaReader::guess_abundances( bool value )
{
    guess_abundances_ = value;
    return *this;
}

bool FastaReader::guess_abundances() const
{
    return guess_abundances_;
}

FastaReader& FastaReader::valid_chars( std::string const& chars )
{
    if( chars.size() == 0 ) {
        lookup_.set_all( true );
        use_validation_ = false;
    } else {
        lookup_.set_all( false );
        lookup_.set_selection( chars, true );
        use_validation_ = true;
    }

    return *this;
}

std::string FastaReader::valid_chars() const
{
    // We need to check the valid chars lookup here, because we don't want to return a string
    // of _all_ chars.
    if( ! use_validation_ || lookup_.all_equal_to( true ) ) {
        return "";
    } else {
        return lookup_.get_chars_equal_to( true );
    }
}

utils::CharLookup<bool>& FastaReader::valid_char_lookup()
{
    return lookup_;
}

} // namespace sequence
} // namespace genesis
