#ifndef GENESIS_SEQUENCE_FORMATS_FASTA_READER_H_
#define GENESIS_SEQUENCE_FORMATS_FASTA_READER_H_

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

#include "genesis/utils/tools/char_lookup.hpp"
#include "genesis/utils/io/input_source.hpp"

#include <iosfwd>
#include <memory>
#include <string>

namespace genesis {

// =================================================================================================
//     Forward Declarations
// =================================================================================================

namespace utils {
    class InputStream;
}

namespace sequence {
    class SequenceSet;
    class Sequence;
}

// =================================================================================================
//     Fasta Reader
// =================================================================================================

namespace sequence {

/**
 * @brief Read Fasta sequence data.
 *
 * This class provides simple facilities for reading Fasta data.
 *
 * Exemplary usage:
 *
 *     std::string infile = "path/to/file.fasta";
 *     SequenceSet sequence_set;
 *
 *     FastaReader()
 *         .site_casing( SiteCasing::kUnchanged )
 *         .valid_chars( nucleic_acid_codes_all() )
 *         .read( utils::from_file( infile ), sequence_set );
 *
 * The expected data format:
 *
 *   1. Has to start with a '>' character, followed by a label, ended by a '\\n'.
 *   2. An arbitrary number of comment lines, starting with ';', can follow, but are ignored.
 *   3. After that, a sequence has to follow, over one or more lines.
 *
 * More information on the format can be found at:
 *
 *    * http://en.wikipedia.org/wiki/FASTA_format
 *    * http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
 *    * http://zhanglab.ccmb.med.umich.edu/FASTA/
 *
 * Using site_casing(), the sequences can automatically be turned into upper or lower case letter.
 * Also, see valid_chars( std::string const& chars ) for a way of checking correct input sequences.
 */
class FastaReader
{
public:

    // ---------------------------------------------------------------------
    //     Typedefs and Enums
    // ---------------------------------------------------------------------

    /**
     * @brief Enumeration of the available methods for parsing Fasta sequences.
     */
    enum class ParsingMethod
    {
        /**
         * @brief Fast method, used by default.
         *
         * This is by far the preferred method, it however has one slight limitation:
         * It only reports errors using the line where the sequence starts.
         * This does not affect most applications, as good data won't produce errors to
         * report. If you however want error reporting at the
         * exact line and column where the error occurs, use kPedantic instead.
         *
         * With this setting, parse_sequence() is used for parsing.
         * In our tests, it achieved ~350 MB/s parsing speed.
         */
        kDefault,

        /**
         * @brief Pedantic method.
         *
         * Compared to the fast method, this one reports errors at the exact line and column
         * where they occur. It is however slower (~3.5x the time of the default method).
         * Apart from that, there are no differences.
         *
         * If you need this method for certain files, it might be useful to use it only once and
         * use a FastaWriter to write out a new Fasta file without errors.
         * This way, for subsequent reading you can then use the faster default method.
         *
         * With this setting, parse_sequence_pedantic() is used for parsing.
         * In our tests, it achieved ~100 MB/s parsing speed.
         */
        kPedantic
    };

    /**
     * @brief Enumeration of casing methods to apply to each site of a Sequence.
     */
    enum class SiteCasing
    {
        /**
         * @brief Do not change the case of the sites.
         */
        kUnchanged,

        /**
         * @brief Make all sites upper case.
         */
        kToUpper,

        /**
         * @brief Make all sites lower case.
         */
        kToLower
    };

    // ---------------------------------------------------------------------
    //     Constructor and Rule of Five
    // ---------------------------------------------------------------------

    /**
     * @brief Create a default FastaReader. Per default, chars are turned upper case, but not validated.
     *
     * See site_casing() and valid_chars() to change this behaviour.
     */
    FastaReader();
    ~FastaReader() = default;

    FastaReader( FastaReader const& ) = default;
    FastaReader( FastaReader&& )      = default;

    FastaReader& operator= ( FastaReader const& ) = default;
    FastaReader& operator= ( FastaReader&& )      = default;

    // ---------------------------------------------------------------------
    //     Reading
    // ---------------------------------------------------------------------

    /**
     * @brief Read all Sequence%s from an input source in Fasta format and return them as a
     * SequenceSet.
     *
     * Use functions such as utils::from_file() and utils::from_string() to conveniently
     * get an input source that can be used here.
     */
    SequenceSet read( std::shared_ptr< utils::BaseInputSource > source ) const;

    /**
     * @brief Read all Sequence%s from an input source in Fasta format into a SequenceSet.
     *
     * The Sequences are added to the SequenceSet, whose existing Sequences are kept. Thus, by
     * repeatedly calling this or similar read functions, multiple input files can easily be read
     * into one SequenceSet.
     *
     * Use functions such as utils::from_file() and utils::from_string() to conveniently
     * get an input source that can be used here.
     */
    void read( std::shared_ptr< utils::BaseInputSource > source, SequenceSet& sequence_set ) const;

    // ---------------------------------------------------------------------
    //     Parsing
    // ---------------------------------------------------------------------

    /**
     * @brief Parse a whole fasta document into a SequenceSet.
     *
     * This function is mainly used internally by the reading functions read().
     * It uses the currently set parsing_method() for parsing the data.
     */
    void parse_document(
        utils::InputStream& input_stream,
        SequenceSet&        sequence_set
    ) const;

    /**
     * @brief Parse a Sequence in Fasta format.
     *
     * This function takes an InputStream and interprets it as a Fasta formatted sequence. It extracts
     * the data and writes it into the given Sequence object. See the class description of FastaReader
     * for the expected data format.
     *
     * The function stops after parsing one such sequence. It returns true if a sequence was extracted
     * and false if the stream is empty. If the input is not in the correct format, an
     * `std::runtime_error` exception is thrown indicating the malicious position in the input stream.
     */
    bool parse_sequence(
        utils::InputStream& input_stream,
        Sequence&           sequence
    ) const;

    /**
     * @brief Parse a Sequence in Fasta format.
     *
     * This function takes an InputStream and interprets it as a Fasta formatted sequence. It extracts
     * the data and writes it into the given Sequence object. See the class description of FastaReader
     * for the expected data format.
     *
     * The function stops after parsing one such sequence. It returns true if a sequence was extracted
     * and false if the stream is empty. If the input is not in the correct format, an
     * `std::runtime_error` exception is thrown indicating the malicious position in the input stream.
     *
     * Compared to parse_sequence(), this function reports errors at the exact line and column
     * where they occur. It is however slower.
     * Apart from that, there are no differences. See FastaReader::ParsingMethod for details.
     */
    bool parse_sequence_pedantic(
        utils::InputStream& input_stream,
        Sequence&           sequence
    ) const;

    // ---------------------------------------------------------------------
    //     Properties
    // ---------------------------------------------------------------------

    /**
     * @brief Set the parsing method.
     *
     * The parsing method is used for all the reader functions and parse_document().
     * See the @link FastaReader::ParsingMethod ParsingMethod enum@endlink for details.
     */
    FastaReader&  parsing_method( ParsingMethod value );

    /**
     * @brief Return the currently set parsing method.
     *
     * See the @link FastaReader::ParsingMethod ParsingMethod enum@endlink for details.
     */
    ParsingMethod parsing_method() const;

    /**
     * @brief Set whether Sequence sites are automatically turned into upper or lower case.
     *
     * Default is SiteCasing::kToUpper, that is, all sites of the read Sequences are turned into
     * upper case letters automatically. This is demanded by the Fasta standard.
     * The function returns the FastaReader object to allow for fluent interfaces.
     */
    FastaReader& site_casing( SiteCasing value );

    /**
     * @brief Return whether Sequence sites are automatically turned into upper or lower case.
     */
    SiteCasing site_casing() const;

    /**
     * @brief Set whether Sequence labels are used to guess/extract Sequence abundances
     *
     * Default is `false`, that is, labels are just taken as they are in the input.
     * If set to `true`, the label is used to guess an abundance count, which is set in the Sequence.
     * See guess_sequence_abundance( Sequence const& ) for the valid formats of such abundances.
     */
    FastaReader& guess_abundances( bool value );

    /**
     * @brief Return whether the label is used to guess/extracat Sequence abundances.
     */
    bool guess_abundances() const;

    /**
     * @brief Set the chars that are used for validating Sequence sites when reading them.
     *
     * When this function is called with a string of chars, those chars are used to validate the
     * sites when reading them. That is, only sequences consisting of the given chars are valid.
     *
     * If set to an empty string, this check is deactivated. This is also the default, meaning that
     * no checking is done.
     *
     * In case that site_casing() is set to a value other than `SiteCasing::kUnchanged`:
     * The validation is done after changing the casing, so that only lower or capital letters have
     * to be provided for validation. In case that site_casing() is set to `SiteCasing::kUnchanged`:
     * All chars that are to be considered valid have to be provided for validation.
     *
     * See `nucleic_acid_codes...()` and `amino_acid_codes...()` functions for presettings of chars
     * that can be used for validation here.
     */
    FastaReader& valid_chars( std::string const& chars );

    /**
     * @brief Return the currently set chars used for validating Sequence sites.
     *
     * An empty string means that no validation is done.
     */
    std::string  valid_chars() const;

    /**
     * @brief Return the internal CharLookup that is used for validating the Sequence sites.
     *
     * This function is provided in case direct access to the lookup is needed. Usually, the
     * valid_chars() function should suffice. See there for details.
     */
    utils::CharLookup<bool>& valid_char_lookup();

    // ---------------------------------------------------------------------
    //     Members
    // ---------------------------------------------------------------------

private:

    ParsingMethod           parsing_method_ = ParsingMethod::kDefault;

    SiteCasing              site_casing_      = SiteCasing::kToUpper;
    bool                    guess_abundances_ = false;
    bool                    use_validation_   = false;
    utils::CharLookup<bool> lookup_;

};

} // namespace sequence
} // namespace genesis

#endif // include guard
