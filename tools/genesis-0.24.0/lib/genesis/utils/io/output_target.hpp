#ifndef GENESIS_UTILS_IO_OUTPUT_TARGET_H_
#define GENESIS_UTILS_IO_OUTPUT_TARGET_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2020 Lucas Czech

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

#include "genesis/utils/io/base_output_target.hpp"
#include "genesis/utils/io/file_output_target.hpp"
#include "genesis/utils/io/gzip_output_target.hpp"
#include "genesis/utils/io/gzip_stream.hpp"
#include "genesis/utils/io/stream_output_target.hpp"
#include "genesis/utils/io/string_output_target.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/std.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace genesis {
namespace utils {

// =================================================================================================
//     Output Target Convenience Functions
// =================================================================================================

/**
 * @brief Obtain an output target for writing to a file.
 *
 * The output target returned from this function can be used in the writer classes, e.g.,
 * placement::JplaceWriter or sequence::FastaWriter.
 *
 * If @p compression_level is set to a compression level other than ::GzipCompressionLevel::kNoCompression
 * (which is the default, which means, no compression by default), the output is compressed using
 * gzip. We recommend to use ::GzipCompressionLevel::kDefaultCompression%.
 *
 * Furthermore, if @p auto_adjust_filename is set to `true` (default), the file name is
 * adjusted according to the compression setting: If compression is used, the file name
 * is appended by the `.gz` extension, if this is not alreay present. (For completeness,
 * the oppositve also works: If the file name ends in `.gz`, but no compression is chosen, the `.gz`
 * extension is removed.)
 *
 * If the file cannot be written to, the function throws an exception.
 * Also, by default, if the file already exists, an exception is thrown. See
 * @link utils::Options::allow_file_overwriting( bool ) Options::allow_file_overwriting()@endlink
 * to change this behaviour.
 */
inline std::shared_ptr<BaseOutputTarget> to_file(
    std::string const& file_name,
    GzipCompressionLevel compression_level = GzipCompressionLevel::kNoCompression,
    bool auto_adjust_filename = true
) {
    auto fn = file_name;
    auto const ext = file_extension( fn );

    // With compression
    if( compression_level != GzipCompressionLevel::kNoCompression ) {
        if( auto_adjust_filename && ext != "gz" ) {
            fn += ".gz";
        }
        return std::make_shared<GzipOutputTarget>(
            std::make_shared< FileOutputTarget >( fn, std::ios_base::out | std::ios_base::binary ),
            compression_level
        );
    }

    // Without compression. Not in the `else` branch, to not confuse old compilers.
    if( auto_adjust_filename && ext == "gz" ) {
        fn.erase( fn.size() - 2 );
    }
    return std::make_shared< FileOutputTarget >( fn );
}

/**
 * @brief Obtain an output target for writing to a stream.
 *
 * The output target returned from this function can be used in the writer classes, e.g.,
 * placement::JplaceWriter or sequence::FastaWriter.
 *
 * If @p compression_level is set to a compression level other than ::GzipCompressionLevel::kNoCompression
 * (which is the default, which means, no compression by default), the output is compressed using
 * gzip. In that case, it is recommended that the @p target_stream uses `std::ios_base::binary`
 * when opening the stream.
 */
inline std::shared_ptr<BaseOutputTarget> to_stream(
    std::ostream& target_stream,
    GzipCompressionLevel compression_level = GzipCompressionLevel::kNoCompression
) {
    if( compression_level != GzipCompressionLevel::kNoCompression ) {
        return std::make_shared<GzipOutputTarget>(
            std::make_shared< StreamOutputTarget >( target_stream ),
            compression_level
        );
    }
    return std::make_shared< StreamOutputTarget >( target_stream );
}

/**
 * @brief Obtain an output target for writing to a string.
 *
 * The output target returned from this function can be used in the writer classes, e.g.,
 * placement::JplaceWriter or sequence::FastaWriter.
 */
inline std::shared_ptr<BaseOutputTarget> to_string(
    std::string& target_string
) {
    return std::make_shared< StringOutputTarget >( target_string );
}

} // namespace utils
} // namespace genesis

#endif // include guard
