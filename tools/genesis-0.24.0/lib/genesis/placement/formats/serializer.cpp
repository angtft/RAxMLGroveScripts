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
 * @brief Implementation of SampleSerializer class.
 *
 * @file
 * @ingroup placement
 */

#include "genesis/placement/formats/serializer.hpp"

#include "genesis/placement/formats/newick_reader.hpp"
#include "genesis/placement/formats/newick_writer.hpp"
#include "genesis/placement/sample.hpp"
#include "genesis/placement/sample_set.hpp"
#include "genesis/utils/io/deserializer.hpp"
#include "genesis/utils/io/serializer.hpp"

#include <stdexcept>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

namespace genesis {
namespace placement {

// =================================================================================================
//     Version
// =================================================================================================

/**
 * @brief Version of this serialization helper. Is written to the stream and read again to make
 * sure that different versions don't crash inexpectedly.
 */
unsigned char SampleSerializer::version = 1;

// =================================================================================================
//     Save
// =================================================================================================

/**
 * @brief Saves the Sample to a binary file that can later be read by using load().
 */
void SampleSerializer::save( Sample const& map, std::string const& file_name )
{
    // Prepare.
    utils::Serializer ser (file_name);
    if (!ser) {
        throw std::invalid_argument("Serialization failed.");
    }

    // Write header.
    char magic[] = "BPLACE\0\0";
    ser.put_raw(magic, 8);
    ser.put_int<unsigned char>(version);

    // Write tree.
    // TODO if there is a tree serialization in the future, this one could be used here, and in
    // addition to edge numbers, the edge indices can be stored, so that deserialization is easier.
    auto nw = PlacementTreeNewickWriter();
    nw.enable_names(true);
    nw.enable_branch_lengths(true);
    std::string tree;
    nw.write( map.tree(), utils::to_string( tree ));
    ser.put_string( std::move( tree ));

    // Write pqueries.
    ser.put_int(map.size());
    for (auto& pqry : map.pqueries()) {

        // Write placements.
        ser.put_int(pqry.placement_size());
        for( auto const& place : pqry.placements() ) {
            // We set the edge index instead of edge num. This is faster, simpler to resorte, and
            // consinstend with Pquery.add_placement() parameters.
            ser.put_int   (place.edge().index());

            ser.put_float( place.likelihood );
            ser.put_float( place.like_weight_ratio );
            ser.put_float( place.proximal_length );
            ser.put_float( place.pendant_length );
        }

        // Write names.
        ser.put_int(pqry.name_size());
        for( auto const& name : pqry.names() ) {
            ser.put_string (name.name);
            ser.put_float  (name.multiplicity);
        }
    }
}

// =================================================================================================
//     Load
// =================================================================================================

/**
 * @brief Loads a Sample from a binary file that was written by using save().
 */
Sample SampleSerializer::load( std::string const& file_name )
{
    // Create returned object.
    Sample map;

    // Prepare, check stream status.
    utils::Deserializer des( file_name );
    if( ! des ) {
        throw std::invalid_argument( "Deserialization failed: Cannot open file." );
    }

    // Read and check header.
    std::string magic = des.get_raw_string(8);
    if (strncmp (magic.c_str(), "BPLACE\0\0", 8) != 0) {
        throw std::invalid_argument("Wrong file format: \"" + magic + "\".");
    }
    auto ver = des.get_int<unsigned char>();
    if (ver != version) {
        throw std::invalid_argument("Wrong serialization version: " + std::to_string(ver));
    }

    // Read and check tree.
    auto tree_string = des.get_string();
    map.tree() = PlacementTreeNewickReader().read( utils::from_string( tree_string ));

    // Read pqueries.
    size_t num_pqueries = des.get_int<size_t>();
    for (size_t i = 0; i < num_pqueries; ++i) {
        Pquery& pqry = map.add();

        // Read placements.
        size_t num_place = des.get_int<size_t>();
        for (size_t p = 0; p < num_place; ++p) {
            // Get edge index, add the placement there.
            size_t edge_idx = des.get_int<size_t>();
            auto&   edge    = map.tree().edge_at( edge_idx );
            auto&  place    = pqry.add_placement( edge );

            place.likelihood        = des.get_float<double>();
            place.like_weight_ratio = des.get_float<double>();
            place.proximal_length   = des.get_float<double>();
            place.pendant_length    = des.get_float<double>();
        }

        // Read names.
        size_t num_names = des.get_int<size_t>();
        for (size_t n = 0; n < num_names; ++n) {
            auto name = pqry.add_name( des.get_string() );
            name.multiplicity = des.get_float<double>();
        }
    }

    if (!des.finished()) {
        throw std::invalid_argument("Deserialization failed: File longer than expected.");
    }

    return map;
}

SampleSet SampleSerializer::load( std::vector<std::string> const& file_names )
{
    SampleSet sample_set;
    load( file_names, sample_set );
    return sample_set;
}

void SampleSerializer::load( std::vector<std::string> const& file_names, SampleSet& sample_set )
{
    #if defined( GENESIS_OPENMP )

        // Make a vector of default-constructed Samples of the needed size.
        // We do this so that the order of input jplace files is kept.
        auto tmp = std::vector<Sample>( file_names.size() );

        // Parallel loading.
        #pragma omp parallel for
        for( size_t i = 0; i < file_names.size(); ++i ) {
            tmp[ i ] = load( file_names[i] );
        }

        // Move to target SampleSet.
        for( size_t i = 0; i < file_names.size(); ++i ) {
            auto const name = utils::file_filename( utils::file_basename( file_names[i] ) );
            sample_set.add( std::move( tmp[i] ), name );
        }

    #else

        for( auto const& fn : file_names ) {
            auto const name = utils::file_filename( utils::file_basename(fn) );
            sample_set.add( load( fn ), name );
        }

    #endif
}

} // namespace placement
} // namespace genesis
