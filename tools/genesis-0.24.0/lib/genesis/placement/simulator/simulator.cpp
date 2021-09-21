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
 * @brief Implementation of Placement Simulator class.
 *
 * @file
 * @ingroup placement
 */

#include "genesis/placement/simulator/simulator.hpp"

#include "genesis/placement/function/functions.hpp"

#include <cassert>

namespace genesis {
namespace placement {

// =================================================================================================
//     Placement Simulator
// =================================================================================================

/**
 * @brief Generate `n` many Pqueries and place them in the Sample.
 */
void Simulator::generate( Sample& sample, size_t n )
{
    // Prepare distributions.
    edge_distribution_.prepare( sample );
    extra_placement_distribution_.prepare( sample );
    like_weight_ratio_distribution_.prepare( sample );
    proximal_length_distribution_.prepare( sample );
    pendant_length_distribution_.prepare( sample );

    // Helper function to add a Placement with random properties.
    auto add_random_placement = [&] ( Pquery& pqry, PlacementTreeEdge& edge ) {
        PqueryPlacement& place  = pqry.add_placement( edge );
        place.proximal_length   = proximal_length_distribution_.generate( edge );
        place.pendant_length    = pendant_length_distribution_.generate( edge );
        place.like_weight_ratio = like_weight_ratio_distribution_.generate();
    };

    for (size_t i = 0; i < n; ++i) {
        // Generate one Pquery.
        Pquery& pqry = sample.add();
        pqry.add_name( "pquery_" + std::to_string( sample.size() - 1 ));

        // Get a random edge.
        size_t edge_idx = edge_distribution_.generate();
        assert( edge_idx < sample.tree().edge_count() );
        auto& edge = sample.tree().edge_at( edge_idx );

        // Add a placement at the edge.
        add_random_placement( pqry, edge );

        // Generate additional placements around that edge.
        auto extra_edge_idcs = extra_placement_distribution_.generate( edge );
        for( auto extra_edge_idx : extra_edge_idcs ) {

            // Make sure that the distribution works correctly and actually gives us only indices
            // that are not our current central one.
            assert( extra_edge_idx != edge_idx );
            assert( extra_edge_idx < sample.tree().edge_count() );

            // Add a placement at the edge.
            auto& extra_edge = sample.tree().edge_at( extra_edge_idx );
            add_random_placement( pqry, extra_edge );
        }

        // Make sure that our generating process ends up with correct like weights ratios.
        // We generated them randomly, so their sum will not be 1.0 in general.
        normalize_weight_ratios( pqry );
    }
}

} // namespace placement
} // namespace genesis
