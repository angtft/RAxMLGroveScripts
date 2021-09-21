/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2018 Lucas Czech and HITS gGmbH

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
 * @ingroup placement
 */

#include "genesis/placement/function/epca.hpp"

#include "genesis/placement/function/functions.hpp"
#include "genesis/placement/function/helper.hpp"
#include "genesis/placement/function/masses.hpp"
#include "genesis/placement/function/sample_set.hpp"

#include "genesis/tree/function/functions.hpp"
#include "genesis/tree/iterator/postorder.hpp"

#include "genesis/utils/containers/matrix.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/math/common.hpp"
#include "genesis/utils/math/matrix.hpp"
#include "genesis/utils/math/pca.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

namespace genesis {
namespace placement {

// =================================================================================================
//     Edge PCA Imbalance Vector
// =================================================================================================

std::vector<double> epca_imbalance_vector( Sample const& sample, bool normalize )
{
    // Result vector: imbalance of masses at each edge of the tree.
    auto vec = std::vector<double>( sample.tree().edge_count(), 0.0 );

    // We need the masses per edge, and their sum, for later.
    auto const masses   = placement_mass_per_edges_with_multiplicities( sample );
    auto const mass_sum = std::accumulate( masses.begin(), masses.end(), 0.0 );

    // Collect the placement masses at each link of the tree.
    // Use init to -1 as indicator for assertions.
    auto link_masses = std::vector<double>( sample.tree().link_count(), -1.0 );

    for( auto tree_it : postorder( sample.tree() )) {

        // Skip the last iteration. We are interested in edges, not in nodes.
        if( tree_it.is_last_iteration() ) {
            continue;
        }

        // Get the indices of the links at both sides of the current edge.
        // cur_idx is the link that points away from the root,
        // out_idx is the link that points towards it (i.e., its subtree contains the root).
        auto const cur_idx = tree_it.link().index();
        auto const out_idx = tree_it.link().outer().index();

        // Assert that we have not processed those links before.
        assert( link_masses[ cur_idx ] < 0.0 );
        assert( link_masses[ out_idx ] < 0.0 );

        // Assert that the cur_idx belongs to the link away from the root.
        // This is the case if the primary link of its node is the link itself,
        // because the node uses this link to point towards the root - thus, the link itself
        // is away from the root, while the out_idx link lies towards it.
        assert( sample.tree().link_at( cur_idx ).node().primary_link().index() == cur_idx );

        // Some more ways to do the same assertion, just to be sure.
        assert( tree_it.edge().index() == sample.tree().link_at( cur_idx ).edge().index() );
        assert( tree_it.edge().primary_link().index() == out_idx );
        assert( tree_it.edge().secondary_link().index() == cur_idx );

        // Leaf links have no mass.
        if( is_leaf( tree_it.link() )) {
            link_masses[ cur_idx ] = 0.0;

        // If the link belongs to an inner node, we calculate its mass as the sum of the masses
        // of the other links of this node. Those have already been processed, as we are doing
        // postorder traversal.
        } else {

            // Collect the mass.
            double round_sum = 0.0;

            // Iterate around all other links of the node that belongs to the cur_idx link.
            auto round_link = &tree_it.link().next();
            while( round_link != &tree_it.link() ) {

                // We are doing postorder traversal, so we should have seen this link before.
                assert( link_masses[ round_link->index() ] >= 0.0 );

                // The mass of the subtree behind this link can be calculated from the total mass
                // minus the mass of the link itself.
                round_sum += mass_sum - link_masses[ round_link->index() ];

                // Next link of the node.
                round_link = &round_link->next();
            }

            // The sum should always be >0, but for numerical reaons, we better make sure it is.
            link_masses[ cur_idx ] = std::max( 0.0, round_sum );
        }

        // Calculate the mass at the other side of the edge. We need to correct negative values,
        // which can occur for numerical reasons (in the order of e-12).
        link_masses[out_idx] = std::max(
            0.0,
            mass_sum - link_masses[ cur_idx ] - masses[ tree_it.edge().index() ]
        );

        // Make sure we have processed all masses that we are going to use.
        assert( link_masses[cur_idx] > -1.0 );
        assert( link_masses[out_idx] > -1.0 );

        // Finally, calculate the imbalance of the current edge,
        // normalized by the total mass on the tree (expect for the mass of the current edge).
        auto const imbalance = link_masses[cur_idx] - link_masses[out_idx];
        if( normalize ) {
            auto const normalizer = mass_sum - masses[ tree_it.edge().index() ];
            assert( normalizer > 0.0 );
            vec[ tree_it.edge().index() ] = imbalance / normalizer;
        } else {
            vec[ tree_it.edge().index() ] = imbalance;
        }
    }

    return vec;
}

// =================================================================================================
//     Edge PCA Imbalance Matrix
// =================================================================================================

utils::Matrix<double> epca_imbalance_matrix(
    SampleSet const& samples,
    bool include_leaves,
    bool normalize
) {
    // If there are no samples, return empty matrix.
    if( samples.size() == 0 ) {
        return utils::Matrix<double>();
    }

    // Check if all trees have the same tpology and edge nums.
    if( ! all_identical_trees( samples )) {
        throw std::runtime_error(
            "Cannot calculate Edge PCA on trees that have a different topology."
        );
    }

    assert( samples.size() > 0 );
    auto const edge_count = samples.at( 0 ).tree().edge_count();

    if( include_leaves ) {

        auto imbalance_matrix = utils::Matrix<double>( samples.size(), edge_count );

        #pragma omp parallel for
        for( size_t s = 0; s < samples.size(); ++s ) {
            auto const& smp = samples[s];
            auto const imbalance_vec = epca_imbalance_vector( smp, normalize );

            // We need to have the right number of imbalance values.
            assert( imbalance_vec.size() == edge_count );

            // Copy imbalance values to the matrix.
            for( size_t i = 0; i < edge_count; ++i ) {

                // Either the edge is an inner edge, or (if not, i.e., it leads to a leaf),
                // it's imbalance is minus 1, as all its mass is on the root side.
                assert(
                    ! is_leaf( samples[s].tree().edge_at(i).secondary_node() ) ||
                    utils::almost_equal_relative( imbalance_vec[ i ], -1.0 )
                );

                imbalance_matrix( s, i ) = imbalance_vec[ i ];
            }
        }

        return imbalance_matrix;

    } else {

        // Get the indices of all edges that do not lead to a tip.
        auto const inner_edge_indices = tree::inner_edge_indices( samples.at( 0 ).tree() );

        // Prepare result
        auto imbalance_matrix = utils::Matrix<double>( samples.size(), inner_edge_indices.size() );

        #pragma omp parallel for
        for( size_t s = 0; s < samples.size(); ++s ) {
            auto const& smp = samples[s];
            auto const imbalance_vec = epca_imbalance_vector( smp, normalize );

            // We need to have the right number of imbalance values, which also needs to be
            // smaller than the number of inner edges (there can be no tree with just inner
            // edges, thus the total number of edges has to be bigger).
            assert( imbalance_vec.size() == edge_count );
            assert( imbalance_vec.size() >  inner_edge_indices.size() );

            // Copy those imbalance values to the matrix that belong to inner edges.
            for( size_t i = 0; i < inner_edge_indices.size(); ++i ) {
                auto idx = inner_edge_indices[i];
                imbalance_matrix( s, i ) = imbalance_vec[ idx ];
            }
        }

        return imbalance_matrix;
    }
}

// =================================================================================================
//     Splitify Transform with Kappa
// =================================================================================================

void epca_splitify_transform( utils::Matrix<double>& imbalance_matrix, double kappa )
{
    // Precondition check.
    if( kappa < 0.0 ) {
        throw std::runtime_error( "Argument for kappa must be non-negative." );
    }

    // Save time if the transformation throws away the actual value.
    // We do not need to calculate the abs and power in this case.
    if( kappa == 0.0 ) {
        for( auto& elem : imbalance_matrix ) {
            elem = static_cast<double>( utils::signum( elem ));
        }
        return;
    }

    // Save time if the transformation does not change anything.
    if( kappa == 1.0 ) {
        return;
    }

    // If neither applies, do the full transformation.
    for( auto& elem : imbalance_matrix ) {
        elem = static_cast<double>( utils::signum( elem )) * std::pow( std::abs( elem ), kappa );
    }
}

// =================================================================================================
//     Edge PCA
// =================================================================================================

EpcaData epca( SampleSet const& samples, double kappa, double epsilon, size_t components )
{
    // If there are no samples, return empty result.
    if( samples.size() == 0 ) {
        return EpcaData();
    }

    // Calculate the imbalance_matrix.
    auto imbalance_matrix = epca_imbalance_matrix( samples, false );
    assert( samples.size() > 0 );
    assert( imbalance_matrix.rows() == samples.size() );
    assert( imbalance_matrix.cols() == tree::inner_edge_count( samples[0].tree() ) );

    // Get the indices of the inner edges.
    auto const inner_edge_indices = tree::inner_edge_indices( samples.at( 0 ).tree() );
    assert( imbalance_matrix.cols() == inner_edge_indices.size() );

    // Filter and transform the imbalance matrix.
    auto const not_filtered_indices = filter_constant_columns( imbalance_matrix, epsilon );
    epca_splitify_transform( imbalance_matrix, kappa );

    // We now use the list of not filtered indices to selected from the list of inner edge indices.
    // The result is just the indices of the edges that are still in the matrix.
    std::vector<size_t> edge_indices;
    for( auto const& not_filt : not_filtered_indices ) {
        edge_indices.push_back( inner_edge_indices[ not_filt ] );
    }
    assert( edge_indices.size() == imbalance_matrix.cols() );

    // Get correct number of pca components.
    if( components == 0 || components > imbalance_matrix.cols() ) {
        components = imbalance_matrix.cols();
    }

    // Run and return PCA.
    auto pca = utils::principal_component_analysis(
        imbalance_matrix, components, utils::PcaStandardization::kCovariance
    );
    assert( pca.eigenvalues.size()  == components );
    assert( pca.eigenvectors.rows() == edge_indices.size() );
    assert( pca.eigenvectors.cols() == components );
    assert( pca.projection.rows()   == samples.size() );
    assert( pca.projection.cols()   == components );

    // Move data.
    EpcaData result;
    result.eigenvalues  = std::move( pca.eigenvalues );
    result.eigenvectors = std::move( pca.eigenvectors );
    result.projection   = std::move( pca.projection );
    result.edge_indices = std::move( edge_indices );
    return result;
}

} // namespace placement
} // namespace genesis
