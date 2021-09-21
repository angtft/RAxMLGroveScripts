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
 * @ingroup tree
 */

#include "genesis/tree/function/functions.hpp"

#include "genesis/tree/function/distances.hpp"
#include "genesis/tree/function/operators.hpp"
#include "genesis/tree/iterator/eulertour.hpp"
#include "genesis/tree/iterator/preorder.hpp"
#include "genesis/tree/tree.hpp"
#include "genesis/tree/tree/subtree.hpp"

#include "genesis/utils/containers/matrix/operators.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <unordered_set>
#include <vector>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

namespace genesis {
namespace tree {

// =================================================================================================
//     Node Properties
// =================================================================================================

bool is_leaf( TreeLink const& link )
{
    return &( link.next() ) == &link;
}

bool is_leaf( TreeNode const& node )
{
    return is_leaf( node.link() );
}

bool is_leaf( TreeEdge const& edge )
{
    return is_leaf( edge.secondary_link() );
}

bool is_inner( TreeLink const& link )
{
    return &( link.next() ) != &link;
}

bool is_inner( TreeNode const& node )
{
    return is_inner( node.link() );
}

bool is_inner( TreeEdge const& edge )
{
    return is_inner( edge.secondary_link() );
}

bool is_root( TreeLink const& link )
{
    return is_root( link.node() );
}

bool is_root( TreeNode const& node )
{
    // The link_ is always the one pointing towards the root. Also, the edge of that link always has
    // the primary link set to that it points towards the root.
    // At the root itself, however, this means we are pointing to ourselves. Use this to check
    // if this node is the root.
    return &( node.link().edge().primary_link() ) == &( node.link() );
}

size_t degree( TreeLink const& link )
{
    return degree( link.node() );
}

size_t degree( TreeNode const& node )
{
    size_t dgr = 0;
    TreeLink const* lnk = &node.link();

    do {
        ++dgr;
        lnk = &lnk->next();
    } while( lnk != &node.link() );

    return dgr;
}

// =================================================================================================
//     Node Count Properties
// =================================================================================================

size_t max_degree( Tree const& tree )
{
    size_t max = 0;
    for( size_t i = 0; i < tree.node_count(); ++i ) {
        max = std::max( max, degree( tree.node_at(i) ) );
    }
    return max;
}

bool is_bifurcating( Tree const& tree, bool loose )
{
    // Iterate all nodes and verify their degree.
    for( size_t i = 0; i < tree.node_count(); ++i ) {
        auto const deg = degree( tree.node_at(i) );

        // Any degree > 3 is multifurcating.
        if( deg > 3 ) {
            return false;
        }

        // A degree of 2 is okay of we are loose, and always okay for the root.
        if( deg == 2 ) {
            auto const isroot = ( tree.node_at(i).index() == tree.root_node().index() );
            assert( isroot == is_root( tree.node_at(i) ));

            if( ! isroot && ! loose ) {
                return false;
            }
        }
    }
    return true;
}

bool is_binary( Tree const& tree, bool loose )
{
    return is_bifurcating( tree, loose );
}

bool is_rooted( Tree const& tree )
{
    return degree( tree.root_node() ) == 2;
}

size_t leaf_node_count( Tree const& tree )
{
    size_t sum = 0;
    for (size_t i = 0; i < tree.node_count(); ++i) {
        auto const& n = tree.node_at(i);
        if( is_leaf(n) ) {
            ++sum;
        }
    }
    return sum;
}

size_t inner_node_count( Tree const& tree )
{
    return tree.node_count() - leaf_node_count( tree );
}

size_t node_count( Tree const& tree )
{
    return tree.node_count();
}

size_t leaf_edge_count(  Tree const& tree )
{
    size_t sum = 0;
    for( auto const& edge : tree.edges() ) {
        if( is_leaf( edge.primary_node() ) || is_leaf( edge.secondary_node() ) ) {
            ++sum;
        }
    }
    return sum;
}

size_t inner_edge_count( Tree const& tree )
{
    size_t sum = 0;
    for( auto const& edge : tree.edges() ) {
        if( is_inner( edge.primary_node() ) && is_inner( edge.secondary_node() ) ) {
            ++sum;
        }
    }
    return sum;
}

size_t edge_count( Tree const& tree )
{
    return tree.edge_count();
}

std::vector<size_t> inner_edge_indices( Tree const& tree )
{
    std::vector<size_t> result;
    for( auto const& edge : tree.edges() ) {
        if( is_inner( edge.secondary_node() ) ) {
            result.push_back( edge.index() );
        }
    }
    return result;
}

std::vector<size_t> leaf_edge_indices( Tree const& tree )
{
    std::vector<size_t> result;
    for( auto const& edge : tree.edges() ) {
        if( is_leaf( edge.secondary_node() ) ) {
            result.push_back( edge.index() );
        }
    }
    return result;
}

std::vector<size_t> inner_node_indices( Tree const& tree )
{
    std::vector<size_t> result;
    for( auto const& node : tree.nodes() ) {
        if( is_inner( node )) {
            result.push_back( node.index() );
        }
    }
    return result;
}

std::vector<size_t> leaf_node_indices( Tree const& tree )
{
    std::vector<size_t> result;
    for( auto const& node : tree.nodes() ) {
        if( is_leaf( node )) {
            result.push_back( node.index() );
        }
    }
    return result;
}

// =================================================================================================
//     Tree Sides
// =================================================================================================

utils::Matrix<signed char> edge_sides( Tree const& tree )
{
    // Get a quadratic matrix: for each edge, it gives a value whether each other edge is
    // proximal (1) or distal (-1) relative to itself (0).
    auto result = utils::Matrix<signed char>( tree.edge_count(), tree.edge_count(), 0 );

    // Helper function that traverses the subtree starting at a link,
    // and for each edge in the subtree, sets its entry in the matrix to the given sign.
    auto traverse = [&]( TreeLink const& start_link, size_t i, signed char sign ){
        auto link = &( start_link.next() );
        while( link != &start_link ) {
            result( i, link->edge().index() ) = sign;
            link = &( link->outer().next() );
        }
    };

    // For each edge, do the traversal in both directions and set the signs.
    // This can probably be done more efficiently with only one smart traversal of the whole tree,
    // but this function is not needed often enough right now.
    for( size_t i = 0; i < tree.edge_count(); ++i ) {
        traverse( tree.edge_at(i).primary_link(), i,   -1 );
        traverse( tree.edge_at(i).secondary_link(), i,  1 );
    }

    return result;
}

utils::Matrix<signed char> node_root_direction_matrix( Tree const& tree )
{
    auto mat = utils::Matrix<signed char>( tree.node_count(), tree.node_count(), 0 );

    // Fill every row of the matrix.
    #pragma omp parallel for
    for( size_t i = 0; i < tree.node_count(); ++i ) {
        auto const& row_node     = tree.node_at(i);
        auto const  row_index    = row_node.index();
        auto const  primary_link = &row_node.primary_link();

        // Fill root side subtree with 1s.
        // We do set the value of inner nodes multiple times, but that's no problem.
        // Also, we need to do an extra check for the root here, in order to set all
        // subtrees of the root to -1.
        signed char const value = is_root( row_node ) ? -1 : 1;
        auto current_link = &( primary_link->outer() );
        while( current_link != primary_link ) {
            mat( row_index, current_link->node().index() ) = value;
            current_link = &( current_link->next().outer() );
        }

        // Fill all non-root side subtrees with -1s.
        // We explicitly go through all non-root links of the node, in order to be really clear
        // about our intentions. It would also work to simply follow the link chain until
        // we reach the primary link again. However, this would also set -1 to our row node,
        // thus we'd have to reset it, making the algorithm a bit messy.
        // So, to be clear and clean, we avoid this.
        auto sub_link = &( primary_link->next() );
        while( sub_link != primary_link ) {

            // Now, for a given non-root subtree, set everything to -1.
            current_link = &( sub_link->outer() );
            while( current_link != sub_link ) {
                mat( row_index, current_link->node().index() ) = -1;
                current_link = &( current_link->next().outer() );
            }

            // Go to next subtree.
            sub_link = &( sub_link->next() );
        }

        // Check that the diagonal element is untouched.
        assert( mat( row_index, row_index ) == 0 );
    }

    return mat;
}

// =================================================================================================
//     Subtrees
// =================================================================================================

size_t subtree_size( Tree const& tree, TreeLink const& link )
{
    if( ! belongs_to( tree, link )) {
        throw std::runtime_error(
            "Cannot caluclate subtree_size, as the given Link does not belong to the Tree."
        );
    }

    // TODO This is a quick and dirty solution. Traverse the whole subtree, add all nodes to a set
    // and simply return the size of that set.

    std::unordered_set< TreeNode const* > visited_nodes;

    auto cur_link = &link.outer();
    while( cur_link != &link ) {
        visited_nodes.insert( &cur_link->node() );
        cur_link = &cur_link->next().outer();
    }

    return visited_nodes.size();
}

std::vector<size_t> subtree_sizes( Tree const& tree, TreeNode const& node )
{
    if( ! belongs_to( tree, node )) {
        throw std::runtime_error(
            "Cannot calculate subtree_sizes(), as the given Node does not belong to the Tree."
        );
    }

    // TODO this is an overly complex and thus error prone solution, maybe there is a better way?!

    // Prepare result vector.
    std::vector<size_t> result;
    result.resize( tree.node_count(), 0 );

    // We use a stack to track the subtree sizes.
    // We store the entry link of the preorder traversal of the nodes. The entry link is the one
    // that is given when visiting the node first while doing a eulertour traversal of the tree.
    // This is always the next() link after the towards-the-starting-node/root link.
    std::vector< TreeLink const* > stack;
    stack.push_back( &node.link() );

    // Traverse the tree.
    for( auto it : eulertour( node )) {

        // If this is the last time we visit that node on our way back up the tree.
        // (The second part of the condition checks whether it is the starting node, because
        // in this case, we do not want to remove it.)
        if( &it.link().next() == stack.back() && stack.back() != &node.link() ) {

            // We finished with a subtree. Add the cummulative number of children of that subtree
            // to the parent node, and remove the parent from the stack (as we are done with it).
            auto st_size = result[ stack.back()->node().index() ];
            stack.pop_back();
            result[ stack.back()->node().index() ] += st_size;

        // If this node is already the current top stack element.
        } else if( &it.node() == &stack.back()->node() ) {

            // Do nothing.

        // If it is a leaf.
        } else if( is_leaf( it.link() )) {

            // Simply increment its parent's counter.
            ++result[ stack.back()->node().index() ];

        // If we will visit that node in the future again.
        } else {

            // Add a count for the immediate child (i.e., the current node) to the current stack
            // end (i.e., increment the counter of children of that node),
            // then add the current node itself to the stack, so that in the next iteration,
            // we will increase its counts.
            ++result[ stack.back()->node().index() ];
            stack.push_back( &it.link() );
        }
    }

    // The stack now should contain only a single node, which is the starting node itself.
    assert( stack.size() == 1 && stack.back() == &node.link() );

    // The size of the subtree of the starting node is always the number of nodes in the tree
    // minus one for that node itself (as it is not counted as part of its subtree).
    assert( result[ node.index() ] == tree.node_count() - 1 );

    return result;
}

std::vector<size_t> subtree_sizes( Tree const& tree )
{
    return subtree_sizes( tree, tree.root_node() );
}

size_t subtree_max_path_height( Tree const& tree, TreeLink const& link )
{
    if( ! belongs_to( tree, link )) {
        throw std::runtime_error(
            "Cannot calculate subtree_max_path_height(), "
            "as the given Link does not belong to the Tree."
        );
    }

    // TODO make more efficient. no need for full dist vector.
    auto dists = node_path_length_vector( tree, link.outer().node() );
    size_t max = 0;

    auto cur_link = &link.outer();
    while( cur_link != &link ) {
        max = std::max( max, dists[ cur_link->node().index() ] );
        cur_link = &cur_link->next().outer();
    }
    return max;
}

std::vector<size_t> subtree_max_path_heights( Tree const& tree, TreeNode const& node )
{
    if( ! belongs_to( tree, node )) {
        throw std::runtime_error(
            "Cannot calculate subtree_max_path_heights(), "
            "as the given Node does not belong to the Tree."
        );
    }

    auto result = std::vector<size_t>( tree.node_count(), 0 );

    // Recursive helper function that evaluates the wanted size for a given subtree,
    // stores the result in the vector and returns it for recursive usage.
    std::function< size_t( TreeLink const* )> rec_subtree_height = [&]( TreeLink const* l ) {
        size_t link_max = 0;
        TreeLink const* cl = &l->next();
        while( cl != l ) {
            link_max = std::max( link_max, 1 + rec_subtree_height( &cl->outer() ));
            cl = &cl->next();
        }

        result[ l->node().index() ] = link_max;
        return link_max;
    };

    // Loop all subtrees of the given node and find the highest.
    // This loop is a bit different from the one in the recursive function, as we need to evaluate
    // all links of the given starting node, instead of just the ones away from the start node.
    size_t node_max = 0;
    TreeLink const* cur_l = &node.link();
    do {
        node_max = std::max( node_max, 1 + rec_subtree_height( &cur_l->outer() ));
        cur_l = &cur_l->next();
    } while( cur_l != &node.link() );
    result[ node.index() ] = node_max;

    return result;
}

std::vector<size_t> subtree_max_path_heights( Tree const& tree )
{
    return subtree_max_path_heights( tree, tree.root_node() );
}

utils::Matrix<signed char> sign_matrix( Tree const& tree, bool compressed )
{
    // Edge cases and input checks.
    if( tree.empty() ) {
        return utils::Matrix<signed char>();
    }
    if( ! is_rooted( tree )) {
        throw std::invalid_argument( "Tree is not rooted. Cannot calculate its sign matrix." );
    }
    if( ! is_bifurcating( tree )) {
        throw std::invalid_argument( "Tree is not bifurcating. Cannot calculate its sign matrix." );
    }

    // Prepare a result matrix of the full size. For the compressed version,
    // we later replate it again.
    auto result = utils::Matrix<signed char>( tree.node_count(), tree.node_count(), 0 );

    // Helper function that fills all columns of a subtree with a given sign.
    auto fill_subtree_indices = [&]( size_t row_idx, Subtree const& st, signed char sign ){
        for( auto const& it : preorder(st) ) {
            result( row_idx, it.node().index() ) = sign;
        }
    };

    // Fill every row of the matrix.
    #pragma omp parallel for
    for( size_t i = 0; i < tree.node_count(); ++i ) {
        auto const& row_node = tree.node_at(i);
        auto const  row_idx  = row_node.index();

        if( row_idx == tree.root_node().index() ) {

            // The root node is special: we use its two subtrees directly.
            assert( &row_node.link().next().next() == &row_node.link() );
            fill_subtree_indices( row_idx, Subtree{ row_node.link().outer() },        +1 );
            fill_subtree_indices( row_idx, Subtree{ row_node.link().next().outer() }, -1 );

        } else if( is_inner( row_node )) {

            // All other inner nodes are filled using their subtrees.
            assert( &row_node.link().next().next().next() == &row_node.link() );
            fill_subtree_indices( row_idx, Subtree{ row_node.link().next().outer() },        +1 );
            fill_subtree_indices( row_idx, Subtree{ row_node.link().next().next().outer() }, -1 );
        }
    }

    // For the compressed version, we re-use the previous result matrix,
    // and simply fill a new one with the needed rows and columns.
    // The data is not too big, and this is way easier and cleaner to implement.
    if( compressed ) {
        // Create a matrix with rows for each inner node and columns for each tip node.
        auto const in_node_idcs = inner_node_indices( tree );
        auto const lf_node_idcs = leaf_node_indices( tree );
        auto result_cmpr = utils::Matrix<signed char>( in_node_idcs.size(), lf_node_idcs.size(), 0 );

        // Fill the matrix at the indices that belong to inner nodes (for rows) and
        // leaf nodes (for columns).
        for( size_t r = 0; r < in_node_idcs.size(); ++r ) {
            for( size_t c = 0; c < lf_node_idcs.size(); ++c ) {
                result_cmpr( r, c ) = result( in_node_idcs[r], lf_node_idcs[c] );
            }
        }

        // Replace the result matrix efficiently.
        using std::swap;
        swap( result, result_cmpr );
    }

    return result;
}

// =================================================================================================
//     Misc
// =================================================================================================

std::vector< TreeLink const* > path_to_root( TreeNode const& node )
{
    std::vector< TreeLink const* > path;

    // Move towards the root and record all links in between.
    TreeLink const* cur_link = &node.primary_link();
    while( &cur_link->edge().secondary_link() == cur_link ) {

        // The above while condition means: is it the root?! Assert, that the default way of
        // checking for the root by using the node gives the same result.
        assert( ! is_root( cur_link->node() ));

        // Assert that the primary direction is correct.
        assert( cur_link == &cur_link->edge().secondary_link() );

        // Add the primary link of the current node to the list.
        path.push_back( cur_link );

        // Move one node towards the root.
        // Assert that the default way of finding the next node towards the root (by using
        // the edge) gives the same result as simply using the link's outer node.
        // This is the case because the cur link is the one that points towards the root
        // (which was asserted above).
        assert( &cur_link->edge().primary_link() == &cur_link->outer() );
        cur_link = &cur_link->outer().node().primary_link();
    }

    // Now finally add the root itself and return the list.
    assert( is_root( cur_link->node() ));
    path.push_back( cur_link );
    return path;
}

TreeNode const& lowest_common_ancestor( TreeNode const& node_a, TreeNode const& node_b )
{
    // Speedup and simplification.
    if( &node_a == &node_b ) {
        return node_a;
    }

    auto path_a = path_to_root( node_a );
    auto path_b = path_to_root( node_b );

    // We must have at least the two original links in the front and the root in the back.
    assert( path_a.size() > 0 && path_b.size() > 0 );
    assert( path_a.front()    == &node_a.link() );
    assert( path_b.front()    == &node_b.link() );
    assert( path_a.back()     == path_b.back() );

    // Remove from back as long as the last two elements are the same.
    // At the end of this, the remaining links are the ones on the path between
    // the two original links.
    while(
        path_a.size() > 1 &&
        path_b.size() > 1 &&
        path_a.at( path_a.size() - 1 ) == path_b.at( path_b.size() - 1 ) &&
        path_a.at( path_a.size() - 2 ) == path_b.at( path_b.size() - 2 )
    ) {
        path_a.pop_back();
        path_b.pop_back();
    }

    // Now, the last elements need to be the same (the LCA of the start and finish node).
    assert( path_a.size() > 0 && path_b.size() > 0 );
    assert( path_a.back()     == path_b.back() );

    return path_a.back()->node();
}

TreeNode&       lowest_common_ancestor( TreeNode& node_a,       TreeNode& node_b )
{
    auto const& c_node_a = static_cast< TreeNode const& >( node_a );
    auto const& c_node_b = static_cast< TreeNode const& >( node_b );
    return const_cast< TreeNode& >( lowest_common_ancestor( c_node_a, c_node_b ));
}

utils::Matrix<size_t> lowest_common_ancestors( Tree const& tree )
{
    auto res = utils::Matrix<size_t>( tree.node_count(), tree.node_count() );

    // This is not the best way to calculate all pairwise LCAs.
    // In the Quartet Scores code, we use range minimum queries and eulertours to achive the
    // same result in less time. But for now, this code is good enough.

    // Parallel specialized code.
    #ifdef GENESIS_OPENMP

        // We only need to calculate the upper triangle. Get the number of indices needed
        // to describe this triangle.
        size_t const max_k = utils::triangular_size( tree.node_count() );

        #pragma omp parallel for
        for( size_t k = 0; k < max_k; ++k ) {

            // For the given linear index, get the actual position in the Matrix.
            auto const rc = utils::triangular_indices( k, tree.node_count() );
            auto const r = rc.first;
            auto const c = rc.second;

            auto const& lca = lowest_common_ancestor( tree.node_at(r), tree.node_at(c) );
            res( r, c ) = lca.index();
            res( c, r ) = lca.index();
        }

        // Lastly, because the trinangular indices exluce the diagonale, we need to fill this
        // by hand. Luckily, those are always the indices themselves, as the LCA of a node
        // and itself is again itself.
        #pragma omp parallel for
        for( size_t d = 0; d < tree.node_count(); ++d ) {
            res( d, d ) = d;
        }

    // If no threads are available at all, use serial version.
    #else

        for( size_t r = 0; r < tree.node_count(); ++r ) {

            // The result is symmetric - we only calculate the upper triangle.
            for( size_t c = r; c < tree.node_count(); ++c ) {

                auto const& lca = lowest_common_ancestor( tree.node_at(r), tree.node_at(c) );
                res( r, c ) = lca.index();
                res( c, r ) = lca.index();
            }

            // See above: the diagonale contains its indices.
            assert( res( r, r ) == r );
        }

    #endif

    return res;
}

} // namespace tree
} // namespace genesis
