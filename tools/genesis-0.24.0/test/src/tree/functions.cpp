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
 * @brief Testing Newick class.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include <string>
#include <vector>

#include "genesis/tree/common_tree/functions.hpp"
#include "genesis/tree/common_tree/newick_reader.hpp"
#include "genesis/tree/common_tree/tree.hpp"
#include "genesis/tree/formats/newick/reader.hpp"
#include "genesis/tree/function/functions.hpp"
#include "genesis/tree/tree.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"

using namespace genesis;
using namespace tree;

// =================================================================================================
//     Tree Sides
// =================================================================================================

void test_print_tree_sides_matrix( utils::Matrix<signed char> const& mat )
{
    // LOG_DBG << mat;
    std::stringstream os;
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            // os << mat(i, j);
            if( mat(i, j) == 0 ) {
                os << " 0";
            } else if( mat(i, j) == 1 ) {
                os << " 1";
            } else if( mat(i, j) == -1 ) {
                os << "-1";
            } else {
                os << " x";
            }

            if (j < mat.cols() - 1) {
                os << " ";
            }
        }
        os << "\n";
    }
    LOG_DBG << os.str();
}

TEST( TreeFunctions, EdgeSides )
{
    std::string const input = "((B,(D,E)C)A,F,(H,I)G)R;";
    Tree const tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto const edge_side_mat = edge_sides( tree );

    auto const exp = utils::Matrix<signed char>( 9, 9, {
        0,  1,  1, -1, -1, -1, -1, -1, -1,
       -1,  0, -1, -1, -1, -1, -1, -1, -1,
       -1, -1,  0, -1, -1, -1, -1, -1, -1,
       -1, -1, -1,  0, -1, -1, -1, -1, -1,
       -1, -1, -1, -1,  0,  1,  1,  1,  1,
       -1, -1, -1, -1, -1,  0,  1,  1, -1,
       -1, -1, -1, -1, -1, -1,  0, -1, -1,
       -1, -1, -1, -1, -1, -1, -1,  0, -1,
       -1, -1, -1, -1, -1, -1, -1, -1,  0
    });

    EXPECT_EQ( exp, edge_side_mat );

    // test_print_tree_sides_matrix( edge_side_mat );
}

TEST( TreeFunctions, NodeRootDirections )
{
    std::string const input = "((B,(D,E)C)A,F,(H,I)G)R;";
    Tree const tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto const node_root_mat = node_root_direction_matrix( tree );

    auto const exp = utils::Matrix<signed char>( 10, 10, {
        0, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        1,  0, -1, -1,  1,  1,  1,  1,  1,  1,
        1,  1,  0,  1,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  0,  1,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  0,  1,  1,  1,  1,  1,
        1,  1,  1,  1,  1,  0, -1, -1, -1, -1,
        1,  1,  1,  1,  1,  1,  0, -1, -1,  1,
        1,  1,  1,  1,  1,  1,  1,  0,  1,  1,
        1,  1,  1,  1,  1,  1,  1,  1,  0,  1,
        1,  1,  1,  1,  1,  1,  1,  1,  1,  0,
    });

    EXPECT_EQ( exp, node_root_mat );

    // test_print_tree_sides_matrix( node_root_mat );
}

// =================================================================================================
//     Subtree Size
// =================================================================================================

void TestSubtreeSize( size_t link_index, size_t out_subtree_size )
{
    std::string input = "((B,(D,E)C)A,F,(H,I)G)R;";

    Tree tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto st_size = subtree_size( tree, tree.link_at( link_index ));
    EXPECT_EQ( out_subtree_size, st_size ) << " with link index " << link_index;
}

TEST( TreeFunctions, SubtreeSize )
{
    TestSubtreeSize(  0, 5 );
    TestSubtreeSize(  1, 1 );
    TestSubtreeSize(  2, 3 );
    TestSubtreeSize(  3, 7 );
    TestSubtreeSize(  4, 1 );
    TestSubtreeSize(  5, 1 );
    TestSubtreeSize(  6, 9 );
    TestSubtreeSize(  7, 9 );
    TestSubtreeSize(  8, 9 );
    TestSubtreeSize(  9, 5 );
    TestSubtreeSize( 10, 1 );
    TestSubtreeSize( 11, 3 );
    TestSubtreeSize( 12, 7 );
    TestSubtreeSize( 13, 1 );
    TestSubtreeSize( 14, 1 );
    TestSubtreeSize( 15, 9 );
    TestSubtreeSize( 16, 9 );
    TestSubtreeSize( 17, 9 );
}

// =================================================================================================
//     Subtree Sizes
// =================================================================================================

void TestSubtreeSizes( std::string node_name, std::vector<size_t> out_sizes )
{
    std::string input = "((B,(D,E)C)A,F,(H,I)G)R;";

    Tree tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto node = find_node( tree, node_name );
    ASSERT_NE( nullptr, node );

    auto st_sizes = subtree_sizes( tree, *node );
    EXPECT_EQ( out_sizes, st_sizes ) << " with start node " << node_name;
}

TEST( TreeFunctions, SubtreeSizes )
{
    TestSubtreeSizes( "R", { 9, 2, 0, 0, 0, 4, 2, 0, 0, 0 } );
    TestSubtreeSizes( "A", { 4, 2, 0, 0, 0, 9, 2, 0, 0, 0 } );
    TestSubtreeSizes( "B", { 4, 2, 0, 0, 0, 8, 2, 0, 0, 9 } );
    TestSubtreeSizes( "C", { 4, 2, 0, 0, 0, 6, 9, 0, 0, 0 } );
    TestSubtreeSizes( "D", { 4, 2, 0, 0, 0, 6, 8, 0, 9, 0 } );
    TestSubtreeSizes( "E", { 4, 2, 0, 0, 0, 6, 8, 9, 0, 0 } );
    TestSubtreeSizes( "F", { 8, 2, 0, 0, 9, 4, 2, 0, 0, 0 } );
    TestSubtreeSizes( "G", { 6, 9, 0, 0, 0, 4, 2, 0, 0, 0 } );
    TestSubtreeSizes( "H", { 6, 8, 0, 9, 0, 4, 2, 0, 0, 0 } );
    TestSubtreeSizes( "I", { 6, 8, 9, 0, 0, 4, 2, 0, 0, 0 } );
}

// =================================================================================================
//     Subtree Max Path Height
// =================================================================================================

void TestSubtreeMaxPathHeight( std::string node_name, size_t out_size )
{
    std::string input = "((B,(D,E)C)A,F,(H,I)G)R;";

    Tree tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto node = find_node( tree, node_name );
    ASSERT_NE( nullptr, node );

    // We are lazy and only evaluate the links towards the root.
    auto st_size = subtree_max_path_height( tree, node->link().outer() );
    EXPECT_EQ( out_size, st_size ) << " with start node " << node_name;
}

TEST( TreeFunctions, SubtreeMaxPathHeight )
{
    // TestSubtreeMaxPathHeight( "R", 3 );
    TestSubtreeMaxPathHeight( "A", 2 );
    TestSubtreeMaxPathHeight( "B", 0 );
    TestSubtreeMaxPathHeight( "C", 1 );
    TestSubtreeMaxPathHeight( "D", 0 );
    TestSubtreeMaxPathHeight( "E", 0 );
    TestSubtreeMaxPathHeight( "F", 0 );
    TestSubtreeMaxPathHeight( "G", 1 );
    TestSubtreeMaxPathHeight( "H", 0 );
    TestSubtreeMaxPathHeight( "I", 0 );
}

TEST( TreeFunctions, SubtreeMaxPathHeights )
{
    std::string input = "((B,(D,E)C)A,F,(H,I)G)R;";

    Tree tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto heights     = subtree_max_path_heights( tree );
    auto exp_heights = std::vector<size_t>({ 3, 1, 0, 0, 0, 2, 1, 0, 0, 0 });
    EXPECT_EQ( exp_heights, heights );
}

// =================================================================================================
//     Misc
// =================================================================================================

void TestTreeLCA( std::string node_name_a, std::string node_name_b, std::string node_name_lca )
{
    std::string input = "((B,(D,E)C)A,F,(H,I)G)R;";

    Tree tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto node_a = find_node( tree, node_name_a );
    auto node_b = find_node( tree, node_name_b );
    ASSERT_NE( nullptr, node_a );
    ASSERT_NE( nullptr, node_b );

    auto& node_lca = lowest_common_ancestor( *node_a, *node_b );
    EXPECT_EQ( node_name_lca, node_lca.data<CommonNodeData>().name )
        << " with nodes " << node_name_a << ", " << node_name_b << " and LCA " << node_name_lca;
}

TEST( TreeFunctions, LCA )
{
    TestTreeLCA( "A", "A", "A" );
    TestTreeLCA( "A", "B", "A" );
    TestTreeLCA( "A", "F", "R" );
    TestTreeLCA( "E", "C", "C" );
    TestTreeLCA( "E", "H", "R" );
    TestTreeLCA( "H", "I", "G" );
}

TEST( TreeFunctions, LCAs )
{
    std::string input = "((B,(D,E)C)A,F,(H,I)G)R;";
    Tree const tree = CommonTreeNewickReader().read( utils::from_string( input ));

    auto const lcas = lowest_common_ancestors( tree );
    auto const exp = utils::Matrix<size_t>( 10, 10, {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 2, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 3, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 4, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 5, 5, 5, 5, 5,
        0, 0, 0, 0, 0, 5, 6, 6, 6, 5,
        0, 0, 0, 0, 0, 5, 6, 7, 6, 5,
        0, 0, 0, 0, 0, 5, 6, 6, 8, 5,
        0, 0, 0, 0, 0, 5, 5, 5, 5, 9
    });

    EXPECT_EQ( exp, lcas );
}

TEST( TreeFunctions, SignMatrix )
{
    // Skip test if no data availabe.
    NEEDS_TEST_DATA;

    // Read and process tree.
    std::string infile = environment->data_dir + "tree/rooted.newick";
    Tree tree = CommonTreeNewickReader().read( utils::from_file( infile ));

    // Full sign matrix
    auto const exf = utils::Matrix<signed char>( 9, 9, {
        0, -1, -1, -1, +1, +1, +1, +1, +1,
        0,  0, -1, +1,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0, -1, -1, -1, +1,
        0,  0,  0,  0,  0,  0, -1, +1,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0
    });
    auto const smf = sign_matrix( tree );
    EXPECT_EQ( exf, smf );

    // LOG_DBG << smf;
    // LOG_DBG;
    // LOG_DBG << exf;

    // Compressed sign matrix
    auto const exc = utils::Matrix<signed char>( 4, 5, {
        -1, -1, +1, +1, +1,
        -1, +1,  0,  0,  0,
         0,  0, -1, -1, +1,
         0,  0, -1, +1,  0
    });
    auto const smc = sign_matrix( tree, true );
    EXPECT_EQ( exc, smc );

    // LOG_DBG << smc;
    // LOG_DBG;
    // LOG_DBG << exc;
}
