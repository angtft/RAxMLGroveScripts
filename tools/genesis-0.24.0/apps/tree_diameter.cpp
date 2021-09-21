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
 * This is the demo "Compare Jplace Files". See the Manual for more information.
 */

#include "genesis/genesis.hpp"

#include <string>
#include <unordered_map>
#include <numeric>
#include <algorithm>

using namespace genesis;

/**
 * @brief Main function that processes two jplace files and compares them.
 *
 * This program is useful to compare two `jplace` files that were computed with different settings
 * or even different programs (EPA or pplacer). It is not meant for evaluating differences in the
 * microbial communities that are represented by the input. Instead, it is meant for files that
 * share @link Pquery Pqueries@endlink (identified by their names), and gives information about
 * differences between the Placements in those Pqueries.
 *
 * The program takes two input `jplace` file paths as input. It compares the
 * @link Pquery Pqueries@endlink and their Placements and prints two tables:
 *
 *   1. An overview table that lists all Pqueries of the two files that have a PqueryName in common.
 *      This table indicates whether the top Placement (the one with the highest `like_weight_ratio`)
 *      of both Pqueries is the same (i.e., is located at the same branch); it furthermore indicates
 *      whether all Placements (sorted by their `like_weight_ratio`) are the same, that is, if they
 *      are located on the same branches. Lastly, the difference in log-likelhood and the
 *      @link earth_movers_distance() Earth Movers Distance@endlink between the Pqueries is printed.
 *   2. A detail table that lists all Placements of the Pqueries that were marked invalid in the
 *      overview table - that is, if either the top rank or any other placement was not equally
 *      placed in a Pquery. This table lists the Placements for such Pqueries, sorted by their
 *      `like_weight_ratio`, and shows on which branches (edge_num) they are placed in the two
 *      Pqueries. If the Placements are on the same branch, they are considered correct.
 *
 * The program expects that the reference Tree%s of the input are topologically identical. In order
 * to compensate for differences in branch lengths, both Trees are normalized in the beginning, so
 * that their length (sum of branch lengths) is 1.0. This also means that the Earth Movers Distance
 * yields comparable values in the range `[ 0.0, 1.0 ]`.
 */
int main( int argc, const char* argv[] )
{
    // -----------------------------------------------------
    //     Init and Settings.
    // -----------------------------------------------------

    using namespace ::genesis::placement;
    using namespace ::genesis::tree;
    using namespace ::genesis::utils;
    
    Tree orig_tree = CommonTreeNewickReader().read(from_file(argv[1]));
    
    auto printer = PrinterCompact();
    
    //std::cout << printer.print(orig_tree) << std::endl;
    //std::cout << "Tree diameter: " <<  << "\n";
    std::cout << length(orig_tree) << " " << diameter(orig_tree) << "\n";

    return 0;
}
