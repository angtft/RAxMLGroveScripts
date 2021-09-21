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
 * @ingroup test
 */

#include "src/common.hpp"

#include "genesis/utils/tools/char_lookup.hpp"

using namespace genesis::utils;

TEST(CharLookup, Simple)
{
    auto cl = CharLookup<bool>();

    cl.set_selection( "abc", true );
    EXPECT_TRUE(  cl.get('a') );
    EXPECT_FALSE( cl.get('A') );

    cl.set_range( 'G', 'L', true );
    EXPECT_FALSE( cl.get('F') );
    EXPECT_TRUE(  cl.get('I') );

    cl.set_char( 'I', false );
    EXPECT_FALSE( cl.get('I') );
}
