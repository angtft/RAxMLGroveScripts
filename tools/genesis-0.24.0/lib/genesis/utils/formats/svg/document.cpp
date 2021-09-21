/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

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

#include "genesis/utils/formats/svg/document.hpp"

#include "genesis/utils/core/options.hpp"
#include "genesis/utils/core/version.hpp"
#include "genesis/utils/formats/svg/helper.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/date_time.hpp"

namespace genesis {
namespace utils {

// =================================================================================================
//     Svg Document
// =================================================================================================

// -------------------------------------------------------------
//     Static Members
// -------------------------------------------------------------

std::string SvgDocument::indentation_string = "    ";

// -------------------------------------------------------------
//     Members
// -------------------------------------------------------------

SvgBox SvgDocument::bounding_box() const
{
    // Get bounding box of all elements and the dimensions of the document.
    SvgBox bbox;
    for( auto const& elem : content_ ) {
        bbox = SvgBox::combine( bbox, elem.bounding_box() );
    }
    return bbox;
}

/**
 * @brief Write the SvgDocument to an output stream.
 */
void SvgDocument::write( std::ostream& out ) const
{
    // Get a box around all elements, and use it to measure doc dimensions and shifting.
    auto bbox = bounding_box();
    double doc_width  = margin.left + bbox.width()  + margin.right;
    double doc_height = margin.top  + bbox.height() + margin.bottom;

    // SVG header.
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg";
    out << svg_attribute( "xmlns", "http://www.w3.org/2000/svg" );
    out << svg_attribute( "xmlns:xlink", "http://www.w3.org/1999/xlink" );
    out << svg_attribute( "width",  doc_width );
    out << svg_attribute( "height", doc_height );
    if( overflow != Overflow::kNone ) {
        out << overflow_to_string( overflow );
    }
    out << ">\n";

    // Some metadata.
    out << svg_comment(
        "Created with genesis " + genesis_version() + " (" + genesis_url() + ") " +
        "on " + current_date() + " at " + current_time()
    ) << "\n";
    if( Options::get().command_line_string() != "" ) {
        out << svg_comment( "Program invocation: " + Options::get().command_line_string() ) << "\n";
    }

    // Gradients and other definitions. Need to come before the content.
    if( ! defs.empty() ) {
        out << SvgDocument::indentation_string << "<defs>\n";
        for( auto const& def : defs ) {
            def.write( out, 2 );
        }
        out << SvgDocument::indentation_string << "</defs>\n";
    }

    // Options to hand over to all elements. Currently not needed, because we do the shifting
    // for the margin by using a group (see immediately below).
    auto options = SvgDrawingOptions();
    // options.offset_x = margin.top;
    // options.offset_y = margin.top;
    // options.offset_x = margin.left - bbox.top_left.x;
    // options.offset_y = margin.top  - bbox.top_left.y;

    // Main group for all elements. We use this to make the handling of the margin
    // easier.
    out << SvgDocument::indentation_string << "<g transform=\"translate( ";
    out << margin.left - bbox.top_left.x << ", " << margin.top - bbox.top_left.y;
    out << ")\" >\n";

    // Different approach, use inner svg element instead of group.
    // Didn't work well with some viewers.
    // out << SvgDocument::indentation_string << "<svg";
    // out << svg_attribute( "x", margin.left - bbox.top_left.x );
    // out << svg_attribute( "y", margin.top  - bbox.top_left.y );
    // out << ">\n";

    // Print content.
    for( auto const& elem : content_ ) {
        elem.write( out, 2, options );

        // Draw bounding boxes around all elements, for testing purposes.
        // auto bb = elem.bounding_box();
        // SvgRect( bb.top_left, bb.size(), SvgStroke(), SvgFill( Color(), 0.0 ) ).write( out, 1 );
    }

    // Close main grouping.
    // out << SvgDocument::indentation_string << "</svg>\n";
    out << SvgDocument::indentation_string << "</g>\n";

    // Finish.
    out << "</svg>\n";
}

SvgDocument& SvgDocument::add( SvgObject const& object )
{
    content_.push_back( object );
    return *this;
}

SvgDocument& SvgDocument::add( SvgObject&& object )
{
    content_.push_back( std::move( object ));
    return *this;
}

SvgDocument& SvgDocument::operator <<( SvgObject const& object )
{
    return add( object );
}

SvgDocument& SvgDocument::operator <<( SvgObject&& object )
{
    return add( std::move( object ));
}

std::string SvgDocument::overflow_to_string( SvgDocument::Overflow value )
{
    // switch( value ) {
    //     case Overflow::kNone:
    //         return std::string();
    //     case Overflow::kVisible:
    //         return svg_attribute( "style", "overflow: visible" );
    //     case Overflow::kHidden:
    //         return svg_attribute( "style", "overflow: hidden" );
    //     case Overflow::kScroll:
    //         return svg_attribute( "style", "overflow: scroll" );
    //     case Overflow::kAuto:
    //         return svg_attribute( "style", "overflow: auto" );
    //     case Overflow::kInherit:
    //         return svg_attribute( "style", "overflow: inherit" );
    //     default:
    //         throw std::invalid_argument(
    //             "Invalid Svg attribute Overflow for Svg Document."
    //         );
    // }
    switch( value ) {
        case Overflow::kNone:
            return std::string();
        case Overflow::kVisible:
            return svg_attribute( "overflow", "visible" );
        case Overflow::kHidden:
            return svg_attribute( "overflow", "hidden" );
        case Overflow::kScroll:
            return svg_attribute( "overflow", "scroll" );
        case Overflow::kAuto:
            return svg_attribute( "overflow", "auto" );
        case Overflow::kInherit:
            return svg_attribute( "overflow", "inherit" );
        default:
            throw std::invalid_argument(
                "Invalid Svg attribute Overflow for Svg Document."
            );
    }
}

} // namespace utils
} // namespace genesis
