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
 * @brief Implementation of Options members.
 *
 * @file
 * @ingroup main
 */

#include "genesis/utils/core/options.hpp"

#include "genesis/utils/core/version.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>

#if defined( _WIN32 ) || defined(  _WIN64  )
#   include <io.h>
#   include <windows.h>
#else
#   include <cpuid.h>
#   include <stdio.h>
#   include <sys/ioctl.h>
#   include <unistd.h>
#endif

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

#ifdef GENESIS_PTHREADS
#    include <thread>
#endif

namespace genesis {
namespace utils {

// =================================================================================================
//     Initialization
// =================================================================================================

#if defined( DEBUG ) && defined( NDEBUG )
    static_assert( false, "Cannot compile with both DEBUG and NDEBUG flags set." );
#endif

#if ! defined( DEBUG ) && ! defined( NDEBUG )
    static_assert( false, "Cannot compile with neiher DEBUG nor NDEBUG flag set." );
#endif

Options::Options()
{
    // Initialize number of threads to hardware cores.
    number_of_threads( guess_number_of_threads() );

    // Initialize random seed with time.
    random_seed( std::chrono::system_clock::now().time_since_epoch().count() );
}

// =================================================================================================
//     Command Line
// =================================================================================================

std::string Options::command_line_string () const
{
    std::string ret = "";
    for (size_t i = 0; i < command_line_.size(); ++i) {
        std::string a = command_line_[i];
        ret += (i==0 ? "" : " ") + a;
    }
    return ret;
}

void Options::command_line( int const argc, char const* const* argv )
{
    // Store all arguments in the array.
    command_line_.clear();
    for (int i = 0; i < argc; i++) {
        command_line_.push_back(argv[i]);
    }
}

// =================================================================================================
//     Number of Threads
// =================================================================================================

void Options::number_of_threads ( unsigned int number )
{
    if( number == 0 ) {
        #ifdef GENESIS_PTHREADS
            number = std::thread::hardware_concurrency();
            if( number == 0 ) {
                number = 1;
            }
        #else
            number = 1;
        #endif
    }
    number_of_threads_ = number;

    #if defined( GENESIS_OPENMP )

        // If we use OpenMp, set the thread number there, too.
        omp_set_num_threads( number );

    #endif
}

bool Options::hyperthreads_enabled() const
{
    // Get CPU info.
    int32_t info[4];
    #ifdef _WIN32
        __cpuid( info, 1 );
    #else
        __cpuid_count( 1, 0, info[0], info[1], info[2], info[3] );
    #endif

    return (bool) (info[3] & (0x1 << 28));
}

unsigned int Options::guess_number_of_threads( bool use_openmp ) const
{
    // Dummy to avoid compiler warnings.
    (void) use_openmp;

    #if defined( GENESIS_OPENMP )

        // Use number of OpenMP threads, which might be set through the
        // `OMP_NUM_THREADS` environment variable.
        if( use_openmp ) {
            return omp_get_max_threads();
        }

    #endif

    #if defined( GENESIS_PTHREADS )

        // Initialize threads with actual number of cores.
        // return std::thread::hardware_concurrency();

        auto const lcores = std::thread::hardware_concurrency();
        auto const threads_per_core = hyperthreads_enabled() ? 2 : 1;
        return lcores / threads_per_core;

    #else

        // Default to single threaded.
        return 1;

    #endif
}

// =================================================================================================
//     Random Seed & Engine
// =================================================================================================

void Options::random_seed(const unsigned long seed)
{
    random_seed_ = seed;
    random_engine_.seed( seed );
}

// =================================================================================================
//     Run Time Environment
// =================================================================================================

bool Options::stdin_is_terminal()
{
    // Using http://stackoverflow.com/a/1312957/4184258
    #if defined( _WIN32 ) || defined(  _WIN64  )
        return _isatty( _fileno( stdin ));
    #else
        return isatty( fileno( stdin ));
    #endif
}

bool Options::stdout_is_terminal()
{
    #if defined( _WIN32 ) || defined(  _WIN64  )
        return _isatty( _fileno( stdout ));
    #else
        return isatty( fileno( stdout ));
    #endif
}

bool Options::stderr_is_terminal()
{
    #if defined( _WIN32 ) || defined(  _WIN64  )
        return _isatty( _fileno( stderr ));
    #else
        return isatty( fileno( stderr ));
    #endif
}

std::pair<int, int> Options::terminal_size()
{
    #if defined( _WIN32 ) || defined(  _WIN64  )

        CONSOLE_SCREEN_BUFFER_INFO csbi;
        int cols, rows;
        GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
        cols = csbi.srWindow.Right - csbi.srWindow.Left + 1;
        rows = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;
        return { cols, rows };

    #else

        struct winsize w;
        ioctl( STDOUT_FILENO, TIOCGWINSZ, &w );
        return { w.ws_col, w.ws_row };

    #endif
}

// =================================================================================================
//     Compile Time Environment
// =================================================================================================

bool Options::is_debug()
{
    #ifdef DEBUG
        return true;
    #else
        return false;
    #endif
}

bool Options::is_release()
{
    #ifdef NDEBUG
        return true;
    #else
        return false;
    #endif
}

std::string Options::build_type()
{
    #if defined( DEBUG )
        return "debug";
    #elif defined( NDEBUG )
        return "release";
    #else
        return "unknown";
    #endif
}

bool Options::is_little_endian()
{
    static const uint16_t e = 0x1000;
    return 0 == *reinterpret_cast< uint8_t const* >( &e );
}

bool Options::is_big_endian()
{
    static const uint16_t e = 0x0001;
    return 0 == *reinterpret_cast< uint8_t const* >( &e );
}

std::string Options::platform()
{
    #if defined _WIN64
        return "Win64";
    #elif defined _WIN32
        return "Win32";
    #elif defined __linux__
        return "Linux";
    #elif defined __APPLE__
        return "Apple";
    #elif defined __unix__
        return "Unix";
    #else
        return "Unknown";
    #endif
}

std::string Options::compiler_family()
{
    #if defined(__clang__)
        return "clang";
    #elif defined(__ICC) || defined(__INTEL_COMPILER)
        return "icc";
    #elif defined(__GNUC__) || defined(__GNUG__)
        return "gcc";
    #elif defined(__HP_cc) || defined(__HP_aCC)
        return "hp";
    #elif defined(__IBMCPP__)
        return "ilecpp";
    #elif defined(_MSC_VER)
        return "msvc";
    #elif defined(__PGI)
        return "pgcpp";
    #elif defined(__SUNPRO_CC)
        return "sunpro";
    #else
        return "unknown";
    #endif
}

std::string Options::compiler_version()
{
    #if defined(__clang__)
        return __clang_version__;
    #elif defined(__ICC) || defined(__INTEL_COMPILER)
        return __INTEL_COMPILER;
    #elif defined(__GNUC__) || defined(__GNUG__)
        return std::to_string(__GNUC__)            + "." +
               std::to_string(__GNUC_MINOR__)      + "." +
               std::to_string(__GNUC_PATCHLEVEL__)
        ;
    #elif defined(__HP_cc) || defined(__HP_aCC)
        return "";
    #elif defined(__IBMCPP__)
        return __IBMCPP__;
    #elif defined(_MSC_VER)
        return _MSC_VER;
    #elif defined(__PGI)
        return __PGI;
    #elif defined(__SUNPRO_CC)
        return __SUNPRO_CC;
    #else
        return "unknown";
    #endif
}

std::string Options::cpp_version()
{
    #ifdef __cplusplus
        return std::to_string(__cplusplus);
    #else
        return "unknown";
    #endif
}

std::string Options::compile_date_time()
{
    return std::string( __DATE__ " " __TIME__ );
}

bool Options::using_pthreads()
{
    #ifdef GENESIS_PTHREADS
        return true;
    #else
        return false;
    #endif
}

bool Options::using_openmp()
{
    #ifdef GENESIS_OPENMP
        return true;
    #else
        return false;
    #endif
}

bool Options::using_zlib()
{
    #ifdef GENESIS_ZLIB
        return true;
    #else
        return false;
    #endif
}

// =================================================================================================
//     Dump & Overview
// =================================================================================================

std::string Options::info() const
{
    std::string res;
    res += genesis_header() + "\n";
    res += info_compile_time() + "\n";
    res += info_run_time() + "\n";
    return res;
}

std::string Options::info_compile_time() const
{
    std::string res;
    res += "Compile Time Options\n";
    res += "=============================================\n\n";
    res += "Platform:          " + platform() + "\n";
    res += "Compiler:          " + compiler_family() + " " + compiler_version() + "\n";
    res += "C++ version:       " + cpp_version() + "\n";
    res += "Build type:        " + build_type()  + "\n";
    res += "Endianness:        " + std::string( is_little_endian() ? "little endian" : "big endian" ) + "\n";
    res += "Using Pthreads:    " + std::string( using_pthreads() ? "true" : "false" ) + "\n";
    res += "Using OpenMP:      " + std::string( using_openmp() ? "true" : "false" ) + "\n";
    return res;
}

std::string Options::info_run_time() const
{
    std::string res;
    res += "Run Time Options\n";
    res += "=============================================\n\n";
    auto const cli_str = command_line_string();
    res += "Command line:      " + ( cli_str.size() > 0 ? cli_str : "(not available)" ) + "\n";
    res += "Number of threads: " + std::to_string( number_of_threads() ) + "\n";
    res += "Random seed:       " + std::to_string( random_seed_ ) + "\n";
    return res;
}

} // namespace utils
} // namespace genesis
