#
# @BEGIN LICENSE
#
# libresponse_psi4 by Eric Berquist, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(libresponse_psi4 CXX)

find_package(psi4 1.1 REQUIRED)

# set(SRC
#     libresponse/configurable.C
#     libresponse/dump_ao_integrals.C
#     libresponse/index_printing.C
#     libresponse/indices.C
#     libresponse/matvec_i.C
#     libresponse/utils.C
#     libresponse/linear/helpers.C
#     libresponse/linear/interface.C
#     libresponse/linear/interface_nonorthogonal.C
#     libresponse/linear/printing.C
#     libresponse/operator_spec.C
#     libresponse/set_defaults.C
#     )

set(SRC_PSI4
    matvec_psi4.cc
    # operatordatatype.cc
    parse_operators.cc
    wrappers.cc
    plugin.cc
    )

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(libresponse)

function(add_psi4_plugin_custom TARGET SOURCES)
    # remove ${TARGET} from ${ARGV} to use ${ARGV} as SOURCES
    list(REMOVE_AT ARGV 0)

    set(psi4_CXX_STANDARD 14)
    # checks compiler (and gcc, if necessary) C++14 compliant
    include(custom_cxxstandard)

    add_library(${TARGET} MODULE ${ARGV})
    target_link_libraries(${TARGET} PRIVATE psi4::core response)
    target_link_libraries(${TARGET} PRIVATE tgt::MathOpenMP)
    set_target_properties(${TARGET} PROPERTIES PREFIX "")
    if(APPLE)
        set_target_properties(${TARGET} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()
endfunction()

add_psi4_plugin_custom(libresponse_psi4 ${SRC_PSI4})
