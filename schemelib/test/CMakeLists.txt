set( ROSETTA_SRC $ENV{CMAKE_ROSETTA_PATH}/source )

link_directories( ${ROSETTA_BIN} )

set( ALL_ROSETTA_LIBS
basic cifparse core.1 core.2 core.3 core.4 core.5 cppdb devel xml2 numeric ObjexxFCL protocols.1 protocols.3 protocols.4 protocols.7 protocols.8 protocols_a.2 protocols_a.5 protocols_a.6 protocols_b.2 protocols_b.5 protocols_b.6 protocols_c.5 protocols_c.6 protocols_d.5 protocols_d.6 protocols_e.5 protocols_e.6 protocols_f.5 protocols_g.5 protocols_h.5 sqlite3 utility
 
#   basic core.1 core.2 core.3 core.4 core.5 devel numeric ObjexxFCL protocols.1 protocols.3 protocols.6 protocols.7 
#   protocols_a.2 protocols_a.4 protocols_a.5 protocols_b.2 protocols_b.4 protocols_b.5 protocols_c.4 protocols_c.5 
#   protocols_d.4 protocols_d.5 protocols_e.4 protocols_e.5 protocols_f.4 protocols_g.4 protocols_h.4 utility 
    cifparse cppdb sqlite3 z
)
include_directories( ${ROSETTA_SRC}/src )
include_directories( ${ROSETTA_SRC}/external/include )
include_directories( ${ROSETTA_SRC}/external/dbio )
include_directories( ${ROSETTA_SRC}/external/boost_1_55_0 )
include_directories( ${ROSETTA_SRC}/src/platform/linux )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCXX11 -DPTR_STD" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMULTI_THREADED" )

# shitty "data" files to compile faster
list( APPEND LIB_CC ../scheme/rosetta/score/EtableParams_init.lib.cc ../scheme/rosetta/atype_map.lib.cc ../scheme/chemical/ligand_factory.lib.cc )
SET_SOURCE_FILES_PROPERTIES( ${LIB_CC} PROPERTIES COMPILE_FLAGS -O0 )
add_library( scheme ${LIB_CC} )

FILE(GLOB L1 "../scheme/[0-9a-zA-Z_]*.cc")
FILE(GLOB L2 "../scheme/*/[0-9a-zA-Z_]*.cc")
FILE(GLOB L3 "../scheme/*/*/[0-9a-zA-Z_]*.cc")
FILE(GLOB L4 "../scheme/*/*/*/[0-9a-zA-Z_]*.cc")
FILE(GLOB L5 "../scheme/*/*/*/*/[0-9a-zA-Z_]*.cc")

#list(APPEND EXTRA_LIBS boost_system  )

include_directories(".")

add_executable( test_libscheme main_test.cc ${L1} ${L2} ${L3} ${L4} ${L5} )
target_link_libraries( test_libscheme scheme ${EXTRA_LIBS} )
# install ( TARGETS test_libscheme RUNTIME DESTINATION bin )

add_executable(quick_test_libscheme quick_test.cc  )
target_link_libraries(quick_test_libscheme scheme ${EXTRA_LIBS})
