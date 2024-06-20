
# attemp to use a scons build...
#set( ROSETTA_SRC /work/sheffler/ssd/rosetta/master/source )
#set( ROSETTA_BIN release/linux/3.10/64/x86/gcc/6.1/cxx11-omp )
#link_directories( ${ROSETTA_SRC}/build/src/${ROSETTA_BIN} )
#link_directories( ${ROSETTA_SRC}/build/external/${ROSETTA_BIN} )


set( ROSETTA_SRC $ENV{CMAKE_ROSETTA_PATH}/source )

link_directories( ${ROSETTA_BIN} )

set( ALL_ROSETTA_LIBS
basic cifparse core.1 core.2 core.3 core.4 core.5 cppdb devel xml2 numeric ObjexxFCL protocols.1 protocols.3 protocols.4 protocols.7 protocols.8 protocols_a.2 protocols_a.5 protocols_a.6 protocols_b.2 protocols_b.5 protocols_b.6 protocols_c.5 protocols_c.6 protocols_d.5 protocols_d.6 protocols_e.5 protocols_e.6 protocols_f.5 protocols_g.5 protocols_h.5 sqlite3 utility
 
#	basic core.1 core.2 core.3 core.4 core.5 devel numeric ObjexxFCL protocols.1 protocols.3 protocols.6 protocols.7 
#	protocols_a.2 protocols_a.4 protocols_a.5 protocols_b.2 protocols_b.4 protocols_b.5 protocols_c.4 protocols_c.5 
#	protocols_d.4 protocols_d.5 protocols_e.4 protocols_e.5 protocols_f.4 protocols_g.4 protocols_h.4 utility 
	cifparse cppdb sqlite3 z
)
include_directories( ${ROSETTA_SRC}/src )
include_directories( ${ROSETTA_SRC}/external/include )
include_directories( ${ROSETTA_SRC}/external/dbio )
include_directories( ${ROSETTA_SRC}/external/boost_1_55_0 )
include_directories( ${ROSETTA_SRC}/src/platform/linux )
include_directories( . )


SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCXX11 -DPTR_STD" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMULTI_THREADED" )

add_subdirectory( python )

# this omp shit breaks python bindings! do it after
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP -fopenmp" )
LIST( APPEND ALL_ROSETTA_LIBS gomp )

if ( ${USEHDF5} )
    FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C CXX HL)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
    link_directories( ${HDF5_LIBRARY_DIRS} )
    LIST( APPEND ALL_ROSETTA_LIBS ${HDF5_LIBRARIES}) #hdf5 hdf5_hl hdf5_cpp hdf5_hl_cp )
endif()


add_subdirectory( riflib )

set( EXES "test_librosetta" "rifgen" "rif_dock_test" "scheme_make_bounding_grids" )
foreach( EXE ${EXES} )
	message( "riflib exe: " ${EXE} )

	add_executable( ${EXE} ${EXE}.cc  )
	target_link_libraries ( ${EXE} ${ALL_ROSETTA_LIBS} riflib )
	install ( TARGETS ${EXE} RUNTIME DESTINATION bin )

endforeach( EXE )
