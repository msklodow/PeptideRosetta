CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(COMPILER xlc)
SET(MODE release_bluegene)
SET(LINK_RULES static)

INCLUDE(../build/libraries.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/static.build.cmake)
INCLUDE(../build/build.settings.cmake)
#INCLUDE(../build/build.post.cmake)
INCLUDE(../build/external.static.cmake)

FOREACH(LIBRARY ${LIBRARIES})
	INCLUDE(../build/${LIBRARY}.cmake)
	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
ENDFOREACH( LIBRARY )

LINK_DIRECTORIES(
	/bgsys/drivers/ppcfloor/comm/xl/include/
	/soft/libraries/alcf/current/xl/ZLIB/lib
)

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DUSEMPI)
ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
SET(COMPILE_FLAGS "-O3")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")

#FOREACH(LIBRARY ${LIBRARIES})
#	INCLUDE(../build/${LIBRARY}.cmake)
#	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
#	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
#		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
#	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
#ENDFOREACH( LIBRARY )
##### Libraries

### Static Libraries (including rosetta)
### These MUST be full paths, otherwise binaries will not recompile when these change
### THIS ORDER IS VERY IMPORTANT -> NEEDS TO BE REVERSE ORDER OF DEPENDENCY GRAPH
SET(LINK_STATIC_LIBS
    ${CMAKE_CURRENT_SOURCE_DIR}/libdevel.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.7.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.6.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_e.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_d.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_c.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_h.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_g.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_f.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_e.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_d.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_c.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.3.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.2.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.2.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.1.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.6.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.4.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.3.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.2.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.1.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libbasic.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libnumeric.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libutility.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libObjexxFCL.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcppdb-static.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libcifparse.a
    -lz)
### libc must be linked dynamically if using static ompi libs
#SET(LINK_DYNAMIC_LIBS "-Wl,-Bdynamic -l c")

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic")

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic -lgcc_eh")
SET(TRAILING_STATIC "-Bstatic -L../../external/lib/")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ${TRAILING_STATIC})

### Change dynamic linker from default /lib/ld64.so.1 to /lib64/ld-2.5.so,
### since that file doesn't exist and path is hardcoded
SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--relax")



##### Top level binaries
# just loophash stuff for now
### loophash_createfiltereddb
#ADD_EXECUTABLE(loophash_createfiltereddb.mpi.linuxxlrelease ../../src/apps/pilot/kenjung/loophash_createfiltereddb.cc)
#TARGET_LINK_LIBRARIES(loophash_createfiltereddb.mpi.linuxxlrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(loophash_createfiltereddb.mpi.linuxxlrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(loophash_createfiltereddb.mpi.linuxxlrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
#ADD_DEPENDENCIES(loophash_createfiltereddb.mpi.linuxxlrelease ObjexxFCL utility numeric basic
#	core.1 core.2 core.3 core.4 core.5 core.6 protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols.6 protocols.7 devel)

### design_cycpeptide_MPI
#ADD_EXECUTABLE(design_cycpeptide_MPI.mpi.linuxxlrelease ../../src/apps/pilot/vmullig/design_cycpeptide_MPI.cc)
#TARGET_LINK_LIBRARIES(design_cycpeptide_MPI.mpi.linuxxlrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(design_cycpeptide_MPI.mpi.linuxxlrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(design_cycpeptide_MPI.mpi.linuxxlrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
#ADD_DEPENDENCIES(design_cycpeptide_MPI.mpi.linuxxlrelease cifparse ObjexxFCL utility numeric
#basic core.1 core.2 core.3 core.4 core.5 core.6 protocols.1 protocols_a.2
#protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4
#protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4
#protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5
#protocols.6 protocols.7 devel)

### rosetta_scripts
ADD_EXECUTABLE(rosetta_scripts.mpi.linuxxlrelease ../../src/apps/public/rosetta_scripts/rosetta_scripts.cc)
TARGET_LINK_LIBRARIES(rosetta_scripts.mpi.linuxxlrelease "${LINK_ALL_LIBS}" )
SET_TARGET_PROPERTIES(rosetta_scripts.mpi.linuxxlrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(rosetta_scripts.mpi.linuxxlrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(rosetta_scripts.mpi.linuxxlrelease ObjexxFCL cifparse utility numeric
basic core.1 core.2 core.3 core.4 core.5 core.6 protocols.1 protocols_a.2
protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4
protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4
protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5
protocols.6 protocols.7 devel)

### loophash_mpi
#ADD_EXECUTABLE(loophash_mpi.mpi.linuxxlrelease ../../src/apps/pilot/mike/loophash_mpi.cc)
#TARGET_LINK_LIBRARIES(loophash_mpi.mpi.linuxxlrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(loophash_mpi.mpi.linuxxlrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(loophash_mpi.mpi.linuxxlrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
#ADD_DEPENDENCIES(loophash_mpi.mpi.linuxxlrelease ObjexxFCL utility numeric basic
#	core.1 core.2 core.3 core.4 core.5 core.6 protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols.6 protocols.7 devel)

#### minirosetta
#ADD_EXECUTABLE(minirosetta.mpi.linuxxlrelease ../../src/apps/public/boinc/minirosetta.cc)
#TARGET_LINK_LIBRARIES(minirosetta.mpi.linuxxlrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(minirosetta.mpi.linuxxlrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(minirosetta.mpi.linuxxlrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
## Dependencies are set so that binaries are created after rosetta libraries
## are built
#ADD_DEPENDENCIES(minirosetta.mpi.linuxxlrelease ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 core.6 protocols devel)
