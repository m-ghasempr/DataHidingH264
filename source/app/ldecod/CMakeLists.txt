# executable
set( EXE_NAME ldecod )

# get source files
file( GLOB ENC_SRC_FILES "*.c" )
file( GLOB COMMON_SRC_FILES "../../lib/lcommon/*.c" )

set ( SRC_FILES ${ENC_SRC_FILES} ${COMMON_SRC_FILES} )

# get include files
file( GLOB ENC_INC_FILES "*.h" )
file( GLOB COMMON_INC_FILES "../../lib/lcommon/*.h" )

set ( INC_FILES ${ENC_INC_FILES} ${COMMON_INC_FILES} )

# get additional libs for gcc on Ubuntu systems
if( CMAKE_SYSTEM_NAME STREQUAL "Linux" )
  if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
    if( USE_ADDRESS_SANITIZER )
      set( ADDITIONAL_LIBS asan )
    endif()
  endif()
endif()

# NATVIS files for Visual Studio
if( MSVC )
  file( GLOB NATVIS_FILES "../../VisualStudio/*.natvis" )
  # extend the stack size on windows to 2MB
  set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} /STACK:0x200000" )
endif()

# add executable
add_executable( ${EXE_NAME} ${SRC_FILES} ${INC_FILES} ${NATVIS_FILES} )
include_directories(${CMAKE_CURRENT_BINARY_DIR} . ../../lib/lcommon)

if( SET_ENABLE_TRACING )
  if( ENABLE_TRACING )
    target_compile_definitions( ${EXE_NAME} PUBLIC ENABLE_TRACING=1 )
  else()
    target_compile_definitions( ${EXE_NAME} PUBLIC ENABLE_TRACING=0 )
  endif()
endif()

if( CMAKE_COMPILER_IS_GNUCC AND BUILD_STATIC )
  set( ADDITIONAL_LIBS ${ADDITIONAL_LIBS} -static -static-libgcc )
  target_compile_definitions( ${EXE_NAME} PUBLIC ENABLE_WPP_STATIC_LINK=1 )
endif()

if(NOT MSVC)
  target_link_libraries( ${EXE_NAME} m Threads::Threads ${ADDITIONAL_LIBS} )
else()
  target_link_libraries( ${EXE_NAME} WS2_32 Threads::Threads ${ADDITIONAL_LIBS} )
endif()

# lldb custom data formatters
if( XCODE )
  add_dependencies( ${EXE_NAME} Install${PROJECT_NAME}LldbFiles )
endif()

if( CMAKE_SYSTEM_NAME STREQUAL "Linux" )
  add_custom_command( TARGET ${EXE_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
                                                          $<$<CONFIG:Debug>:${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}/ldecod>
                                                          $<$<CONFIG:Release>:${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}/ldecod>
                                                          $<$<CONFIG:RelWithDebInfo>:${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO}/ldecod>
                                                          $<$<CONFIG:MinSizeRel>:${CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL}/ldecod>
                                                          $<$<CONFIG:Debug>:${CMAKE_SOURCE_DIR}/bin/ldecod_static>
                                                          $<$<CONFIG:Release>:${CMAKE_SOURCE_DIR}/bin/ldecod_static>
                                                          $<$<CONFIG:RelWithDebInfo>:${CMAKE_SOURCE_DIR}/bin/ldecod_static>
                                                          $<$<CONFIG:MinSizeRel>:${CMAKE_SOURCE_DIR}/bin/ldecod_static> )
endif()

# example: place header files in different folders
source_group( "Natvis Files" FILES ${NATVIS_FILES} )

# set the folder where to place the projects
set_target_properties( ${EXE_NAME}  PROPERTIES FOLDER app LINKER_LANGUAGE C )

