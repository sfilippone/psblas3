include(CMakeParseArguments)

# Function to parse version info from git and/or .VERSION file
function(set_version)
  set(options "")
  set(oneValueArgs VERSION_VARIABLE GIT_DESCRIBE_VAR CUSTOM_VERSION_FILE CUSTOM_VERSION_REGEX )
  set(multiValueArgs "")
  cmake_parse_arguments(set_version "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Algorithm:
  # 1. Get first line of .VERSION file, which will be set via `git archive` so long as
  #    
  # 2. If not a packaged release check if this is an active git repo
  # 3. Get version info from `git describe`
  # 4. First the most recent tag is fetched if available
  # 5. Then the full `git describe` output is fetched


  if(NOT set_version_CUSTOM_VERSION_REGEX)
    set(_VERSION_REGEX "[vV]*[0-9]+\\.[0-9]+\\.[0-9]+")
  else()
    set(_VERSION_REGEX ${set_version_CUSTOM_VERSION_REGEX})
  endif()
  if(NOT set_version_CUSTOM_VERSION_FILE)
    set(_VERSION_FILE "${CMAKE_SOURCE_DIR}/.VERSION")
  else()
    set(_VERSION_FILE "${set_version_CUSTOM_VERSION_FILE}")
  endif()

  file(STRINGS "${_VERSION_FILE}" first_line
    LIMIT_COUNT 1
    )

  string(REGEX MATCH ${_VERSION_REGEX}
    _package_version "${first_line}")

  if((NOT (_package_version MATCHES ${_VERSION_REGEX})) AND (EXISTS "${CMAKE_SOURCE_DIR}/.git"))
    message( STATUS "Build from git repository detected")
    find_package(Git)
    if(GIT_FOUND)
      set(GIT_FOUND "${GIT_FOUND}" PARENT_SCOPE)
      execute_process(COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
	RESULT_VARIABLE _git_status
	OUTPUT_VARIABLE _git_output
	OUTPUT_STRIP_TRAILING_WHITESPACE)
      if((_git_status STREQUAL "0") AND (_git_output MATCHES ${_VERSION_REGEX}))
	set(_package_version "${_git_output}")
      endif()
      execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
	RESULT_VARIABLE _git_status
	OUTPUT_VARIABLE _full_git_describe
	OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(NOT (_git_status STREQUAL "0"))
	set(_full_git_describe NOTFOUND)
      endif()
    else()
      message( WARNING "Could not find git executable!")
    endif()
  endif()

  if(NOT (_package_version MATCHES ${_VERSION_REGEX}))
    message( WARNING "Could not extract version from git, falling back on ${_VERSION_FILE}.")
    file(STRINGS ".VERSION" _package_version
      REGEX ${_VERSION_REGEX}
      )
  endif()

  if(NOT _full_git_describe)
    set(_full_git_describe ${_package_version})
  endif()

  # Strip leading "v" character from package version tags so that
  # the version string can be passed to the CMake `project` command
  string(REPLACE "v" "" _package_version "${_package_version}")
  string(REPLACE "V" "" _package_version "${_package_version}")

  if(set_version_VERSION_VARIABLE)
    set(${set_version_VERSION_VARIABLE} ${_package_version} PARENT_SCOPE)
  else()
    set(PROJECT_VERSION ${_package_version} PARENT_SCOPE)
  endif()
  if(set_version_GIT_DESCRIBE_VAR)
    set(${set_version_GIT_DESCRIBE_VAR} ${_full_git_describe} PARENT_SCOPE)
  else()
    set(FULL_GIT_DESCRIBE ${_full_git_describe} PARENT_SCOPE)
  endif()

endfunction()
