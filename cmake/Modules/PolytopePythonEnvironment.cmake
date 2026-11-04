#[[
# CMake helper for creating Polytope Python virtual environments.
#
# This creates a target named polytope_python_env_build. The environment gets a
# polytope.pth file containing a relative path from the venv site-packages
# directory to the requested pth target.
#
# param[in] env_root: Directory to create the virtual environment.
# param[in] pth_target: Related lib directory. For the build, this is just lib.
#                       For the install, this is lib/python/site-packages.
#
# Named arguments:
#   ALL: If set, target is added to the default build target.
#   TARGET_NAME: Name of the target associated with the env.
#   DEPENDS: List of dependencies.
#   PYTHON_FILES: List of python files to install into the env.
]]
function(polytope_add_python_environment env_root pth_target)
  set(options ALL)
  set(oneValueArgs TARGET_NAME)
  set(multiValueArgs DEPENDS PYTHON_FILES)

  cmake_parse_arguments(arg
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT Python3_EXECUTABLE)
    message(FATAL_ERROR "polytope_add_python_environment requires Python3_EXECUTABLE")
  endif()

  if(NOT POLYTOPE_SITE_PACKAGES_PATH)
    message(FATAL_ERROR "polytope_add_python_environment requires POLYTOPE_SITE_PACKAGES_PATH")
  endif()

  if(NOT arg_TARGET_NAME)
    message(FATAL_ERROR "Must specify a TARGET_NAME")
  endif()
  set(target "${arg_TARGET_NAME}")

  if(NOT IS_ABSOLUTE "${env_root}")
    get_filename_component(env_root "${env_root}" ABSOLUTE
      BASE_DIR "${CMAKE_BINARY_DIR}")
  endif()

  if(NOT IS_ABSOLUTE "${pth_target}")
    get_filename_component(pth_target "${pth_target}" ABSOLUTE
      BASE_DIR "${CMAKE_BINARY_DIR}")
  endif()

  set(env_site_packages "${env_root}/${POLYTOPE_SITE_PACKAGES_PATH}")
  set(stamp "${env_root}/polytope_venv.stamp")
  set(script "${CMAKE_CURRENT_BINARY_DIR}/${target}.cmake")

  file(RELATIVE_PATH pth_entry "${env_site_packages}" "${pth_target}")

  file(WRITE "${script}" "
execute_process(
  COMMAND \"${Python3_EXECUTABLE}\" -m venv --symlinks \"${env_root}\"
  RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR \"Failed to create Python virtual environment: ${env_root}\")
endif()

set(HIBIT_ENABLED ${POLYTOPE_ENABLE_HIBIT2D})
if(NOT ${HIBIT_ENABLED})
  execute_process(
    COMMAND ${env_root}/bin/python3 -m pip install -r ${POLYTOPE_RUNTIME_REQ}
    RESULT_VARIABLE result
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR \"Failed to install runtime requirements into Python virtual environment: ${env_root}\")
  endif()
endif()

file(MAKE_DIRECTORY \"${env_site_packages}\")
file(WRITE \"${env_site_packages}/polytope.pth\" \"${pth_entry}\\n\")
file(WRITE \"${stamp}\" \"ok\\n\")
")

  add_custom_command(
    OUTPUT "${stamp}"
    COMMAND "${CMAKE_COMMAND}" -P "${script}"
    DEPENDS ${arg_DEPENDS}
    VERBATIM
  )

  if(arg_ALL)
    add_custom_target("${target}" ALL
      DEPENDS "${stamp}"
    )
    set(_output_python_files )
    foreach(_pyfile ${arg_PYTHON_FILES})
      cmake_path(GET _pyfile FILENAME _pyfilename)
      set(_cur_pyfile "${pth_target}/${_pyfilename}")
      list(APPEND _output_python_files ${_cur_pyfile})
      add_custom_command(
        TARGET ${target} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${_pyfile}" "${pth_target}"
        DEPENDS ${arg_DEPENDS} ${_pyfile}
      )
    endforeach()
  else()
    add_custom_target("${target}"
      DEPENDS "${stamp}"
    )
    install(FILES ${arg_PYTHON_FILES}
      DESTINATION ${pth_target}
    )
  endif()
endfunction()
