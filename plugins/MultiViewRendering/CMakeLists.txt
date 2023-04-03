cmake_minimum_required(VERSION 3.24)

project(multi_view_rendering)

# Fetch the CGV Framework if it's not obviously present on the system
if (DEFINED ENV{CGV_DIR})
	set(CGV_DIR $ENV{CGV_DIR})
	message(STATUS "Using CGV Framework source found at:")
	message(STATUS "  - ${CGV_DIR}")
	add_subdirectory(${CGV_DIR} ./framework)
else()
	include(FetchContent)
	FetchContent_Declare(
		cgv GIT_REPOSITORY https://github.com/sgumhold/cgv.git
		GIT_TAG 2157a5f2eed33ec8da39d8b6f4b355c98647dc2d
	)
	FetchContent_MakeAvailable(cgv)
endif()

# compile a list of source files for each specific source type the CGV CMake build system knows about
file(GLOB_RECURSE SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cxx")
file(GLOB_RECURSE HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.h")
file(GLOB_RECURSE SHADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.gl*")

# add our target to the CGV CMake build system
cgv_add_target(multi_view_rendering
	TYPE plugin     # if we didn't want a single executable target to be generated, we could add the NO_EXECUTABLE flag
	SOURCES         ${SOURCES}
	HEADERS         ${HEADERS}
	SHADER_SOURCES  ${SHADERS}
	# RESOURCES     if we wanted to include any resources we could include them by listing them here
	DEPENDENCIES
		cgv_viewer glew cgv_utils cgv_type cgv_reflect cgv_data cgv_signal cgv_base cgv_media cgv_gui cgv_render cgv_os
		fltk cgv_gl cgv_g2d cgv_post cgv_proc cgv_app cg_fltk cg_ext crg_light cg_meta cmi_io crg_grid cg_gamepad
		cmf_tt_gl_font cgv_gpgpu

	# we are not using symbols from other plugins directly anywhere in the examples code, but if we did, we could specify them
	# here so the linker would know to look for these symbols in the import libs / shared objects of the listed plugins
	# LINKTIME_PLUGIN_DEPENDENCIES crg_bg

	ADDITIONAL_CMDLINE_ARGS
		"config:\"${CMAKE_CURRENT_LIST_DIR}/config.def\""

	# we are not overriding the default preprocesser define that indicates we are building the shared library and symbols
	# should be exported (in this case, the default is EXAMPLES_EXPORTS), but if we wanted to, we could
	# OVERRIDE_SHARED_EXPORT_DEFINE I_AM_EXPORTING

	# we are not overriding the default preprocesser define that indicates we are building the static version for use by
	# single executables (in this case, the default is EXAMPLES_FORCE_STATIC), but if we wanted to, we could
	# OVERRIDE_FORCE_STATIC_DEFINE  I_AM_STATIC

	# we don't need any additional defines private or public to the build of the 'examples' plugin, but if we did,
	# we could declare them here
	# ADDITIONAL_PRIVATE_DEFINES  SOME_INTERNAL_PREPROCESSOR_MACRO
	# ADDITIONAL_PUBLIC_DEFINES   SOME_PUBLICALLY_VISIBLE_PREPROCESSER_MACRO
)

# Visual Studio fluff
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT multi_view_rendering)

# FIXME: CMake installs are currently still unsupported and need updating for the new CMake build system
install(TARGETS multi_view_rendering EXPORT cgv_plugins DESTINATION ${CGV_BIN_DEST})
install(FILES config.def DESTINATION ${CGV_PLUGINS_INCLUDE_DEST}/multi_view_rendering)
install(DIRECTORY . DESTINATION ${CGV_PLUGINS_INCLUDE_DEST}/multi_view_rendering FILES_MATCHING PATTERN "*.gl*")
