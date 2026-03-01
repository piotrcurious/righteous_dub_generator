cmake_minimum_required(VERSION 3.16)
project(harmonia VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Find packages ──────────────────────────────────────────────────────────────
find_package(PkgConfig REQUIRED)

# FLTK
find_package(FLTK REQUIRED)
if(NOT FLTK_FOUND)
  message(FATAL_ERROR "FLTK not found. Install: apt install libfltk1.3-dev")
endif()

# ALSA
pkg_check_modules(ALSA REQUIRED alsa)

# OpenGL
find_package(OpenGL REQUIRED)
find_package(GLU REQUIRED)

# ── Sources ────────────────────────────────────────────────────────────────────
set(SOURCES
    src/main.cpp
    src/audio_engine.cpp
    src/tonal_space_widget.cpp
    src/theory_bridge.cpp
)

add_executable(harmonia ${SOURCES})

target_include_directories(harmonia PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${FLTK_INCLUDE_DIRS}
    ${ALSA_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIR}
)

target_link_libraries(harmonia PRIVATE
    ${FLTK_LIBRARIES}
    ${FLTK_GL_LIBRARY}
    ${ALSA_LIBRARIES}
    ${OPENGL_LIBRARIES}
    GLU
    pthread
)

target_compile_options(harmonia PRIVATE
    -Wall -Wextra
    -O2
    ${ALSA_CFLAGS_OTHER}
)

# ── Install ────────────────────────────────────────────────────────────────────
install(TARGETS harmonia DESTINATION bin)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/theory/
        DESTINATION share/harmonia/theory)

# ── Post-build: copy theory dir next to binary ─────────────────────────────────
add_custom_command(TARGET harmonia POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/theory
        ${CMAKE_CURRENT_BINARY_DIR}/theory
    COMMENT "Copying theory scripts alongside binary"
)
