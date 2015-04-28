# Top-Level CMakeList.txt

cmake_minimum_required (VERSION 2.6)
project ( Rasterization )

set(EXECUTABLE_NAME "Rasterization")
add_executable( Rasterization skeleton.cpp)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules" ${CMAKE_MODULE_PATH})

set(SDL_LIBRARY_TEMP "C:/Coding/SDL/SDL-1.2.15/lib/x86/SDL.lib;C:/Coding/SDL/SDL-1.2.15/lib/x86/SDLmain.lib")
set(SDL_INCLUDE_DIR "C:/Coding/SDL/SDL-1.2.15/include")

find_package (SDL REQUIRED)

if ( NOT SDL_FOUND )
   message ( FATAL_ERROR "SDL not found!" )
endif ( NOT SDL_FOUND )

if(SDL_FOUND)
	include_directories(
		${SDL_INCLUDE_DIR}
		${PROJECT_SOURCE_DIR}/glm
	)
	#link_libraries(${SDL_LIBRARY})
	target_link_libraries(Rasterization ${SDL_LIBRARY})
endif(SDL_FOUND)

