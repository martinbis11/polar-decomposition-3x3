set (PROJECT_NAME polar-decomp-3x3)

cmake_minimum_required (VERSION 2.8)
project (${PROJECT_NAME})

# Add the files from the project.
file (GLOB exec_files
           "include/*.h"
           "*.cpp"
           )

# Link everything.
add_executable (${PROJECT_NAME} ${exec_files})
