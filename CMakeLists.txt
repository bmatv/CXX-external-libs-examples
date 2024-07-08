cmake_minimum_required(VERSION 3.25)
project(netcdf4Examples VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(
"$<$<COMPILE_LANGUAGE:CXX>:-Wall;-pedantic>;-O3;-march=native;-mavx2;-ffast-math"
)

#find_package(netCDF REQUIRED)
find_package(MPI REQUIRED)

add_executable(netcdf4Write netcdf4Write.cpp)
target_link_libraries(netcdf4Write PRIVATE netcdf) # this way -lnetcdf flag will be added which is enough for legacy netCDF

add_executable(mpiExample mpiExample.cpp)
target_link_libraries(mpiExample PRIVATE MPI::MPI_CXX)

add_executable(netcdf4Read netcdf4Read.cpp)
target_link_libraries(netcdf4Read PRIVATE netcdf_c++4) # this way -lnetcdf_c++4 flag will be added
# but then why would one need to find netCDF? to expose netCDF dir?

add_executable(testVec vec_test.cpp)