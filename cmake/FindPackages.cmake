# In this CMake file, we will find all required packages


# Find the Boost package and some required components
find_package(Boost REQUIRED COMPONENTS system thread program_options)

# Find Eigen3
find_package(Eigen3 3.3.4 REQUIRED)

# Find hf
find_package(hf 3.0.0 REQUIRED)

# Find the libint integral wrapper (also includes support for integral transformations)
find_package(libwint 3.0.0 REQUIRED)

