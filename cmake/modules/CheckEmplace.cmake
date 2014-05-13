#
# Module that checks whether std::map::emplace is available.
#
# Sets the follwing variable:
#
# HAVE_EMPLACE                     True if map::emplace is available

include(CMakePushCheckState)
cmake_push_check_state()

# check whether map::emplace exists
check_cxx_source_compiles("
  #include <map>

  int main(int argc, char * argv[])
  {
    std::map<int, int> myMap;
    myMap.emplace(0, 0);
    return 0;
  }
" HAVE_EMPLACE
)

cmake_pop_check_state()
