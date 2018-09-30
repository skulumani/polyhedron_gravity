~~~
brew install cgal cmake libomp
~~~

CMake > 3.12

Eigen > 3.3.5

~~~
git submodule update --init --recursive 
~~~

CMake don't search for matlab in libigl

Might need to edit 

~~~
extern/libigl/shared/cmake/libigl.cmake
~~~
