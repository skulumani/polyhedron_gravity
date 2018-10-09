~~~
brew install cgal cmake libomp
~~~

## Requirements

1. CMake > 3.12

2. Eigen > 3.3.5

3. CGAL > 4.12-1

4. Boost > 1.68

~~~
git submodule update --init --recursive 
~~~

CMake doesn't search for matlab in libigl

Might need to edit 

~~~
extern/libigl/shared/cmake/libigl.cmake
~~~

## TODO

* Remove dependency on libigl
* Bundle CGAL requirements with app (similar to how igl works)
* Remove boost dependency (probably easy)
* Analytical verification of potential https://arxiv.org/pdf/1206.3857.pdf
