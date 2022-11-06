# clothSim3D
## Build
Make sure you have [CMake](https://cmake.org/)
```
git clone git@github.com:Aviii06/clothSim3D.git
cd clothSim3D
git submodule update --init --recursive
mkdir build && cd build
cmake ..
cmake --build . -j8
./clothSim 300
```

Then run the output file through [polybin](https://github.com/Aviii06/polybin)
