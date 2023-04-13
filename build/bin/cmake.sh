rm ../../*/._*.cpp
rm ../../src/Lineshapes/._*.cpp
rm ../../src/HyperPlot/._*.cpp
#cmake .. -DCMAKE_CXX_STANDARD=17 -DUSE_SIMD=AVX2d
cmake .. -DCMAKE_CXX_STANDARD=17 -DUSE_SIMD=0
#cmake .. -DUSE_SIMD=0
cd ..
make 
cd bin/
