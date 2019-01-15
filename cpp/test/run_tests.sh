#/usr/bin/env sh

c++ -std=c++11 -O3 -fPIC -shared \
    `root-config --cflags` \
    -o libbayesian_blocks.so ../bayesian_blocks.cc \
    `root-config --libs`

c++ -std=c++11 -O3 -L. -I.. \
    `root-config --cflags` \
    -o test test.cc \
    -lbayesian_blocks `root-config --libs`

./test
