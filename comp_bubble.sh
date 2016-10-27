#!/bin/sh
gcc -c adj_cellwise_parcel.cpp -std=c++11
gfortran -c bubble.2d.Gr.Cl.JAS.1999_blk_1m.for 
gfortran -o bubble_blk_1m adj_cellwise_parcel.o bubble.2d.Gr.Cl.JAS.1999_blk_1m.o -lstdc++
