#!/bin/sh
gcc -c adj_cellwise_parcel.cpp -std=c++11
gfortran -c parcel_GC.JAS99_libcloud_blk_1m.for 
gfortran -o parcel_blk_1m adj_cellwise_parcel.o parcel_GC.JAS99_libcloud_blk_1m.o -lstdc++
