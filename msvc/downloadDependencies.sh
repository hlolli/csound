#!/bin/sh

mkdir -p cache 
mkdir -p deps 

cd cache 

wget -nc http://www.mega-nerd.com/libsndfile/files/libsndfile-1.0.27-w64.zip 


cd ../deps

unzip ../cache/libsndfile-1.0.27-w64.zip

