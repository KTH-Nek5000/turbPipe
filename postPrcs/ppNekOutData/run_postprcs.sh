#!/bin/bash 

./makenek clean
echo "turbPipe_post" >  SESSION.NAME
echo `pwd`'/' >>  SESSION.NAME
./makenek turbPipe_post
echo "**** Running nek5000 in post-processing mode"
./nek5000 >post2d.log

