#!/bin/bash

# do mcscanx
mcscanx=$1
# mcscanx=/home/yutachen/public/Yutangchen/mcscanx/MCScanX-master/MCScanX

if [ -d mymcscanx ]
then
    rm -r mymcscanx
    mkdir mymcscanx
else
    mkdir mymcscanx
fi

cp mcscanx.gff mymcscanx/
cp mcscanx.blast mymcscanx/

$mcscanx mymcscanx/mcscanx
