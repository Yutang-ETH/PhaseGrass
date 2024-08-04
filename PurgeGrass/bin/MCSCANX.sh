#!/bin/bash

# do mcscanx

mcscanx=/home/yutachen/public/Yutangchen/mcscanx/MCScanX-master/MCScanX

if [ -d mymcscanx ]
then
    rm -r mymcscanx
    mkdir mymcscanx
else
    mkdir mymcscanx
fi

mv mcscanx.gff mymcscanx/
mv mcscanx.blast mymcscanx/

$mcscanx mymcscanx/mcscanx
