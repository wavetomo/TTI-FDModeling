# !/usr/bin/sh
###
# @Author: Pingmin Cheung
# @Date: 2025-08-01 11:13:58
 # @LastEditTime: 2026-06-16 16:04:06
 # @LastEditors: Pingmin Cheung
# Copyright (c) 2025 by Wavetomo, All Rights Reserved.
###
path=$(pwd)
filepath=$path/test_bp02
exe=$path/../bin

rm -rf $filepath
mkdir -p $filepath
cd $filepath

nz=1126 nx=1251 nt=8001
dx=10 dz=10 dt=0.001

stype=1 fpeak=10 sdepth=5 sx=625

gdel=5 offsetmin=0 dtr=1 offsetmax=2000
seistype=0 mt=2

freesurface=0

$exe/cutbin2d <../file/vp_bp_1126_7872_10m.bin \
    n1=1126 n2=7872 fx1=0 fx2=6621 nx1=$nz nx2=$nx \
    >vp.bin
$exe/cutbin2d <../file/delta_bp_1126_7872_10m.bin \
    n1=1126 n2=7872 fx1=0 fx2=6621 nx1=$nz nx2=$nx \
    >delta.bin
$exe/cutbin2d <../file/epsilon_bp_1126_7872_10m.bin \
    n1=1126 n2=7872 fx1=0 fx2=6621 nx1=$nz nx2=$nx \
    >epsilon.bin
$exe/cutbin2d <../file/theta_bp_1126_7872_10m.bin \
    n1=1126 n2=7872 fx1=0 fx2=6621 nx1=$nz nx2=$nx \
    >theta.bin

$exe/sumksrc2d sdepth=$sdepth sx=$sx fpeak=$fpeak A=0.001 \
    stype=1 dt=$dt nt=$nt >srcfile.su
$exe/sumkrcev2d sx=$sx gdel=$gdel \
    offsetmin=$offsetmin dtr=$dtr offsetmax=$offsetmax \
    nx=$nx >recepos.su

$exe/sufdmod2dTTIelastic_tensor_rotation \
    vpfile=./vp.bin \
    epsilonfile=./epsilon.bin \
    deltafile=./delta.bin \
    thetafile=./theta.bin scalar=-1 \
    nz=$nz nx=$nx nt=$nt \
    dx=$dx dz=$dz dt=$dt \
    stype=$stype fpeak=$fpeak sourcefile=srcfile.su \
    seistype=$seistype mt=$mt recposfile=recepos.su \
    savew=0 wsave=2 snap_sample=1000 \
    mode=1 freesurface=$freesurface \
    shotID=$i scalar=-1
cd $path