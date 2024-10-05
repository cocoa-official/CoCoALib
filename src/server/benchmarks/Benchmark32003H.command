#!/bin/bash

# this file should be saved into you CoCoALib directory

BENCHMARK_PATH=`dirname $0`

INPUT_PATH=$BENCHMARK_PATH/inputs
SERVER="$BENCHMARK_PATH/../CoCoAServer -t"

date "+%n=====  DATE: %Y-%m-%d TIME: %H:%M  =====  char 32003, homog  ====="

$SERVER <  $INPUT_PATH/t51_h.cocoa5;
$SERVER <  $INPUT_PATH/dense755_h.cocoa5;
$SERVER <  $INPUT_PATH/kin1_h.cocoa5;
$SERVER <  $INPUT_PATH/twomat3_h.cocoa5;
$SERVER <  $INPUT_PATH/intersection-c4-c5_h.cocoa5;
$SERVER <  $INPUT_PATH/c7_h.cocoa5;
$SERVER <  $INPUT_PATH/6x7-4_h.cocoa5;
$SERVER <  $INPUT_PATH/gaukwa4_h.cocoa5;
$SERVER <  $INPUT_PATH/wang_hDL.cocoa5;
$SERVER <  $INPUT_PATH/homog_gonnet_h.cocoa5;
$SERVER <  $INPUT_PATH/t51_hDL.cocoa5;
$SERVER <  $INPUT_PATH/7x8-2_h.cocoa5;
$SERVER <  $INPUT_PATH/mora1_h.cocoa5;
$SERVER <  $INPUT_PATH/kin1_hDL.cocoa5;
$SERVER <  $INPUT_PATH/hairer2_h.cocoa5;
#$SERVER <  $INPUT_PATH/colon-c4-c5_h.cocoa5;


## 4x4.cocoa5
## prova12.cocoa5
## virasoro_h.cocoa5
