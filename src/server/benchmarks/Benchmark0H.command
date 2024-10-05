#!/bin/bash

# this file should be saved into your CoCoALib directory

BENCHMARK_PATH=`dirname $0`

INPUT_PATH=$BENCHMARK_PATH/inputs
SERVER="$BENCHMARK_PATH/../CoCoAServer -t"

date "+%n=====  DATE: %Y-%m-%d TIME: %H:%M  =====  char 0, homog  ====="

$SERVER <  $INPUT_PATH/t51_hQ.cocoa5;
$SERVER <  $INPUT_PATH/alex3_hQ.cocoa5;
$SERVER <  $INPUT_PATH/schwartz7_hQDL.cocoa5;
$SERVER <  $INPUT_PATH/katsura7_hQ.cocoa5;
#$SERVER <  $INPUT_PATH/katsura7_hQDL.cocoa5;

$SERVER <  $INPUT_PATH/6x7-4_hQ.cocoa5;
$SERVER <  $INPUT_PATH/twomat3_hQ.cocoa5;
$SERVER <  $INPUT_PATH/c7_hQ.cocoa5;

$SERVER <  $INPUT_PATH/homog_gonnet_hQ.cocoa5;
$SERVER <  $INPUT_PATH/virasoro_hQ.cocoa5;

$SERVER <  $INPUT_PATH/kin1_hQ.cocoa5;

$SERVER <  $INPUT_PATH/hairer2_hQ.cocoa5;
$SERVER <  $INPUT_PATH/mora9_hQ.cocoa5;

## 4x4.cocoa5
## prova12.cocoa5
## virasoro_h.cocoa5
