#!/bin/bash

workDir="/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib"
inputDir=$workDir"/MC"
outputDir=$workDir"/result"
mkdir -p $outputDir
#source $workDir"/run/pcm.sh"

option=ooc_ecal_filter_reverse_spectrum # east or west
inputFile=$inputDir/${option}.root
outputFile=$outputDir/ProtonCalib_${option}.root
echo "======= inputDir: $inputFile"
../build/VLAST_ProtonCalib $inputFile $outputFile
echo "======= outputDir: $outputFile"
