#!/bin/bash

workDir="/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib"
inputDir=$workDir"/MC"
outputDir=$workDir"/result"
mkdir -p $outputDir
#source $workDir"/run/pcm.sh"

option1=ooc_ecal_filter # east or west
option2=ooc_sphere_secondary_proton # east or west

inputFile1=$inputDir/${option1}.root
inputFile2=$inputDir/${option2}.root
outputFile=$outputDir/ProtonCalib_${option1}.root
echo "======= inputDir: $inputFile1 $inputFile2"
../build/VLAST_ProtonCalib $inputFile1 $inputFile2 $outputFile
echo "======= outputDir: $outputFile"
