#!/bin/bash
workDir="/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib"
inputDir=$workDir"/MC"
outputDir=$workDir"/result"
mkdir -p $outputDir

option=east # east or west
inputFile=$inputDir/${option}.root
outputFile=$outputDir/ProtonCalib_${option}.root
echo "======= inputDir: $inputFile"
../build/VLAST_ProtonCalib $inputFile $outputFile
echo "======= outputDir: $outputFile"
