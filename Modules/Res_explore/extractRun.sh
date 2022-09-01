#!/bin/bash

# Results extraction script
#
# @zgr2788

sampleName=$1
dir="runSum_$sampleName"
mkdir -p $dir

# Pseudobulks
mkdir -p "${dir}/Pseudobulks"
mode=$2
cellCount=$3
nSamples=$4
propVar=$5
sampleCT=$6

files=$(ls Input/Psuedobulks | grep $sampleName)
for f in $files
do
  fname=$(echo "$f" | sed "s@.rds@@g")
  cp "Input/Psuedobulks/$f" "${dir}/Pseudobulks/${fname}_${mode}_${cellCount}_${nSamples}_${propVar}_${sampleCT}.rds"
done


# Normtables
mkdir -p "${dir}/NormalizedTables"
scaleTMethod=$7
scaleCMethod=$8
transformMethod=$9

files=$(ls Input/Normalized_tables | grep $sampleName)
for f in $files
do
  fname=$(echo "$f" | sed "s@.rds@@g")

  if [[ $fname =~ "pbulks" ]]
  then
    cp "Input/Normalized_tables/$f" "${dir}/NormalizedTables/${fname}_${mode}_${cellCount}_${nSamples}_${propVar}_${sampleCT}_scale=${scaleTMethod}_transform=${transformMethod}.rds"
  else
    cp "Input/Normalized_tables/$f" "${dir}/NormalizedTables/${fname}_scale=${scaleCMethod}_transform=${transformMethod}.rds"
  fi
done


# Metrics
mkdir -p "${dir}/Metrics"

files=$(ls Metrics | grep $sampleName)
for f in $files
do
  fname=$(echo "$f" | sed "s@.rds@@g")
  cp  "Metrics/$f" "${dir}/Metrics/${fname}_${mode}_${cellCount}_${nSamples}_${propVar}_${sampleCT}_scale=${scaleTMethod}_transform=${transformMethod}.rds"
done


# Plots
mkdir -p "${dir}/Plots"

files=$(ls Plots | grep $sampleName)
for f in $files
do
  fname=$(echo "$f" | sed "s@.png@@g")
  cp  "Plots/$f" "${dir}/Plots/${fname}_${mode}_${cellCount}_${nSamples}_${propVar}_${sampleCT}_scale=${scaleTMethod}_transform=${transformMethod}.png"
done


# Raw outputs
mkdir -p "${dir}/Outputs"

files=$(ls Output | grep $sampleName)
for f in $files
do
  fname=$(echo "$f" | sed "s@.rds@@g")
  cp  "Output/$f" "${dir}/Outputs/${fname}_${mode}_${cellCount}_${nSamples}_${propVar}_${sampleCT}_scale=${scaleTMethod}_transform=${transformMethod}.rds"
done


# Benchmarks
mkdir -p "${dir}/CompBenchmarks"

files=$(ls Benchmarks | grep $sampleName)
for f in $files
do
  fname=$(echo "$f" | sed "s@.txt@@g")
  cp  "Benchmarks/$f" "${dir}/CompBenchmarks/${fname}_${mode}_${cellCount}_${nSamples}_${propVar}_${sampleCT}_scale=${scaleTMethod}_transform=${transformMethod}.txt"
done


# Get the config as well
cp config.yaml ${dir}/
