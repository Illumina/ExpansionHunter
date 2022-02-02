#!/usr/bin/env bash

scriptDir=$(cd $(dirname $0) && pwd -P)
cd $scriptDir

../bin/ExpansionHunter \
  --reads input/variants.bam \
  --reference input/reference.fa \
  --variant-catalog input/variants.json \
  --output-prefix output/repeats
