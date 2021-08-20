#!/bin/bash

../bin/ExpansionHunter \
  --reads input/variants.bam \
  --reference input/reference.fa \
  --variant-catalog input/variants.json \
  --output-prefix output/repeats
  
