#!/usr/bin/env bash

## run_myLLannotator.sh by Rohan Maddamsetti.
## install myLLannotator from here: https://github.com/alyssa-lee/myLLannotator

myllannotator ../results/myLLannotator-inputs/valid_categories.txt ../results/myLLannotator-inputs/system_prompt.txt ../results/myLLannotator-inputs/per_sample_prompt.txt ../results/myLLannotator-inputs/gbk-annotation-table.csv ../results/labeled-gbk-annotation-table.csv
