#!/usr/bin/env bash

## run_myLLannotator.sh by Rohan Maddamsetti.
## install myLLannotator from here: https://github.com/alyssa-lee/myLLannotator

myllannotator myLLannotator-inputs/valid_categories.txt myLLannotator-inputs/system_prompt.txt myLLannotator-inputs/per_sample_prompt.txt myLLannotator-inputs/gbk-annotation-table.csv ../results/labeled-gbk-annotation-table.csv
