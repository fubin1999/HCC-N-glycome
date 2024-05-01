#!/bin/bash

for i in {1..8}
do
  glyhunter run \
    data/mass_list_combined/plate${i}.xlsx \
    -o results/data/glyhunter_results/plate${i}_results
done