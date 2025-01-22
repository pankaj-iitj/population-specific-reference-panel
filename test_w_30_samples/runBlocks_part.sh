#!/bin/bash
mkdir example_data && cd example_data
mkdir bed cov_matrix minima vector && cd cov_matrix
mkdir scripts
for i in {1..22}; do mkdir chr${i}; done
cd ../../
