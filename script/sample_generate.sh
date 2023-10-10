#!/bin/bash

# shut up tensorflow!
export TF_CPP_MIN_LOG_LEVEL=3

# change ./models/ accordingly
# samples 10 random molecules
molecule_generation sample ./models/ 10