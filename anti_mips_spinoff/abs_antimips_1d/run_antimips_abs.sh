#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <parameter_file.dat>"
    exit 1
fi

# Check if the parameter file exists
if [ ! -f "$1" ]; then
    echo "Error: Parameter file '$1' not found."
    exit 1
fi

# Read parameters from the .dat file and store them in an associative array
declare -A parameters
while read -r label value; do
    parameters["$label"]=$value
done < "$1"

# execute with parameters
./antimips_1d \
    "${parameters[N]}" \
    "${parameters[Lx]}" \
    "${parameters[R]}" \
    "${parameters[eta]}" \
    "${parameters[T]}" \
    "${parameters[Dt]}" \
    "${parameters[Vo]}" \
    "${parameters[gamma]}"
exit 0
