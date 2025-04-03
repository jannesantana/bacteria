#!/bin/bash

# Define parameters as variables
Dt=0.01
Nparticles=200
T=1000
box=10
cutoff=4
epsilon=1
kspring=4
lo=0.8
rate=20
sigma=0.15

# Create folder name based on parameters
folder_name="Dt_${Dt}_Nparticles_${Nparticles}_T_${T}_box_${box}_cutoff_${cutoff}_epsilon_${epsilon}_kspring_${kspring}_lo_${lo}_rate_${rate}_sigma_${sigma}"

# Create the directory
mkdir -p "$folder_name"

# Create and populate params.dat file
echo "box $box" > "$folder_name/params.dat"
echo "cutoff $cutoff" >> "$folder_name/params.dat"
echo "Nparticles $Nparticles" >> "$folder_name/params.dat"
echo "T $T" >> "$folder_name/params.dat"
echo "Dt $Dt" >> "$folder_name/params.dat"
echo "kspring $kspring" >> "$folder_name/params.dat"
echo "rate $rate" >> "$folder_name/params.dat"
echo "lo $lo" >> "$folder_name/params.dat"
echo "epsilon $epsilon" >> "$folder_name/params.dat"
echo "sigma $sigma" >> "$folder_name/params.dat"

# Copy the executable to the directory
cp twitching_sim_2 "$folder_name/"

# Change to the directory and execute the program
cd "$folder_name" || exit
./twitching_sim_2 params.dat

rm twitching_sim_2 

echo "Simulation ${folder_name} completed."