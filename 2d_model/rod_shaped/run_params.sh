#!/bin/bash

# Define parameters as variables
box=20
cutoff=4
Nparticles=1
T=800
Dt=0.01
rod_length=2
rod_radius=0.2
khardcore=0
kspring=2
rate=0.3
lo=10
kalign=12

# Create folder name based on parameters
folder_name="Dt_${Dt}_Nparticles_${Nparticles}_T_${T}_box_${box}_cutoff_${cutoff}_kalign_${kalign}_khardcore_${khardcore}_kspring_${kspring}_lo_${lo}_rate_${rate}_rod_length_${rod_length}_rod_radius_${rod_radius}"

# Create the directory
mkdir -p "$folder_name"

# Create and populate params.dat file
echo "box $box" >> "$folder_name/params.dat"
echo "cutoff $cutoff" >> "$folder_name/params.dat"
echo "Nparticles $Nparticles" >> "$folder_name/params.dat"
echo "T $T" >> "$folder_name/params.dat"
echo "Dt $Dt" >> "$folder_name/params.dat"
echo "kspring $kspring" >> "$folder_name/params.dat"
echo "rate $rate" >> "$folder_name/params.dat"
echo "lo $lo" >> "$folder_name/params.dat"
echo "rod_length $rod_length" >> "$folder_name/params.dat"
echo "khardcore $khardcore" >> "$folder_name/params.dat"
echo "rod_radius $rod_radius" >> "$folder_name/params.dat"
echo "kalign $kalign" >> "$folder_name/params.dat"


# Copy the executable to the directory
cp rods_twitch_sim_torque "$folder_name/"

# Change to the directory and execute the program
cd "$folder_name" || exit
./rods_twitch_sim_torque params.dat

rm rods_twitch_sim_torque 

echo "Simulation ${folder_name} completed."