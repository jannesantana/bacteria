import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.patches import Rectangle,Circle


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 25
import os 

import matplotlib.transforms as transforms

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.patches import Rectangle, Circle
import matplotlib.transforms as transforms
import os
import subprocess
from tqdm import tqdm

def movie_maker_capsules(folder_name, X, Y, Theta, T0, Tf, BOX_SIZE, step,
                         ROD_LENGTH=0.6, ROD_WIDTH=0.1, 
                         make_ffmpeg_movie=False, del_pngs=False):
    
    os.makedirs(folder_name, exist_ok=True)  # Ensure folder exists

    fig, ax = plt.subplots(figsize=(5, 5))
    radius = ROD_WIDTH / 2
    rect_length = ROD_LENGTH - 2 * radius
    particles = []  # Store (rectangle, cap1, cap2) for each particle

    # Setup plot
    ax.set_xlim(0, BOX_SIZE)
    ax.set_ylim(0, BOX_SIZE)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_aspect('equal', adjustable='box')

    # Create all patches initially
    Nparticles = X.shape[1]
    for _ in range(Nparticles):
        rect = Rectangle((-rect_length/2, -ROD_WIDTH/2), rect_length, ROD_WIDTH,
                         facecolor='red', edgecolor='none')
        cap1 = Circle((0, 0), radius=radius, facecolor='red', edgecolor='none')
        cap2 = Circle((0, 0), radius=radius, facecolor='red', edgecolor='none')
        ax.add_patch(rect)
        ax.add_patch(cap1)
        ax.add_patch(cap2)
        particles.append((rect, cap1, cap2))

    plt.tight_layout()

    # Loop over frames
    count = 0
    for time in tqdm(range(T0, Tf,step), desc="Generating frames"):
        Xnow, Ynow, Thetanow = X[time], Y[time], Theta[time]
        ux, uy = np.cos(Thetanow), np.sin(Thetanow)

        for i, (rect, cap1, cap2) in enumerate(particles):
            posx, posy, theta = Xnow[i], Ynow[i], Thetanow[i]
            # Update rectangle transform
            t = transforms.Affine2D().rotate_around(0, 0, theta).translate(posx, posy) + ax.transData
            rect.set_transform(t)

            # Update cap positions
            cap1_center = (posx - (ROD_LENGTH/2 - radius)*ux[i], posy - (ROD_LENGTH/2 - radius)*uy[i])
            cap2_center = (posx + (ROD_LENGTH/2 - radius)*ux[i], posy + (ROD_LENGTH/2 - radius)*uy[i])
            cap1.center = cap1_center
            cap2.center = cap2_center

        ax.set_title(f'Time {time}', fontsize=12, loc='left')
        fig.savefig(f"{folder_name}/positions_{count:03d}.png", dpi=200, bbox_inches='tight')
        count +=1

    plt.close(fig)

    # Make movie if requested
    if make_ffmpeg_movie:
        subprocess.run([
            "ffmpeg", "-y", "-r", "45", "-start_number", f"{0:03d}",
            "-i", f"{folder_name}/positions_%03d.png",
            "-vf", "crop=trunc(iw/2)*2:trunc(ih/2)*2",
            "-pix_fmt", "yuv420p", f"{folder_name}/movie.mp4"
        ], check=True)

    # Optionally delete pngs
    if del_pngs:
        for time in range(count):
            os.remove(f"{folder_name}/positions_{time:03d}.png")

box=30
cutoff=5
Nparticles=300
T=12000
Dt=0.005
rod_length=2
rod_radius=0.2
khardcore=11
kspring=5
rate=0.5
lo=2
kalign=4
noise_pos=0.01
noise_thetab=0.01

folder_name=f"Dt_{Dt}_Nparticles_{Nparticles}_T_{T}_box_{box}_cutoff_{cutoff}_kalign_{kalign}_khardcore_{khardcore}_kspring_{kspring}_lo_{lo}_noise_pos_{noise_pos}_noise_thetab_{noise_thetab}_rate_{rate}_rod_length_{rod_length}_rod_radius_{rod_radius}"
positions_file = f"{folder_name}/particle_positions_{folder_name}.dat"
squared_disp_file = f"{folder_name}/squared_disp_{folder_name}.dat"
X,Y,Theta_b,Theta_p = np.loadtxt(positions_file, unpack=True)
Steps = len(X)//Nparticles
X = X.reshape((Steps, Nparticles))
Y = Y.reshape((Steps, Nparticles))
Theta_b = Theta_b.reshape((Steps,Nparticles))

T0 = 0
Tf = Steps
movie_maker_capsules(folder_name,X,Y,Theta_b,T0,Tf,box,2,ROD_LENGTH=rod_length,ROD_WIDTH=2*rod_radius,make_ffmpeg_movie=True,del_pngs=True)
# ffmpeg -r 5 -start_number 000 -i 'positions_%03d.png' -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2"  -pix_fmt yuv420p movie.mp4
