import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import glob
import os
import re
from tqdm import tqdm




# Print a message to be patient with a cow

# Print the initial message
print("""
    ________________________________________________________
   /                                                        \\
  |    Your movie is getting ready. Please be patient...     |
  |    -------------------------------------------------     |
  |         \\   ^__^                                         |
  |          \\  (oo)\\_______                                 |
  |             (__)\\       )\\/\\                             |
  |                 ||----w |                                |
  |                 ||     ||                                |
   \\                                                         /
    --------------------------------------------------------
""")



# Define custom color map
colors = ['#0000FF', '#FFFF00', '#00FF00']
cmap = mcolors.ListedColormap(colors)

# Set the path for the files
path = 'Output/3D/'

# Define marker styles
marker_dict = {0: 'v', 1: 'o', 2: 'o'}

# Create output directory if it does not exist
output_dir = 'frames_py'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


file_list = sorted(glob.glob('Output/3D/visualisation_total_*.txt'), key=lambda x: int(x.split('_')[-1].split('.')[0]))



# Set figure size and resolution
width=1080
height=1080
dpi = 100

i = int(0)

print("Making digestion frames")
# Loop over files and create frames
for file_name in tqdm(file_list, unit='file'):
    # Load data from file
    data = np.loadtxt(file_name, usecols=(0, 1, 2, 3, 4))

    # Extract data columns
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    color_val = data[:, 3]
    timer = data[0, 4]


    fig = plt.figure(figsize=(width/dpi, height/dpi),dpi=dpi)
    ax = fig.add_subplot(projection='3d')

    ax.scatter3D(z, x, y, c=color_val, s=300, cmap=cmap, depthshade=True, marker='H', edgecolor='black')

    ax.set_xlim(-5, 205)
    ax.set_ylim(-5, 9)
    ax.set_zlim(-7, 9)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    # Rotate the plot and set the camera angle
    ax.view_init(elev=10, azim=i*0.25)
    # Set title

    plt.title(f'Time = {timer} hours',fontsize=22)
    plt.tight_layout()
    fig.savefig(f'frames_py/frame_{i+1:05d}.png', dpi=100)
    i=i+1
     # Save frame
    frame_num = int(os.path.splitext(os.path.basename(file_name))[0].split('_')[-1])
    file_path = os.path.join(output_dir, f'frame_{frame_num:05d}.png')
    plt.close(fig)
