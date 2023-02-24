import os
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio.v2 as imageio

# Set the path for the files
path = 'Output/3D/'

# Get the list of files for the "total" set in the proper order
total_files = sorted([f for f in os.listdir(path) if re.match(r'visualisation_total_1_\d+\.txt', f)], key=lambda x: int(x.split('.')[0].split('_')[-1]))


# Initialize lists to store data for each frame
x_vals = []
y_vals = []
z_vals = []
color_vals = []

# Loop through each file and extract data for the "total" set
for total in total_files:
    # Extract x, y, z, and color data
    with open(path + total, 'r') as f:
        lines = f.readlines()
        x_vals.append([float(line.split()[0]) for line in lines])
        y_vals.append([float(line.split()[1]) for line in lines])
        z_vals.append([float(line.split()[2]) for line in lines])
        color_vals.append([float(line.split()[3]) for line in lines])


# Create the "frames" directory if it does not exist
frames_folder = 'frames_py'
if not os.path.exists(frames_folder):
    os.makedirs(frames_folder)

print ("Python is creating the frames for your animation. Please be patient.")

# Create the frames
for i in range(len(total_files)):
    fig = plt.figure(figsize=(12, 16))
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect([1,1,1])  # This line removes the borders
    ax.set_axis_off()  # This line removes the box outline
    ax.scatter(x_vals[i], y_vals[i], z_vals[i], c=color_vals[i], s=100, cmap='brg', marker='o')
    ax.set_xlim(-5, 9)
    ax.set_ylim(-5, 8)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    #plt.title(f'frame {i+1}')
    plt.tight_layout()
    fig.savefig(f'frames_py/frame_{i+1}.png', dpi=100)
    plt.close(fig)

print("Your animation is almost ready.")

# Create the GIF
frames = []
for i in range(len(total_files)):
    filename = f'frames_py/frame_{i+1}.png'
    frames.append(imageio.imread(filename))
imageio.mimsave('digest_animation.gif', frames, fps=30)

# remove frames folder
for file_name in os.listdir(frames_folder):
    file_path = os.path.join(frames_folder, file_name)
    os.remove(file_path)
os.rmdir(frames_folder)
