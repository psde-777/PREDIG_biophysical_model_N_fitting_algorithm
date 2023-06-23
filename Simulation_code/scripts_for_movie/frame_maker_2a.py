import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


# create output directory if it does not exist
if not os.path.exists('frames_py'):
    os.makedirs('frames_py')

# find all files in Output/3D/ that match the pattern conversion_1_*.txt and sort them
filepaths = sorted(glob.glob('Output/3D/conversion_1_*.txt'), key=lambda x: int(re.findall(r'\d+(?=\.txt)', x)[0]))


# Set figure size and resolution
width=840
height=1080
dpi = 100

print("Making saccharification plots")
# loop over all files and plot the data
with tqdm(total=len(filepaths), unit='file') as pbar:
    for i, filepath in enumerate(filepaths):
        # load data from file into a numpy array
        data = np.loadtxt(filepath)

        # extract the x, y1, and y2 data columns
        x = data[:, 0]
        y1 = data[:, 1]
        y2 = data[:, 3]

        # create two subplots, one for glucose and one for xylose
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(width/dpi, height/dpi), dpi=dpi)

        # plot the glucose data on the first subplot
        axs[0].plot(x, y1, linewidth=4, label='simulated', color='green')
        axs[0].set_xlabel('time (hours)', fontsize=18)
        axs[0].set_ylabel('cellulose hydrolysed %', fontsize=18)
        axs[0].legend(loc='lower right', fontsize=18)
        axs[0].set_xlim(0, )
        axs[0].set_ylim(0, 100)
        axs[0].tick_params(axis='both', which='major', labelsize=16)

        # plot the xylose data on the second subplot
        axs[1].plot(x, y2, linewidth=4, label='simulated', color='orange')
        axs[1].set_xlabel('time (hours)', fontsize=18)
        axs[1].set_ylabel('xylose hydrolysed %', fontsize=18)
        axs[1].legend(loc='lower right', fontsize=18)
        axs[1].set_xlim(0, )
        axs[1].set_ylim(0, 100)
        axs[1].tick_params(axis='both', which='major', labelsize=16)

        # save the plot to the output directory with a name matching the input file
        output_filename = os.path.join('frames_py', f'plot_{i+1:05d}.png')
        plt.savefig(output_filename)

        # close the figure to free up memory
        plt.close(fig)

        # update progress bar
        pbar.update(1)

