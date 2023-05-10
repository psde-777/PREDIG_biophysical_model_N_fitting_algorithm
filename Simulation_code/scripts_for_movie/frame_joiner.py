import os
from PIL import Image
import math
from tqdm import tqdm

# create output directory if it does not exist
if not os.path.exists("snapshots"):
    os.mkdir("snapshots")

# find all plot files in frames_py and sort them
plot_files = sorted([f for f in os.listdir("frames_py") if f.startswith("plot_")])

# find all frame files in frames_py and sort them
frame_files = sorted([f for f in os.listdir("frames_py") if f.startswith("frame_")])

# get the number of files
num_files = len(plot_files)

# loop over the files and create snapshots
print("Appending frames")
for i in tqdm(range(num_files),unit='file'):
    # get the filenames for this pair of files
    plot_file = os.path.join("frames_py", plot_files[i])
    frame_file = os.path.join("frames_py", frame_files[i])

    # set the output filename with 5 significant digits
    output_file = os.path.join("snapshots", f"snap_{i:05d}.png")

    # open the images and combine them side-by-side
    with Image.open(frame_file) as frame_img, Image.open(plot_file) as plot_img:
        width, height = frame_img.size
        new_width = width + plot_img.size[0]
        new_height = max(height, plot_img.size[1])
        new_img = Image.new("RGB", (new_width, new_height), color="white")
        new_img.paste(frame_img, (0, 0))
        new_img.paste(plot_img, (width, 0))

        # save the combined image
        new_img.save(output_file, dpi=(160, 160))

print("Done!")

