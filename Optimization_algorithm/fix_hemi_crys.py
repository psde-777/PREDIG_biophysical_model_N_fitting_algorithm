import os

# Define a function to check and update the numbers in a file
def update_file_numbers(filename):
    try:
        with open(filename, 'r') as file:
            data = file.read().strip().split()
            if len(data) != 18:
                print(f"Skipping {filename}: Incorrect number of values.")
                return

            # Convert the 14th and 15th numbers to floats
            cellu_cry, hemi_cry = float(data[13]), float(data[14])

            # Check if hemi_cry is larger than cellu_cry
            if hemi_cry > cellu_cry:
                # Update the value of hemi_cry to be the same as cellu_cry
                data[14] = str(cellu_cry)

                # Write the updated data back to the file
                with open(filename, 'w') as updated_file:
                    updated_file.write('\t'.join(data))

                print(f"Updated {filename}. Hemicellulose crystallinity cannot exceed cellulose crystallinity. Changed value to {cellu_cry}")
            else:
                print("All good. Hemicellulose crystallinity is lower than cellulose crystallinity.")
    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")

# Read the keywords from the keywords.txt file
try:
    with open('keywords.txt', 'r') as keywords_file:
        keywords = keywords_file.read().strip().split('\n')
except FileNotFoundError:
    print("Error: keywords.txt file not found.")
    exit(1)

# Define the directory where your files are located
directory = 'BEST_FIT/best_Run/Params/'

# Loop through the keywords and process the corresponding files
for keyword in keywords:
    filename = os.path.join(directory, f"initial_configuration_parameters_{keyword}.txt")
    if os.path.exists(filename):
        update_file_numbers(filename)
    else:
        print(f"File not found: {filename}")

print("Processing complete.")
