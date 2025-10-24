from PIL import Image
import sys
import os

############################################################################################################
# Figure S2

# Retrieve directory paths
script_dir = os.path.dirname(os.path.abspath(__file__))
figures_dir = os.path.join(script_dir, '../figures')

# Load the images
image1path = os.path.join(figures_dir, 'OH2_Figure_S2A.png')
image2path = os.path.join(figures_dir, 'OH2_Figure_S2B.png')

try:
    with Image.open(image1path) as image1, Image.open(image2path) as image2:

        # Get dimensions of the images
        width1, height1 = image1.size
        width2, height2 = image2.size

        # Calculate the combined width and height
        combined_width = max(width1, width2)
        combined_height = height1 + height2

        # Create a new image with the combined dimensions
        combined_image = Image.new(
            'RGB', (combined_width, combined_height), "white")

        # Paste the images one below the other
        combined_image.paste(image1, (0, 0))
        combined_image.paste(image2, (0, height1))

        combined_image.save(os.path.join(figures_dir, 'OH2_Figure_S2.png'))

        # Delete the original images
        os.remove(image1path)
        os.remove(image2path)
        
except (FileNotFoundError):
    print(f"Files not found at {image1path} or {image2path}. Please ensure the files exist.")
    pass
