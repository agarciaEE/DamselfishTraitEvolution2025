import cv2
import numpy as np
from skimage import exposure
import os
import argparse

def underwater_image_correction(image, clahe_clip_limit, clahe_grid_size, saturation_factor, denoise_h):
    # Convert to LAB color space
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)

    # Split the LAB image into L, A, and B channels
    l, a, b = cv2.split(lab)

    # Apply CLAHE to L channel
    clahe = cv2.createCLAHE(clipLimit=clahe_clip_limit, tileGridSize=(clahe_grid_size, clahe_grid_size))
    cl = clahe.apply(l)

    # Merge the CLAHE enhanced L-channel back with A and B channels
    limg = cv2.merge((cl,a,b))

    # Convert back to BGR color space
    enhanced = cv2.cvtColor(limg, cv2.COLOR_LAB2BGR)

    # White balancing
    grayworld = exposure.equalize_adapthist(enhanced)

    # Convert grayworld to 8-bit unsigned integer format
    grayworld_8bit = (grayworld * 255).astype(np.uint8)

    # Increase saturation
    hsv = cv2.cvtColor(grayworld_8bit, cv2.COLOR_BGR2HSV)
    hsv[:,:,1] = np.clip(hsv[:,:,1] * saturation_factor, 0, 255).astype(np.uint8)
    saturated = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)

    # Denoise
    denoised = cv2.fastNlMeansDenoisingColored(saturated, None, denoise_h, denoise_h, 7, 21)

    return denoised

def process_folder(input_folder, output_folder, clahe_clip_limit, clahe_grid_size, saturation_factor, denoise_h):
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Process all PNG images in the input folder
    for filename in os.listdir(input_folder):
        if filename.lower().endswith('.png'):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)

            # Read the image
            img = cv2.imread(input_path)

            # Apply correction
            corrected_img = underwater_image_correction(img, clahe_clip_limit, clahe_grid_size, saturation_factor, denoise_h)

            # Save the corrected image
            cv2.imwrite(output_path, corrected_img)
            print(f"Processed: {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct underwater images in a folder")
    parser.add_argument("input_folder", help="Path to the input folder containing PNG images")
    parser.add_argument("output_folder", help="Path to the output folder for corrected images")
    parser.add_argument("--clahe_clip_limit", type=float, default=2.0, help="CLAHE clip limit (default: 2.0)")
    parser.add_argument("--clahe_grid_size", type=int, default=8, help="CLAHE grid size (default: 8)")
    parser.add_argument("--saturation_factor", type=float, default=1.0, help="Saturation increase factor (default: 1.0)")
    parser.add_argument("--denoise_h", type=float, default=10, help="Denoising strength (default: 10)")
    args = parser.parse_args()

    process_folder(args.input_folder, args.output_folder, args.clahe_clip_limit, args.clahe_grid_size, args.saturation_factor, args.denoise_h)
    print("All images processed successfully.")
