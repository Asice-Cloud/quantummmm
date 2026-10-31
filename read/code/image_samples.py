"""Read calibrated curve pixels from the original solution.png (not raw data)."""
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import center_of_mass, label

ROOT = Path(__file__).resolve().parents[2]
# Coordinates of the 0/12 horizontal ticks and 0/1 vertical ticks.
CALIBRATIONS = ((108, 964, 556, 108), (1100, 1957, 556, 108),
                (108, 964, 1161, 714), (1100, 1957, 1161, 714))


def image_samples(index):
    arr = np.asarray(Image.open(ROOT/"solution.png").convert("RGB"))
    left, right, zero, one = CALIBRATIONS[index]
    rows, cols = np.indices(arr.shape[:2])
    inside = ((cols > left+4)&(cols < right-4)
              &(rows >= one-2)&(rows <= zero+2))
    # Exclude the first-panel legend rather than interpret text as curve data.
    if index == 0:
        inside &= ~((cols >= 365)&(cols <= 730)&(rows >= 280)&(rows <= 385))
    r, g, b = arr[..., 0].astype(float), arr[..., 1].astype(float), arr[..., 2].astype(float)
    black = inside&(r < 60)&(g < 60)&(b < 60)
    orange = inside&(r > 150)&(g > 45)&(g < 160)&(b < 80)&(r > 1.5*g)
    out = {}
    for color_name, mask in (("black", black), ("orange", orange)):
        samples = []
        for column in range(left+5, right-4):
            ys = np.flatnonzero(mask[:, column])
            # Steep early transients have multiple pixel values at a column.
            # Omit broad strokes from quantitative comparison.
            if len(ys) and ys.max()-ys.min() <= 10:
                samples.append(((column-left)/(right-left)*1200,
                                (zero-float(np.median(ys)))/(zero-one)))
        out[color_name] = np.asarray(samples)
    blue = inside&(b > r+70)&(g > r+45)&(b > 100)&(g < 170)
    components, count = label(blue)
    peaks = []
    for component in range(1, count+1):
        if np.count_nonzero(components == component) < 30:
            continue
        y, x = center_of_mass(blue, components, component)
        peaks.append(((x-left)/(right-left)*1200, (zero-y)/(zero-one)))
    out["peak_markers"] = np.asarray(peaks)
    return out

