from src.generate_bdf_files import generate_bdf_files
from pathlib import Path

path = Path(__file__).parent
generate_bdf_files(path, control_mass=0.5, treatment_mass=1.27, xmin=-9, xmax=9, dx=0.1, ymin=-8.5, ymax=8, dy=0.1)
