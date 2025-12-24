import os
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy
import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from multiprocessing import Pool

pixel_1d_num = 77


def create_equal_area_grid(center_l, center_b, grid_size=pixel_1d_num, cell_size=17.2):
    
    # cell_size is the size for SINGLE pixel, 1/3 of 51.5 arcsec

    # create xy grid
    x_array = np.linspace(0 - grid_size/2, 0 + grid_size/2, grid_size)
    y_array = np.linspace(0 - grid_size/2, 0 + grid_size/2, grid_size)

    # Create a mesh grid from 2 list to a 2d array
    x_grid, y_grid = np.meshgrid(x_array, y_array)
    
    
    # cell_size in arcsec, independent to map
    cell_size_deg = cell_size / 3600 #change to degree
    
    # Create WCS with ZEA projection
    w = WCS(naxis=2)
    w.wcs.ctype = ['GLON-ZEA', 'GLAT-ZEA']  # ZEA = zenithal equal-area projection
    w.wcs.crpix = [0., 0.]  # Reference pixel
    w.wcs.crval = [center_l, center_b]  # Center of your region in degrees
    w.wcs.cdelt = [-cell_size_deg, cell_size_deg]  # Pixel scale in degrees, The negative value indicates that pixel coordinates increase in the opposite direction of increasing l

    # Transform coordinates
    lb_gird = w.all_pix2world(x_grid, y_grid,0)

    # Separate l and b components
    l_grid = lb_gird[0]  # All rows, all columns, first element l
    b_grid = lb_gird[1]  # All rows, all columns, first element b


    return l_grid, b_grid 


def sample_map(map_data, l_grid, b_grid):
    # change the array into np_id
    sampled_grid = np.zeros_like(l_grid)
    for i in range(l_grid.shape[0]):
        # len = grid_size
        for j in range(l_grid.shape[1]):
            pixel_num = hp.ang2pix(nside = 4096, theta = l_grid[i, j], phi = b_grid[i, j], lonlat = True, nest = False)
            sampled_grid[i, j] = map_data[pixel_num] 

    return sampled_grid


def stack_galaxies_3d_array(galaxy_coords, map_data):
    # unzip the coords
    galaxy_coords_array = np.array(galaxy_coords)

    galaxy_l = galaxy_coords_array[:, 0]
    galaxy_b = galaxy_coords_array[:, 1]
    
    stacked_grid = []
    
    for galaxy_num in range(len(galaxy_coords)):
        print('now processing ', galaxy_num, 'th galaxy')
        l_grid, b_grid = create_equal_area_grid(galaxy_l[galaxy_num], galaxy_b[galaxy_num])
        sampled_grid = sample_map(map_data, l_grid, b_grid)
        # Append the new data to the stacked grid
        stacked_grid.append(sampled_grid)
    
    # Convert the list of arrays into a 3D NumPy array
    stacked_grid_array = np.array(stacked_grid)


    return  stacked_grid_array

#=============================================================================================


from astropy.io import fits
# Paths
fits_path = "/Users/po-hanchen/Desktop/ASIAA/Fullsky_Map/IHL/desi_glx_gp/desi_gp_z05_gcoord.fits"
npz_path = '/Users/po-hanchen/Desktop/ASIAA/Fullsky_Map/IHL/fullsky_masking_0630_MJr_sr-1.npz'

# Load FITS data
hdul = fits.open(fits_path)
data = hdul[1].data

# Extract columns
glon_all = data['GLON']
glat_all = data['GLAT']
redshift = data['GRP_Z']
mass = data['GRP_LOGM']
rich = data['RICH']

r_mask = (rich > 1) & (rich < 11)

r_select_glon = glon_all[r_mask]
r_select_glat = glat_all[r_mask]
r_select_redshift = redshift[r_mask]
r_select_mass = mass[r_mask]



hdul.close()

# Load the TESS map
loaded_array = np.load(npz_path)
tess_hp4096 = loaded_array['pix_value']

#=============================================================================================

z_bins = np.arange(0.00, 0.50, 0.05)  # e.g., 0.00, 0.05, ...
m_bins = np.arange(11.0, 14.0, 0.5)  # e.g., 11.00, 11.25, ...
low_limit = 1000  # skip bins with fewer galaxies than this
n_cpus = 4

output_dir = '/Users/po-hanchen/Desktop/ASIAA/Fullsky_Map/IHL/gp_stack/stack_gp_0813'
os.makedirs(output_dir, exist_ok=True)

def process_bin(args):
    z_min, z_max, m_min, m_max = args
    bin_label = f"z_{z_min:.2f}_{z_max:.2f}__m_{m_min:.2f}_{m_max:.2f}"

    mask = (
        (r_select_redshift >= z_min) & (r_select_redshift < z_max) &
        (r_select_mass >= m_min) & (r_select_mass < m_max)
    )
    glon_bin = r_select_glon[mask]
    glat_bin = r_select_glat[mask]

    if len(glon_bin) < low_limit:
        print(f"Skipped {bin_label} — {len(glon_bin)} galaxies.")
        return None

    galaxy_coords = np.array(list(zip(glon_bin, glat_bin)))
    print(f"Stacking {len(galaxy_coords)} galaxies in bin {bin_label}...")

    sample_3d_array = stack_galaxies_3d_array(galaxy_coords, tess_hp4096)

    output_file = os.path.join(output_dir, f'stacked_gp_{bin_label}_mask.npz')
    np.savez_compressed(output_file, sample_3d_array=sample_3d_array)
    print(f"Saved {bin_label}")
    return bin_label

# Step 1: generate task list
tasks = []
for z_min in z_bins:
    z_max = z_min + 0.05
    for m_min in m_bins:
        m_max = m_min + 0.5
        tasks.append((z_min, z_max, m_min, m_max))

# Step 2: run in parallel
if __name__ == '__main__':
    with Pool(processes=n_cpus) as pool:
        pool.map(process_bin, tasks)

print("All bins processed!")