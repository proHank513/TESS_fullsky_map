import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from astropy.table import Table
import matplotlib.pyplot as plt
import os


def load_fits_data(fits_file):
    with fits.open(fits_file) as hdul:
        data = hdul[1].data
        header = hdul[1].header
    return data, WCS(header)

def load_star_catalog(catalog_file):
    star_catalog = Table.read(catalog_file)
    return star_catalog

def healpix_filter_stars(data, wcs, star_catalog, nside=8):
    
    # Step 1: Find the center of the TESS image and its HEALPix pixel ID
    center_x, center_y = data.shape[1] // 2, data.shape[0] // 2
    center_world = wcs.pixel_to_world(center_x, center_y)
    center_galactic = center_world.galactic
    center_l, center_b = center_galactic.l.deg, center_galactic.b.deg
    
    center_pixel = hp.ang2pix(nside, center_l, center_b, lonlat=True)
    
    star_pixels = star_catalog['hpid']
    
    # Step 3: Find the nearest 9 pixels
    neighbors = hp.get_all_neighbours(nside, center_pixel)
    nearby_pixels = np.concatenate(([center_pixel], neighbors))
    
    # Filter stars
    mask = np.isin(star_pixels, nearby_pixels)
    filtered_catalog = star_catalog[mask]
    
    print(f"Number of stars after HEALPix filtering: {len(filtered_catalog)}")
    print(f"Center RA, Dec: {center_l}, {center_b}")
    print(f"Center HEALPix pixel: {center_pixel}")
    print(f"Number of stars in catalog: {len(star_catalog)}")
    print(f"Number of stars in nearby pixels: {np.sum(mask)}")

    return filtered_catalog

def galactic_to_cartesian(l, b):
    """Convert galactic coordinates to 3D Cartesian coordinates on a unit sphere."""
    l_rad = np.radians(l)
    b_rad = np.radians(b)
    x = np.cos(b_rad) * np.cos(l_rad)
    y = np.cos(b_rad) * np.sin(l_rad)
    z = np.sin(b_rad)
    return np.column_stack((x, y, z))

def create_pixel_kdtree(data_shape, wcs):
    # use galactic sphere in 3D Cartesian coord
    y, x = np.indices(data_shape)
    pixel_coords = np.column_stack((x.ravel(), y.ravel()))
    world_coords = wcs.pixel_to_world(pixel_coords[:, 0], pixel_coords[:, 1])
    galactic_coords = world_coords.galactic
    cartesian_coords = galactic_to_cartesian(galactic_coords.l.deg, galactic_coords.b.deg)
    return cKDTree(cartesian_coords)

def calculate_radius(magnitude):
    # Convert degrees to radians for spherical calculations
    return np.where(magnitude < 10,
                    np.radians((-16.212 * (magnitude-1) + 231) / 3600),
                    np.radians((-8.757 * (magnitude-1) + 154.14) / 3600))

def mask_stars_chunk(args):
    star_chunk, pixel_tree, image_shape = args

    masking_radius = calculate_radius(star_chunk['g_dered'])

    star_coords = SkyCoord(l=star_chunk['l'], b=star_chunk['b'], unit='deg', frame='galactic')
    star_cartesian = galactic_to_cartesian(star_coords.l.deg, star_coords.b.deg)
    
    indices = pixel_tree.query_ball_point(star_cartesian, r=masking_radius)
    
    mask = np.zeros(image_shape, dtype=bool)
    mask_flat = mask.ravel()

    for idx_list in indices:
        if idx_list:
            mask_flat[idx_list] = True

    return mask

def parallel_mask_stars(data, wcs, star_catalog, num_processes=None):
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()

    filtered_catalog = healpix_filter_stars(data, wcs, star_catalog)

    pixel_tree = create_pixel_kdtree(data.shape, wcs)
    
    chunk_size = len(filtered_catalog) // num_processes
    star_chunks = [filtered_catalog[i:i+chunk_size] for i in range(0, len(filtered_catalog), chunk_size)]
    
    args = [(chunk, pixel_tree, data.shape) for chunk in star_chunks]
    
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        chunk_masks = list(executor.map(mask_stars_chunk, args))
    
    final_mask = np.any(chunk_masks, axis=0)
    
    masked_data = data.copy()
    masked_data[final_mask] = np.nan
    
    return masked_data

def plot_gaia_on_tess(data, wcs, star_catalog, isec, icam, iccd, savepath = None):
    filtered_catalog = healpix_filter_stars(data, wcs, star_catalog)
    filtered_catalog = filtered_catalog[(filtered_catalog['g_dered'] > 5) & (filtered_catalog['g_dered'] < 7)]

    # Convert the filtered star catalog coordinates (l, b) to pixel (x, y) in the TESS image
    star_coords = SkyCoord(l=filtered_catalog['l'], b=filtered_catalog['b'], unit='deg', frame='galactic')
    x, y = wcs.world_to_pixel(star_coords)

    # Filter out stars that are outside the TESS image (x, y) range
    in_image_mask = (x >= 0) & (x < data.shape[1]) & (y >= 0) & (y < data.shape[0])
    x_in_image = x[in_image_mask]
    y_in_image = y[in_image_mask]

    # Plot and save as PNG
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    
    # Plot the TESS image
    ax.imshow(data, origin='lower', cmap='gray', vmin=np.nanpercentile(data, 5), vmax=np.nanpercentile(data, 95))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    # Plot Gaia stars as red circles
    ax.scatter(x_in_image-43, y_in_image, s=10, edgecolor='red', facecolor='none', label='Gaia Stars')
    
    # get the center g_coord of tess image
    center_x, center_y = data.shape[1] // 2, data.shape[0] // 2
    center_world = wcs.pixel_to_world(center_x, center_y)
    center_galactic = center_world.galactic
    center_l, center_b = center_galactic.l.deg, center_galactic.b.deg  # Galactic coordinates (l, b)

    plt.legend()
    plt.title(f'TESS Sector{isec}-{icam}-{iccd} with Gaia Stars. (l,b) = ({int(center_l)},{int(center_b)})', fontsize = 20)

    
    if savepath == None:
        # Get the current working directory
        current_dir = os.getcwd()

        # Combine it with the file name
        savepath = f'{current_dir}/{isec}-{icam}-{iccd}_(l,b)=({int(center_l)},{int(center_b)})_2.png'

    # Save the figure as PNG
    plt.savefig(savepath, format='png')
    #plt.show()
    plt.close()


def main(fits_file, catalog_file, output_file):
    # Load FITS data
    data, wcs = load_fits_data(fits_file)
    
    # Load star catalog
    star_catalog = load_star_catalog(catalog_file)
    
    # Mask stars
    masked_data = parallel_mask_stars(data, wcs, star_catalog, num_processes = 2)
    
    # Save masked data
    with fits.open(fits_file) as hdul:
        hdul[1].data = masked_data
        hdul.writeto(output_file, overwrite=True)

if __name__ == "__main__":
    fits_file = "/data/star_masking/tess49-1-1.fits"
    catalog_file = "/data/star_masking/Gaia_subset_dered_l_b_hpid.fits"
    output_file = "/data/star_masking/masking_healpy/49-1-1_mask_dyn_hp_0815.fits"
    
    main(fits_file, catalog_file, output_file)
