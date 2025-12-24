# TESS_fullsky_map
combine with [Ethan Kruse's code](https://github.com/ethankruse/tess_fullsky), I provided a code to create fullsky TESS map with scientific units. The map can be further used for extra-galactic background light (EBL) and low surface brightness objects research.

## Image- Fullsky_original

This map illustrates the data prior to star masking, which is dominated by high-intensity signals from stellar sources. The script "fullsky.py" is developed based on Ethan's code.

## Image- Fullsky_star_sm

This map displays the result after masking bright stars(gmag>18). Once the stars are removed, diffuse structures such as the Zodiacal light, M31, and the Magellanic Clouds become clearly visible. This map was generated using "fullsky.py" combined with a custom parallel processing package I developed for efficient masking.

## Image- Fullsky_5deg & residual map

To efficiently remove zody, I applied a high-pass filter (smoothing) to the star-removed map. The algorithm was from my supervisor YK, therefore not showing in the repository. But the residual map is the field that I am using to stack objects.

## Script- Star masking

This script aligns Gaia star catalogs with the TESS HEALPix map and performs the masking process. It is called and executed within "full_sky.py"


## Script- Full_sky
This script generates several .npz files for each TESS sector with HealPix ID, counts(weight), and value of the pixel. Since TESS observing strategy leads to large part of overlapping, so the weighted average is necessary to be done.
I altered Ehthan's code mainly between line 1595-1663, adding star masking and weighting information to later conbine together by weight(not included in this repo).


## Script- Stacking

Using the high-pass filtered map, I can stack whatever target I want, for example
- point source for PSF
- Legacy galaxies for flux calibration
- DESI DR1 group catalog

This code stack the DESI groups to investigate the light profile at different reshift and mass bins.


## ASROC poster

The poster is a [publication](https://indico.phys.nthu.edu.tw/event/143/contributions/615/contribution.pdf) in annual meeting of *The Astronomical Society of the Republic of China(ASROC)*. In the poster I gave stacking results of low luminosity galaxy groups.


