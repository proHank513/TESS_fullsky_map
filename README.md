# TESS_fullsky_map
combine with [Ethan Kruse's code](https://github.com/ethankruse/tess_fullsky), I provided a code to create fullsky TESS map with scientific units. The map can be further used for extra-galactic background light (EBL) and low surface brightness objects research.

## Fullsky_original

This is the map before star masking, showing a high dependency on where stars are. The code "fullsky.py" is developed based on Ethan's code.

## Fullsky_star_sm

This map is after Gaia bright stars (m>18) removal. One can easily observe diffuse structures like zodiacal light, M31, Magellanic Clouds, etc. This is in the same code "fullsky.py" but import parallel processing package written by myself.

## Fullsky_5deg & residual map

To efficiently remove zody, I applied a high-pass filter (smoothing) to the star-removed map. The code is from my supervisor YK, hence not showing here. But the residual map is the map that I am using to stack objects.

## Star masking CODE

This code aligns the Gaia stars with TESS HealPix map then remove them. This code is called and run in "full_sky.py"


## Full_sky CODE
This code generates several npz files for each TESS sector with HealPix ID, counts(weight), and value of the pixel. Since TESS observing strategy gives large part of overlapping, so the weighted average is necessary.
I altered Ehthan's code mainly between line 1595-1663, adding star masking and weighting information to later conbine together by weight(I did not put the code here).


## Stacking CODE

Using the high-pass filtered map, I can stack whatever target I want, for example
- point source for PSF
- Legacy galaxies for flux calibration
- DESI DR1 group catalog

This code stack the DESI groups to see the light profile at different reshift and mass bins.


## ASROC poster

The poster is a [publication](https://indico.phys.nthu.edu.tw/event/143/contributions/615/contribution.pdf) in annual meeting of *The Astronomical Society of the Republic of China(ASROC)*. In the poster I gave stacking results of low luminosity galaxy groups.