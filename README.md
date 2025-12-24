# TESS_fullsky_map
combine with [Ethan Kruse's code](https://github.com/ethankruse/tess_fullsky), I provided a code to create fullsky TESS map with scientific units

## Fullsky_original

This is the map before star masking, showing a high dependency on where stars are. The code "fullsky.py" is developed based on Ethan's code.



## Fullsky_star_sm

This map is after Gaia bright stars (m>18) removal. One can easily observe diffuse structures like zodiacal light, M31, Magellanic Clouds, etc. This is in the same code "fullsky.py" but import parallel processing package written by myself.

## Fullsky_5deg

To efficiently remove zody, I applied a high-pass filter (smoothing) to the star-removed map. The code is from my supervisor YK, hence not showing here. But the residual map is the map that I am using to stack objects.