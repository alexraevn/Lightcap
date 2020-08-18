# Lightcap
A proto-library of functions for calculating astronomical differential magnitude values and producing light curves.

Usage example for Lx Ser using coordinates from frames taken 07/31/2020:
```
    import lightcap
```
Instantiate a Lightcurve class:
```
    l = lightcap.Lightcurve(path)
```
Set the target (your variable star) and reference object(s) (non-variable stars) coordinates, aperture radius, and name(s) (optional).
Pass the coordinates using tuples (x,y). Use list of tuples for multiple references [(x1, y1), (x2, y2)].
Pass a radius in pixels. Use int type.
Pass a name for the target and reference(s) (optional). Use str type.
```
    l.set_target((1108.81, 1015.97), 8, "Lx Ser")
    l.set_reference([(1341.0, 938.0), (1230.0, 905.0)], 8, ["Reference a", "Reference b"])
```
Run the read_apertures method to calculate the raw luminosity of the target and reference(s).
This will populate the following attributes:

**l.target**: a list of target luminosity values, one value for each frame in path.

**l.reference**: a list of values (for one reference) or a list of lists (one for each reference) of luminosity values, one value per frame per target.
```
    l.read_apertures()
```
Calculate the differential magnitude by running the differential_magnitude method.
By default, method='average'. This will populate the following attributes:

**l.target_magnitude**: a list of magnitude values for target. For multiple references, each target magnitude is equal to:
> (-2.5)*log_10(target_luminosity/average_of_reference_luminosities)

**l.reference_magnitude**:
  - If using a single reference:
    - This will be a list of 0s, since log_10(1) is 0.
  - Else if using multiple references:
    - This will be a list of lists [[],[],...], with each sublist containing reference magnitude values. Each reference magnitude is equal to:
> (-2.5)*log_10(reference_luminosity/average_of_reference_luminosities)
```
    l.differential_magnitude()
```
Contact email: alex.re95@gmail.com
