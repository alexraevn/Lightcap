
import os, glob
from photutils import aperture_photometry, CircularAperture
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

class Lightcurve:

    def __init__(self, path):
        """
        Lightcurve class containing methods for differential magnitude computation.

        Attributes:
        self.target & self.reference: list of aperture photometry luminosity values
        self.target_magnitude: list of differential magnitude values for target.
        self.reference_magnitude: list (for a single reference) or list of lists (one list for each reference) of differential magnitude values.
        self.target_pos & self.reference_pos: indicate x and y coordinates of given target and reference(s) objects.

        Arguments:
        path: str, path to a directory which contains all the fits files to be read.

        Example:
        l = lightcap.Lightcurve("/path/to/all/files/to/be/used") # Will read all files in this directory
        l = lightcap.Lightcurve("/path/to/many/files/a-wcs-reduced-*.fit") # Will only read the specified files
        """

        if "*" in path:
            fits_path = path
        else:
            fits_path = os.path.join(self.path, "*.fit")

        os.chdir(os.path.dirname(fits_path))
        self.target = [] #target star aperture counts
        self.target_pos = None
        self.reference = [] #reference and checks aperture counts
        self.reference_pos = None
        self.fits_list = []

        for i in glob.glob(fits_path):
            self.fits_list.append(i)
        self.fits_list.sort()

        self.hdu_list = []
        self.jd_axis = []
        for i in self.fits_list:
            hdu = fits.open(i)[0]
            self.hdu_list.append(hdu)
            self.jd_axis.append(hdu.header["JD"])


    def set_target(self, position, radius, name=None):
        """
        Configure the target object\'s position, aperture radius, and name.

        Arguments:
        position: tuple of len 2, x and y coordinates of target.
        radius: int, aperture radius in pixels.
        name: optional, str, name of the target object

        Example:
        set_target((1108.81, 1015.97), 8, "Lx Ser")
        """

        if type(position) == tuple and len(position) == 2:
            self.target_pos = position
            self.target_aperture = CircularAperture(position, r=radius)
            if type(name) == str:
                self.target_name = name
            self.target = [] #running this function resets the aperture photometry values of target.
        else:
            print("Position must be provided as a tuple of len 2.\nExample: (1108.81, 1015.97)")

        return

    def set_reference(self, position, radius, name=None):
        """
        Configure reference or list of references to use in differential magnitude calculation.

        Arguments:
        position: tuple or list of tuples of len 2, x and y coordinates of references.
        radius: int, aperture radius in pixels.
        name: optional, str or list of str, name(s) of reference(s).

        Example:
        set_reference([(1341.0, 938.0), (1230.0, 905.0)], 8, ["a", "b"])
        """

        if type(position) == tuple and len(position) == 2:
            self.number_of_references = 1
            self.reference_pos = position
            self.reference_aperture = CircularAperture(position, r=radius)
            if type(name) == str or type(name) == list:
                self.reference_name = name
            self.reference = []
        elif type(position) == list and len(position[0]) == 2:
            self.number_of_references = len(position)
            self.reference_pos = position
            self.reference_aperture = CircularAperture(position, r=radius)
            self.reference = [[] for x in range(len(position))]
            if type(name) == str or type(name) == list:
                self.reference_name = name
        else:
            print("Position(s) must be provided as a tuple or list of tuples of len 2.\nExample: (1108.81, 1015.97)")

        return

    def read_apertures(self):
        """
        Read the luminosity of the target and reference objects using the aperture_photometry method of the photutils package.
        List of luminosity values stored in self.target and self.reference attributes.

        Arguments:
        None
        """

        if self.target_pos == None or self.reference_pos == None:

            print("You must set a target and reference(s) before performing aperture photometry.\nUse set_target and set_reference methods.")

        elif self.target != []:

            print("Aperture photometry has been performed already, and target and/or reference(s) were updated.\nPlease update both target and reference before computing new photometry values.")

        else:

            for i in self.hdu_list:
                target_phot_table = aperture_photometry(i, self.target_aperture) #reading target luminosity
                self.target.append(target_phot_table[0][3])

                phot_table = aperture_photometry(i, self.reference_aperture)
                if self.number_of_references == 1: #reading only one reference object
                    self.reference.append(phot_table[0][3])
                elif self.number_of_references > 1: #reading multiple references
                    for r in range(len(self.reference)):
                        self.reference[r].append(phot_table[r][3])

        return

    def differential_magnitude(self, method='average'):
        """
        Compute the differential magnitude using target and reference(s).

        Arguments:
        method='average', (default) average of the luminosity of all references in each frame is subtracted from the target.
        """

        if method == 'average':

            target_diff_mag = self.target.copy()
            reference_diff_mag = self.reference.copy()

            if self.number_of_references == 1:
                print("Average is not computed. Using one reference.")

                for i in range(len(self.target)):
                    target_diff_mag[i] = (-2.5)*np.log10(self.target[i]/self.reference[i])
                    reference_diff_mag[i] = (-2.5)*np.log10(self.reference[i]/self.reference[i]) #0 if using single reference

            if self.number_of_references > 1: #if more than 1 reference is used
                for i in range(len(self.target)): #iterate through each frame
                    values = [] #list of values in a frame (i), one per reference obj

                    for r in range(len(self.reference)): #iterate through each reference object
                        values.append(self.reference[r][i]) #append a value for each reference object

                    avg_value = sum(values)/len(values) #compute the average counts per frame
                    target_diff_mag[i] = (-2.5)*np.log10(self.target[i]/avg_value) #populate the list of target diff magnitudes

                    for r in range(len(self.reference)):
                        reference_diff_mag[r][i] = (-2.5)*np.log10(self.reference[r][i]/avg_value) #populate the list of reference diff magnitudes


        self.target_magnitude = target_diff_mag
        self.reference_magnitude = reference_diff_mag

        return

#--- config path of lights
# path = ("/home/raevn/Documents/CTMO/data/2020-07-31/lights")
# awcsreduced = path+"/a-wcs-reduced-*.fit"
#---

#--- configure aperture
# pos = [(1108.81, 1015.97), (1341.0, 938.0), (1230.0, 905.0), (1235, 1239), (1350, 1212)]
# aperture = CircularAperture(pos, r=8)
# [(LxSer), ...]

#--- read fits
# object_list = []
# os.chdir(path)
#
# for i in glob.glob(awcsreduced):
#     object_list.append(i)
#
# print(len(object_list))
#
# hdu_list = []
# for i in object_list:
#     hdu_list.append(fits.open(i)[0])
# #---
#
# #--- perform aperture aperture photometry
# master = [[],[],[],[],[]]
#
# for i in hdu_list:
#     phot_table = aperture_photometry(i, aperture)
#     master[0].append(phot_table[0][3])
#     master[1].append(phot_table[1][3])
#     master[2].append(phot_table[2][3])
#     master[3].append(phot_table[3][3])
#     master[4].append(phot_table[4][3])
# #---
#
# #--- subtract magnitudes
# lx = master[0]
# a = master[1]
# b = master[2]
# c = master[3]
# d = master[4]
#
# #configure
# subtract = None
# logsubtract = a
# magnitude = None
#
# if subtract == a:
#     for i in range(len(lx)):
#         lx[i] = lx[i]-a[i]
#         b[i] = b[i]-a[i]
#         a[i] = a[i]-a[i]
#
# if subtract == b:
#     for i in range(len(lx)):
#         lx[i] = lx[i]-b[i]
#         a[i] = a[i]-b[i]
#         b[i] = b[i]-b[i]
#
#
# if magnitude == a:
#     for i in range(len(lx)):
#
#         lx[i] = (-2.5)*np.log10(lx[i])
#         b[i] = (-2.5)*np.log10(b[i])
#         c[i] = (-2.5)*np.log10(c[i])
#         a[i] = (-2.5)*np.log10(a[i])
#
# if logsubtract == a:
#     for i in range(len(lx)):
#         avg = (a[i]+b[i]+d[i])/3.0
#         lx[i] = (-2.5)*np.log10(lx[i]/avg)
#         a[i] = (-2.5)*np.log10(a[i]/avg)
#         #c[i] = (-2.5)*np.log10(c[i]/(0.5*(a[i]+b[i])))
#         b[i] = (-2.5)*np.log10(b[i]/avg)
#         d[i] = (-2.5)*np.log10(d[i]/avg)
#
# if logsubtract == c:
#     for i in range(len(lx)):
#         lx[i] = (-2.5)*np.log10(lx[i]/c[i])
#         a[i] = (-2.5)*np.log10(a[i]/c[i])
#         b[i] = (-2.5)*np.log10(b[i]/c[i])
#         c[i] = (-2.5)*np.log10(c[i]/c[i])
#
# #--- plot
# x = range(len(object_list))
# print(len(object_list))
#offset, slope = np.polynomial.polynomial.polyfit(x, lx, 1)

# if __name__ == "__main__":
#
#     plt.plot(x, lx, 'bo', label="Lx Ser")
#     plt.plot(x, a, 'ro', label="Reference") #15:38:16.47 +18:48:32:97
#     plt.plot(x, b, 'mo', label="Check a")
#     plt.plot(x, d, 'o', label="Check d")
#     #plt.plot(x, offset + slope*x, 'c-', label="Linear Fit")
#     #plt.plot(x, c, 'go', label="Check b")
#     plt.legend()
#     plt.xlabel('Image Number')
#     plt.ylabel('Differential Magnitude')
#     plt.title('Lx Ser 2020-08-03 ~03:05-5:00 UTC')
#     plt.show()
