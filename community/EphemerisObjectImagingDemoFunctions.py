#!/usr/bin/env python
# coding: utf-8

# # Utility functions for Ephemeris Object Imaging
# This nootebook defines various utility functions used in the demo notebooks for the ephemeris object imaging.

# In[1]:


# extra packages needed
# import os
# os.system('pip install astropy')
# for plotting images
# os.system('pop install aplpy')


# In[11]:


from casatasks import getephemtable, tclean, exportfits, listobs, clearstat, imstat
from casatools import measures, quanta, table, image, msmetadata
import os
import aplpy
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot
import numpy as np
import pylab as pl
import shutil

me = measures()
qa = quanta()
tb = table()
ia = image()
msmd = msmetadata()


# ## Functions


def delete_tcleanimages(imagename):
    """clean up images from a previous run"""
    imtypes = ["residual", "image", "psf", "pb", "model", "sumwt", "mask", "weight"]
    for type in imtypes:
        if os.path.exists(f"{imagename}.{type}"):
            os.system(f"rm -rf {imagename}.{type}")


# def get_attachedEphemtablepath(msfile, fieldid):
#    ''' returns full path of attached eph table for a specified field ID'''
#    import glob
#    from casatools import table
#    _tb = table()
#    _tb.open(msfile+'/FIELD')
#    ephemid = _tb.getcell('EPHEMERIS_ID',fieldid)
#    _tb.close()
#    tabname = glob.glob(f'{msfile}/FIELD/EPHEM{ephemid}*/')
#    return os.path.abspath(tabname[0])


def get_attachedEphemtablepath(msfile, fieldname):
    """returns full path of attached eph table for a specified field ID"""
    from casatools import table, msmetadata
    import glob

    _tb = table()
    _msmd = msmetadata()
    _msmd.open(msfile)
    fieldid = _msmd.fieldsforname(fieldname)[0]
    _msmd.close()
    _tb.open(msfile + "/FIELD")
    ephemid = _tb.getcell("EPHEMERIS_ID", fieldid)
    _tb.close()
    tabname = glob.glob(f"{msfile}/FIELD/EPHEM{ephemid}*/")
    return tabname[0]


# Convert a CASA image to FITS image for displaying the image via aplpy


# Define functions to use in displaying the image
def displayImage(imgname, markers={}, title=""):
    """convert casa image to fits and display the fits image"""
    # export to a FITS image (overwrite if the output exists)
    from casatools import quanta

    _qa = quanta()
    fitsimagename = imgname + "fits"
    exportfits(imagename=imgname, fitsimage=fitsimagename, overwrite=True)
    fits.setval(fitsimagename, "TIMESYS", value="utc")
    fig = pl.figure(figsize=(15, 15))
    img = aplpy.FITSFigure(fitsimagename, subplot=[0.1, 0.1, 0.5, 0.5])
    img.show_colorscale()
    if markers != {}:
        for ky, pos in markers.items():
            if "ra" in pos and "dec" in pos:
                print(f"ky={ky}")
                if "ext" in ky:
                    color = "red"
                else:
                    color = "yellow"
                img.show_markers(
                    pos["ra"], pos["dec"], edgecolor=color, marker="o", s=10, alpha=0.5
                )
                marker_ra = _qa.time(_qa.quantity(pos["ra"], "deg"), prec=9)[0]
                marker_dec = _qa.angle(_qa.quantity(pos["dec"], "deg"), prec=9)[0]
                print(f"{ky} at: {marker_ra}, {marker_dec}")
    pl.title(title)
    fig.canvas.draw()
    print("image center:", printImageCenter(imgname))


# def displayImage(imgname, markers={}, title=''):
#    ''' convert casa image to fits and display the fits image'''
#    # export to a FITS image (overwrite if the output exists)
#    fitsimagename = imgname+'fits'
#    exportfits(imagename=imgname, fitsimage=fitsimagename, overwrite=True)
#    fits.setval(fitsimagename, 'TIMESYS', value='utc')
#    fig = pl.figure(figsize=(15,15))
#    img = aplpy.FITSFigure(fitsimagename, subplot=[0.1,0.1,0.5,0.5])
#    img.show_colorscale()
#    if markers!={}:
#       img.show_markers(markers['ra'], markers['dec'], edgecolor='red', marker='o',s=10, alpha=0.5)
#       marker_ra = qa.time(qa.quantity(markers['ra'],'deg'),prec=9)[0]
#       marker_dec = qa.angle(qa.quantity(markers['dec'],'deg'),prec=9)[0]
#       print(f'marker at: {marker_ra}, {marker_dec}')
#    pl.title(title)
#    fig.canvas.draw()
#    print('image center:',printImageCenter(imgname))


def printImageCenter(imgname):
    from casatools import image, quanta

    _ia = image()
    _qa = quanta()
    _ia.open(imgname)
    shape = _ia.shape()
    csys = _ia.coordsys()
    center = csys.toworld([shape[0] / 2.0, shape[1] / 2.0])
    _ia.done()
    return (
        _qa.time(_qa.quantity(center["numeric"][0], "rad"), prec=9)[0],
        _qa.angle(_qa.quantity(center["numeric"][1], "rad"), prec=9)[0],
    )


def ephem_dir(ephemtab, refep, observatory):
    from casatools import measures

    _me = measures()
    _me.framecomet(ephemtab)
    _me.doframe(_me.observatory(observatory))
    _me.doframe(_me.epoch("utc", refep))
    return _me.measure(_me.measure(_me.direction("COMET"), "AZELGEO"), "ICRS")


def dispAstropy(imname="sim_onepoint_true.im", chanslice=0):
    exportfits(imagename=imname, fitsimage=imname + ".fits", overwrite=True)
    hdu = fits.open(imname + ".fits")[0]
    # print(hdu.header)
    wcs = WCS(hdu.header, naxis=2)
    fig = pl.figure()
    fig.add_subplot(121, projection=wcs)
    pl.imshow(hdu.data[0, chanslice, :, :], origin="lower", cmap=pl.cm.viridis)
    pl.xlabel("RA")
    pl.ylabel("Dec")


def dispImage(imname="sim_onepoint_true.im", chanslice=0, useAstropy=False):
    from casatools import image

    _ia = image()
    _ia.open(imname)
    pix = _ia.getchunk()
    shp = _ia.shape()
    _ia.close()
    pl.figure(figsize=(20, 10))
    pl.clf()
    # if shp[3]>1:
    if useAstropy == False:
        pl.subplot(121)
        pl.imshow(pix[chanslice, 0, :, :])
        pl.title(f"Image from channel {chanslice}")
    else:
        dispAstropy(imname, chanslice)
        # displayImage(imname)
    if shp[3] > 1:
        pl.subplot(122)
        ploc = np.where(pix == pix.max())
        pl.plot(pix[ploc[0][0], ploc[1][0], 0, :])
        pl.title("Spectrum at source peak")
        pl.xlabel("Channel")


def casadatetoastropydate(casadatestring):
    datelist = casadatestring.split("/")
    retstr = ""
    if len(datelist) == 4:
        retstr = datelist[0] + "-" + datelist[1] + "-" + datelist[2] + "T" + datelist[3]
    return retstr
