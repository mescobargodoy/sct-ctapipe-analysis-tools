"""
Author
@ Miguel Escobar Godoy

Useful ctapipe based functions for opening and manipulating
SCT prod3b simtel.io (.simtel) files.

"""
from ctapipe.io import EventSource
from ctapipe.visualization import CameraDisplay
from ctapipe.calib import CameraCalibrator
from ctapipe.image import dilate
from ctapipe.io import read_table
import matplotlib.pyplot as plt
#import numpy as np
#import sys


def open_simtel(file_name, focal_length='EQUIVALENT', back_seeking=False, calibration=False, data_reducer=False):
    """ 
    Opens .simtel file using ctapipe with optional calibrator class and optional data
    reducing based on 2 pass cleaning using ctapipe tail_cuts and dilation modules.
    Calibrator fills DL0 and DL1 containers. 
    It loads all events in simulation file.

    Parameters

    :param file_name: path to .simtel file
    :type file_name: event.io format but other formats possible

    :param focal_length: Focal lenght. Ideally it should be EFFECTIVE but prod3 can only take EQUIVALENT, defaults to 'EQUIVALENT'
    :type focal_length: str, optional

    :param back_seeking: require the event source to be backwards seekable, defaults to False
    :type back_seeking: bool, optional

    :param calibration: True to apply calibration to events, defaults to False
    :type calibration: bool, optional

    :param data_reducer: optional data volume reducer (not implemented yet), defaults to False
    :type data_reducer: bool, optional

    Returns

    :return: EventSource class from ctapipe
             Camera Calibrator class from ctapipe (optional)
    :rtype: EventSource, CameraCalibrator class

    """

    source = EventSource(
        file_name, focal_length_choice=focal_length, back_seekable=back_seeking)
    if calibration == True:
        subarray = source.subarray
        if data_reducer == True:
            calibrator = CameraCalibrator(
                subarray=subarray, data_volume_reducer="TailCutsDataVolumeReducer")
            return source, calibrator
        else:
            calibrator = CameraCalibrator(subarray=subarray)
            return source, calibrator
    return source


def dilation(geom, mask, N=1):
    """    
    Repeats the process of dilate N times. Default is to do once.
    dilate: Add one row of neighbors to the True values of a pixel mask and return
    the new mask.This can be used to include extra rows of pixels in a mask that was
    pre-computed, e.g. via `tailcuts_clean`.

    Parameters

    :param geom: Camera geometry information
    :type geom: `~ctapipe.instrument.CameraGeometry`

    :param mask: ndarray input mask to be dilated
    :type mask: boolean array

    :param N: number of times to do dilation process, defaults to 1
    :type N: int, optional

    Returns

    array of size number of pixels in camera
    :rtype: boolean array
    """

    assert N > 0, f"number greater than 0 expected, got: {N}"
    i = 0
    while i < N:
        mask = dilate(geom, mask)
        i += 1
    return mask


def disp_images(simtel, hdf5, path_in_table, index_to_disp=(0, 4), two_rows=False):
    """
    Displays 4 or 8 DL1 images with no cleaning from an hdf5 file.

    Parameters
    :param simtel: simulation file to retrieve camera geometry information
    :type simtel: string (simtel.io)

    :param hdf5: path to hdf5 DL1 container file
    :type hdf5: h5

    :param path_in_table: path to table within h5 file ie could write 'dl1/images'
                          will depend on how user defined it
    :type path_in_table: string

    :param im_index_to_disp: indices of images to set with a step of 4
    :type im_index_to_disp: tuple

    :param two_rows: whether two display 4 or 8 images
    :type two_rows: boolean

    """

    a, b = index_to_disp[0], index_to_disp[1]
    assert (abs(
        a-b) == 4), f"Can only display 4 images followed from one another. ie (a-b)%4=0"

    # Indices of images to display
    indices = range(a, b)
    # Retrieving geometry of telescope camera. Assumes prod3b file is being opened.
    source = EventSource(simtel, focal_length_choice='EQUIVALENT',
                         back_seekable=True, max_events=1)
    subarray = source.subarray
    geom = subarray.tel[1].camera.geometry

    # Opening hdf5 file
    dl1_ = read_table(hdf5, path_in_table)
    dl1_im = dl1_['image']

    # Displaying 8 images - 4 per row
    if two_rows == True:
        indices_ = range(a, b+4)
        f, axs = plt.subplots(2, 4, figsize=(12, 5), dpi=150)
        for image, ax in zip(dl1_im[indices_], axs.flatten()):
            square_image = geom.image_to_cartesian_representation(image)
            center = [-square_image.shape[1]/2., square_image.shape[1] /
                      2., -square_image.shape[0]/2., square_image.shape[0]/2.]
            ax.imshow(square_image, extent=center)
        f.tight_layout()

    # Displaying 4 images - single row
    else:
        f, axs = plt.subplots(1, 4, figsize=(12, 5), dpi=150)
        for image, ax in zip(dl1_im[indices], axs.flatten()):
            square_image = geom.image_to_cartesian_representation(image)
            center = [-square_image.shape[1]/2., square_image.shape[1] /
                      2., -square_image.shape[0]/2., square_image.shape[0]/2.]
            ax.imshow(square_image, extent=center)
        f.tight_layout()


def disp_clean(simtel, hdf5, path_in_table, index_to_disp=(0, 4), thresh=0., disp_pixels=False):
    """
    Displays 4 unclean images vs 4 corresponding images with photoelectron cut.
    It will display 4 images based on the index set to display. So if index_to_disp=6
    it will display first images 2-6. If disp_pixels is enabled, it will display
    ctapipe-based plots showing pixels that survived pe cut.

    Parameters
    :param simtel: simulation file to retrieve camera geometry information
    :type simtel: string (simtel.io)

    :param hdf5: path to hdf5 DL1 container file
    :type hdf5: h5

    :param path_in_table: path to table within h5 file ie could write 'dl1/images'
                          will depend on how user defined it
    :type path_in_table: string

    :param im_index_to_disp: indices of images to set with a step of 4
    :param im_index_to_disp: tuple

    :param thresh: pe cut
    :type thresh: float, defaults to 0
    """
    a, b = index_to_disp[0], index_to_disp[1]
    assert (abs(
        a-b) == 4), f"Can only display 4 images followed from one another. ie (a-b)%4=0"

    # Indices of images to display
    indices = range(a, b)
    # Retrieving geometry of telescope camera. Assumes prod3b file is being opened.
    source = EventSource(simtel, focal_length_choice='EQUIVALENT',
                         back_seekable=True, max_events=1)
    subarray = source.subarray
    geom = subarray.tel[1].camera.geometry

    # Need 2 copies of file to not overwrite one another
    dl1_ = read_table(hdf5, path_in_table)
    dl1__ = read_table(hdf5, path_in_table)

    # Assumption here hdf5 files are Dl1 containers only so  there should be an image column
    dl1_im = dl1_['image']

    dl1_thresh = dl1__['image']

    # Applying pe cut
    for image in dl1_thresh:
        mask = image > thresh
        image[~mask] = 0

    # Empty arrays to store images
    no_clean = []
    pe_thresh = []

    # Storing images in arrays
    for im_1 in dl1_im[indices]:
        no_clean.append(im_1)
    for im_2 in dl1_thresh[indices]:
        pe_thresh.append(im_2)
    # Array of images to display
    to_display = no_clean+pe_thresh

    # Plotting starts here
    if disp_pixels == True:
        f, axs = plt.subplots(2, 4, figsize=(12, 5), dpi=150)
        cumulative = 0

    # Displaying ctapipe-based plots
        for image, ax in zip(to_display, axs.flatten()):
            disp = CameraDisplay(geom, image=image, ax=ax, cmap='viridis')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
            if cumulative > 3:
                disp.highlight_pixels(image > thresh, color='green')
                disp.cmap = 'gray'
            cumulative += 1
        f.tight_layout()

    # Displaying more numpy-like 2D histogram plots
    else:
        f, axs = plt.subplots(2, 4, figsize=(12, 5), dpi=150)

        for image, ax in zip(to_display, axs.flatten()):
            square_image = geom.image_to_cartesian_representation(image)
            center = [-square_image.shape[1]/2., square_image.shape[1] /
                      2., -square_image.shape[0]/2., square_image.shape[0]/2.]
            ax.imshow(square_image, extent=center)
        f.tight_layout()
