"""
Useful ctapipe based functions for opening and manipulating
SCT prod3b simtel.io (.simtel) files.
"""
from ctapipe.io import EventSource
from ctapipe.visualization import CameraDisplay
from ctapipe.calib import CameraCalibrator
from ctapipe.image import dilate
from ctapipe.io import read_table
import matplotlib.pyplot as plt


def open_simtel(file_name, focal_length='EQUIVALENT', back_seeking=True, calibration=True, data_reducer=False):
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
