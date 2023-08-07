"""
ctapipe based Functions for writing different data level 
containers to hdf5 files.
"""

from ctapipe.io import EventSource
from ctapipe.calib import CameraCalibrator
from ctapipe.io import HDF5TableWriter
import numpy as np
from ctapipe.image import hillas_parameters, tailcuts_clean
from ctapipe_funcs import open_simtel
from helpers import name_file

def write_sim_params(path):
    """
    Opens simtel file and writes event simulation parameters to hdf5.
    It specifically writes truth information of triggered telescopes
    per event and truth shower information. 

    To read table you will need to know the names of tables i.e.
    tel_sim_image_params, tel_sim_core, shower_sim_params

    :param path: path to simtel file
    :type path: string
    """
    # Opening simtel file
    source, __ = open_simtel(path, back_seeking=True, calibration=True)
    
    # Output file name
    filename = name_file(path, 'simulation_container.h5')
    
    # Class for writing hdf5 files
    with HDF5TableWriter(filename=filename,group_name='simulation',overwrite=True) as writer:
    
        for event in source:
            
            # Triggered telescope indices
            triggered_tel = event.trigger.tels_with_trigger
            # Simulation shower params stored here (energy, coords, etc)
            sim_shower = event.simulation.shower
            # Looping over triggered telescopes
            for tel in triggered_tel:
                sim_tel_image_params = event.simulation.tel[tel]
                sim_tel_core = event.simulation.tel[tel].impact
                # writing to hdf5
                writer.write('tel_sim_image_params',sim_tel_image_params)
                writer.write('tel_sim_core',sim_tel_core)
                writer.write('shower_sim_params', sim_shower)
