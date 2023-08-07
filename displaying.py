"""
Useful ctapipe based functions for displaying SCT images. 
"""

from ctapipe.io import EventSource
from ctapipe.visualization import CameraDisplay
from ctapipe.io import read_table
from ctapipe_funcs import open_simtel
import matplotlib.pyplot as plt


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


def disp_event(simtel, lower_energy_bound=0., upper_energy_bound=1000., disp_timing=True):
    """
    Will display two images SCT camera as a function of charge and as a function of
    arrival time optionally. It will also show telescope index as well as true energy 
    of gamma ray and core distance. 

    :param simtel: simulation file to retrieve camera geometry information
    :type simtel: string (simtel.io)

    :param lower_energy_bound: lower energy bound of gamma ray in TeV to look for, defaults to 0.
    :type lower_energy_bound: float, optional

    :param higher_energy_bound: upper energy bound of gamma ray in TeV to look for, defaults to 1000.
    :type higher_energy_bound: float, optional

    :param disp_timing: Whether to display timing information as well, defaults to True
    :type disp_timing: bool, optional

    """
    print("Out of the events (energies) printed out, it will only store the last event in memory. \n"
          "If you wanted to select one specific event, set more stringent energy bounds. \n"
          "These are the events within {} to {} TeV energy range:".format(lower_energy_bound, upper_energy_bound))
    source, calib = open_simtel(simtel)
    subarray = source.subarray

    # Looping over events in simulation file and finding one with specific energy criteria.
    for event_ in source:
        energy = event_.simulation.shower.energy.value

        if lower_energy_bound < energy < upper_energy_bound:
            calib(event_)
            event = event_
            print(event.simulation.shower.energy)

    # Storing triggered telescope indices
    triggered_tels = event.trigger.tels_with_trigger

    # Calling helper functions that will display the images
    for triggered_tel in triggered_tels:
        # Set the telescope and Geometry
        tel = subarray.tel[triggered_tel]
        geometry = tel.camera.geometry
        charge = event.dl1.tel[triggered_tel].image
        time = event.dl1.tel[triggered_tel].peak_time

        # Displaying charge and timing information
        if disp_timing == True:
            f, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=150)
            disp1 = CameraDisplay(geometry, image=charge,
                                  ax=axs[0], cmap='viridis')
            disp1.add_colorbar(label='charge [pe]')
            axs[0].set_title('SCT {}, {:.1f} TeV event, core {:.0f} m'.format(triggered_tel,
                                                                              event.simulation.shower.energy.value, 
                                                                              event.simulation.tel[triggered_tel].impact.distance.value))
            disp2 = CameraDisplay(geometry, image=time,
                                  ax=axs[1], cmap='viridis')
            axs[1].set_ylabel('')
            axs[1].set_title('SCT {}, {:.1f} TeV event, core {:.0f} m'.format(triggered_tel,
                                                                              event.simulation.shower.energy.value, 
                                                                              event.simulation.tel[triggered_tel].impact.distance.value))
            disp2.add_colorbar(label='arrival time [ns]')
            f.tight_layout()

        # Only displaying charge information
        else:
            f, axs = plt.subplots(1, 1, figsize=(5, 4), dpi=120)
            disp1 = CameraDisplay(geometry, image=charge,
                                  ax=axs, cmap='viridis')
            disp1.add_colorbar(label='charge [pe]')
            axs.set_title('SCT {}, {:.1f} TeV event, core {:.0f} m'.format(triggered_tel,
                                                                           event.simulation.shower.energy.value, event.simulation.tel[triggered_tel].impact.distance.value))
            f.tight_layout()
