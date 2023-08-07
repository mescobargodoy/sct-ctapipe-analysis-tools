"""
Basic example of reconstruction using ctapipe.
"""
from ctapipe_funcs import open_simtel
from helpers import name_file
from sct_imageprocessor import config_2pass
from ctapipe.io import DataWriter
from ctapipe.image import ImageProcessor
from ctapipe.reco import ShowerProcessor

path = 'simtel_files/gamma_20deg_0deg_run838___cta-prod3-sct_desert-2150m-Paranal-SCT.simtel.gz'
output = name_file(path,'reco_test.h5')

source, calibrator = open_simtel(path)

sct_2pass_config = config_2pass()

image_processor = ImageProcessor(
    subarray=source.subarray, config=sct_2pass_config
)
shower_processor = ShowerProcessor(subarray=source.subarray)

with DataWriter(source, output_path=output, overwrite=True, write_showers=True) as writer:
    for event in source:
        energy = event.simulation.shower.energy
        n_telescopes_r1 = len(event.r1.tel)
        event_id = event.index.event_id
        print(f"Id: {event_id}, E = {energy:1.3f}, Telescopes (R1): {n_telescopes_r1}")

        calibrator(event)
        image_processor(event)
        shower_processor(event)

        stereo = event.dl2.stereo.geometry["HillasReconstructor"]
        if stereo.is_valid:
            print("  Alt: {:.2f}°".format(stereo.alt.deg))
            print("  Az: {:.2f}°".format(stereo.az.deg))
            print("  Hmax: {:.0f}".format(stereo.h_max))
            print("  CoreX: {:.1f}".format(stereo.core_x))
            print("  CoreY: {:.1f}".format(stereo.core_y))
            print("  Multiplicity: {:d}".format(len(stereo.telescopes)))
            writer(event)