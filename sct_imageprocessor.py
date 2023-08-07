from traitlets.config import Config

def config_2pass(image_pixel=4., boundary_pixel=2., min_neighbors=0,keep_isolated_pixels=False):
    """
    Generates ctapipe configuration file for performing 2 pass cleaning and hillas parameterization. 

    :param image_pixel: threshold above which all pixels are retained, defaults to 4
    :type image_pixel: float, optional

    :param boundary_pixel: threshold above which pixels are retained if they have a 
                           neighbor already above the picture_thresh, defaults to 2
    :type boundary_pixel: float, optional
    
    :param min_neighbors: _description_, defaults to 0
    :type min_neighbors: int, optional
    
    :param keep_isolated_pixels: If True, pixels above the picture threshold will be 
                                included always, if not they are only included if a 
                                neighbor is in the picture or boundary, defaults to False.
    :type keep_isolated_pixels: bool, optional
    
    :return: returns configuration file for ctapipe based image processing
    :rtype: ctapipe Config class
    """
    
    image_processor_config = Config(
        {
            "ImageProcessor": {
                "image_cleaner_type": "TailcutsImageCleaner",
                "TailcutsImageCleaner": {
                    "picture_threshold_pe": [
                        ("type", "MST_SCT_SCTCam",image_pixel),
                    ],
                    "boundary_threshold_pe": [
                        ("type", "MST_SCT_SCTCam", boundary_pixel),
                        ],
                    "keep_isolated_pixels": keep_isolated_pixels,
                    "min_picture_neighbors": min_neighbors,
                }
            }
        }
    )

    return image_processor_config

def config_3pass(image_pixel=4., boundary_pixel=2., min_neighbors=1,keep_isolated_pixels=False):
    """
    Generates ctapipe configuration file for performing 3 pass cleaning and hillas parameterization. 

    :param image_pixel: threshold above which all pixels are retained, defaults to 4
    :type image_pixel: float, optional

    :param boundary_pixel: threshold above which pixels are retained if they have a 
                           neighbor already above the picture_thresh, defaults to 2
    :type boundary_pixel: float, optional
    
    :param min_neighbors: _description_, defaults to 0
    :type min_neighbors: int, optional
    
    :param keep_isolated_pixels: If True, pixels above the picture threshold will be 
                                included always, if not they are only included if a 
                                neighbor is in the picture or boundary, defaults to False.
    :type keep_isolated_pixels: bool, optional
    
    :return: returns configuration file for ctapipe based image processing
    :rtype: ctapipe Config class
    """
    
    image_processor_config = Config(
        {
            "ImageProcessor": {
                "image_cleaner_type": "MARSImageCleaner",
                "MARSImageCleaner": {
                    "picture_threshold_pe": [
                        ("type", "MST_SCT_SCTCam",image_pixel),
                    ],
                    "boundary_threshold_pe": [
                        ("type", "MST_SCT_SCTCam", boundary_pixel),
                        ],
                    "keep_isolated_pixels": keep_isolated_pixels,
                    "min_picture_neighbors": min_neighbors,
                }
            }
        }
    )

    return image_processor_config

def config_fact(image_pixel=4., boundary_pixel=2., min_neighbors=1, time_lim=3, keep_isolated_pixels=False):
    """
    Generates ctapipe configuration file for performing 3 pass cleaning and hillas parameterization. 

    :param image_pixel: threshold above which all pixels are retained, defaults to 4 pe
    :type image_pixel: float, optional

    :param boundary_pixel: threshold above which pixels are retained if they have a 
                           neighbor already above the picture_thresh, defaults to 2 pe
    :type boundary_pixel: float, optional
    
    :param min_neighbors: _description_, defaults to 0
    :type min_neighbors: int, optional
    
    :param time_lim: arrival time limit for neighboring pixels
    :type float, defaults to 3 ns 
    
    :param keep_isolated_pixels: If True, pixels above the picture threshold will be 
                                included always, if not they are only included if a 
                                neighbor is in the picture or boundary, defaults to False.
    :type keep_isolated_pixels: bool, optional
    
    :return: returns configuration file for ctapipe based image processing
    :rtype: ctapipe Config class
    """
    
    image_processor_config = Config(
        {
            "ImageProcessor": {
                "image_cleaner_type": "MARSImageCleaner",
                "MARSImageCleaner": {
                    "picture_threshold_pe": [
                        ("type", "MST_SCT_SCTCam",image_pixel),
                    ],
                    "boundary_threshold_pe": [
                        ("type", "MST_SCT_SCTCam", boundary_pixel),
                        ],
                    "keep_isolated_pixels": keep_isolated_pixels,
                    "min_picture_neighbors": min_neighbors,
                    "time_limit_ns": time_lim,
                }
            }
        }
    )

    return image_processor_config