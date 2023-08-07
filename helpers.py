"""
Will contain functions to parse filenames, management of files
and other related functionalities.
"""

def find_characters(word, character):
    """
    Finds the indices of a character in a string.
    https://stackoverflow.com/questions/21199943/index-of-second-repeated-character-in-a-string
    author: mhlester

    :param word: filename to parse
    :type word: string
    
    :param character: _description_
    :type character: string
    
    :return: indices where character shows up
    :rtype: integer or array of integers
    """
    found = []
    last_index = -1
    while True:
        try:
            last_index = word.index(character, last_index+1)
        except ValueError:
            break
        else:
            found.append(last_index)
    return found

def name_file(filepath, filename_ending):
    """
    Returns a newly named file while removing any path dependencies.
    F
    or example:
    filepath = 'simtel_files/gamma_20deg_0deg_run10.simtel.gz'
    filename_ending = 'sim_container.h5'
    name_file(filepath,filename_ending) will return
    >> gamma_20deg_0deg_run10_sim_container.h5

    :param file_path: path to simtel file
    :type file_path: string
    
    :param filename_ending: ending to be specified for file
    :type filename_ending: string
    
    :return: new file name with removal of '/' and '.simtel.gz'
    :rtype: string
    """
    # Finds the index of every instance of '/'
    ind1 = find_characters(filepath, '/')
    # Finds the starting index of '.simtel.gz'
    ind2 = find_characters(filepath, '.simtel.gz')
    # removes all characters up to and including the last '/' as well as .simtel.gz
    new_filename = '{}_{}'.format(filepath[ind1[-1]+1:ind2[0]],filename_ending)
    
    return new_filename
