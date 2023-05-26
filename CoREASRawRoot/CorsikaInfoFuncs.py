import sys
from re import search
import io
import numpy as np
sys.path.append("/home/grand/") # !TODO: make sure to update this path once RawROOT is included in grandlib!!
from grand.io.root_trees import * # this is home/grand/grand (at least in docker) or ../../grand


# read values from SIM.reas or RUN.inp
def find_input_vals(line):
  return search(r'[-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?', line)



def read_params(input_file, param):
  # works for both SIM.reas and RUN.inp, as long as you are looking for numbers
  val = "1111"
  with open(input_file, "r") as datafile:
    for line in datafile:
      if param in line:
        line = line.lstrip()
        if find_input_vals(line):
          val = find_input_vals(line).group()
          print(param, "=", val)
  return float(val)



def read_atmos(input_file):
    # RUN.inp only
    with open(input_file, mode="r") as datafile:
        for line in datafile:
            if "ATMFILE" in line:
                atmos = line.split("/")[-1]
    return atmos



def read_site(input_file):
    # RUN.inp only
    atmos = read_atmos(input_file)
    if "Dunhuang" in atmos:
        site = "Dunhuang"
    elif "Lenghu" in atmos:
        site = "Lenghu"
    else:
        site = atmos
    return site



def read_first_interaction(log_file):
    """
    For now, get this from the log.
    (Yes, this is ugly and unreliable, but it's good enough for now)
    I want to change the Coreas output so this is included in the reas file instead.
    """
    with open(log_file, mode="r") as datafile:
        for line in datafile:
            if "height of first interaction" in line:
                val = line.split("interaction")[-1]
                first_interaction = find_input_vals(val)
    return first_interaction



def antenna_positions_dict(pathAntennaList):
    """
    get antenna positions from SIM??????.list and store in a dictionary
    .list files are structured like "AntennaPosition = x y z name"

    """
    antennaInfo = {} # store info in a dict

    # get antenna positions from file
    file = np.genfromtxt(pathAntennaList, dtype = "str")
    # file[:,0] and file[:,1] are useless (they are simply "AntennaPosition" and "=")
    
    # get the x, y and z positions
    antennaInfo["x"] = file[:,2].astype(float) * 100 # convert to m
    antennaInfo["y"] = file[:,3].astype(float) * 100 # convert to m
    antennaInfo["z"] = file[:,4].astype(float) * 100 # convert to m
    # get the IDs of the antennas
    antennaInfo["ID"] = file[:,5]

    return antennaInfo



def get_antenna_position(pathAntennaList, antenna):
    """
    get the position for one antenna from SIM??????.list
    .list files are structured like "AntennaPosition = x y z name"

    """
    with open(pathAntennaList, "r") as datafile:
        for line in datafile:
            if antenna in line:
                positions = line
    """
    >>> positions
    'AntennaPosition = 39.0 39.0 1000 16\n'
    >>> positions.split(" ")
    ['AntennaPosition', '=', '39.0', '39.0', '1000', '16\n']
    """
    x = np.float32(positions.split(" ")[2]) * 100 # convert to m
    y = np.float32(positions.split(" ")[3]) * 100 # convert to m
    z = np.float32(positions.split(" ")[4]) * 100 # convert to m

    return x, y, z



def read_long(pathLongFile):
    """
    read the .long Corsika output file, which contains the "longitudinal profile"
    more specifically, it contains the energy deposit and particle numbers for different particles
    since the files are set up in two blocks, this function helps read the data from it
    
    this function is mostly taken from corsika_long_parser.py in the coreasutilities module by Felix Schlüter

    TODO: fix hillas_parameter - something's not working yet
    """
    with open(pathLongFile, mode="r") as file:
        # create a temporary file to write the corrected contents
        temp_file = io.StringIO()

        for line in file:
            # use a regex to search for a minus sign that is not part of an exponent
            if search(r"(?<!e)(-)(?=\d)", line):
                # if the minus sign is not part of an exponent, replace it with a space and a minus sign
                line = line.replace("-", " -")
            # write the corrected line to the temporary file
            temp_file.write(line)

        # set the file pointer to the beginning of the temporary file
        temp_file.seek(0)

        # read the contents of the temporary file into a list of strings
        lines = temp_file.readlines()



    n_steps = int(lines[0].rstrip().split()[3])

    # store n table
    n_data_str = io.StringIO()
    n_data_str.writelines(lines[2:(n_steps + 2)])
    n_data_str.seek(0)
    n_data = np.genfromtxt(n_data_str)

    # store dE table
    dE_data_str = io.StringIO()
    dE_data_str.writelines(lines[(n_steps + 4):(2 * n_steps + 4)])
    dE_data_str.seek(0)
    dE_data = np.genfromtxt(dE_data_str)

    # read out hillas fit
    hillas_parameter = []
    # for line in lines:
    #     if bool(search("PARAMETERS", line)):
    #         hillas_parameter = [float(x) for x in line.split()[2:]]
    #     if bool(search("CHI", line)):
    #         hillas_parameter.append(float(line.split()[2]))


    print("The file", pathLongFile, "has been separated into energy deposit and particle distribution.")
    return n_data, dE_data, hillas_parameter
