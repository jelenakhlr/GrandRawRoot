from re import search
from grand.io.root_trees import *

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

def get_antenna_positions(pathAntennaList):
    """
    get antenna positions from SIM??????.list
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











