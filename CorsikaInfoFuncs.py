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
    x = positions.split(" ")[2] * 100 # convert to m
    y = positions.split(" ")[3] * 100 # convert to m
    z = positions.split(" ")[4] * 100 # convert to m

    return x, y, z


def read_long(pathLongFile):
    """
    read the .long Corsika output file, which contains the "longitudinal profile"
    more specifically, it contains the energy deposit and particle numbers for different particles
    since the files are set up in two blocks, this function helps read the data from it
    
    WARNING: this reader only works for long files which were NOT created with mpi!

    """
    # get the DAT?????? name of the file
    longName = pathLongFile.split("/")[-1].split(".")[0]

    # open the file
    with open(pathLongFile, mode="r") as file:
        """
        this is most likely the ugliest function i have ever written, but as long
        as it works, i guess it's okay

        there's one space between the columns, but if there's a negative value
        the minus sign is put in the place of the space, which can cause problems
        
        so just to be safe, add spaces before negative signs like this:
        (but not before the minus in e.g. e-02)
        """
        for line in file:
            if search(r"(?<!e)(-)(?=\d)", line):
                # if the minus is actually a negative sign, replace:
                line.replace("-", " -")

        # now read the blocks from the file and save to new separate files
        reader = file.read()
        for i, block in enumerate(reader.split(' LONGITUDINAL')):

            # the block for i=0 is empty, so skip that
            if i == 0:
                pass

            # block 1 will be the particle distribution
            elif i== 1:
                with open(longName + "_particledist.txt", "w") as newfile:
                    newfile.write(block)
                
            # block 2 will be the energy deposit
            elif i==2:
                with open(longName + "_energydep.txt", "w") as newfile:
                    newfile.write(block)

            # block 3 will be parameters for the hillas curve
            elif i==3:
                with open(longName + "_hillasparams.txt", "w") as newfile:
                    newfile.write(block)

            # there should not be any more blocks, but if they are, print them out here
            else:
                print(i, block)

    # read the separate files using numpy, because numpy makes life easier
    # particle dist and energy deposit have columns, the hillas file is different
    # TODO: the hillas params file is not setup with columns, so do this later
    
    return print("The file", longName, "has been separated into energy deposit and particle distribution.")
