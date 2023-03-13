from optparse import OptionParser
import glob
import time
from grand.io.root_trees import *
from CorsikaInfoFuncs import *
import raw_root_trees as RawTrees

# TODO: check if the IDs in SIM.reas and RUN.inp match


def CoreasRawToRawROOT(path):
  """
  put meaningful comments here

  """
  OutputFileName = "Test.root"
  #######################
  # load Corsika files  #
  #######################

  print("Checking", path, "for *.reas and *.inp files (shower info).")

  # ********** load SIM.reas and RUN.inp **********
  # TODO: maybe specify RUN number as input value - or as choice when there is more than one 

  # find reas files
  if glob.glob(path + "SIM??????-*.reas"):
      available_reas_files = glob.glob(path + "SIM??????-*.reas") # these are from parallel runs
  else:
      available_reas_files = glob.glob(path + "SIM??????.reas") # these are from normal runs


  # reas status messages
  if len(available_reas_files) == 0:
    print("[ERROR] No reas file found in this directory. Please check directory and try again.")
    quit()
  elif len(available_reas_files) > 1:
    print("Found", available_reas_files)
    print("[WARNING] More than one reas file found in directory. Only reas file", available_reas_files[0], "will be used.")
    reas_input = available_reas_files[0]
  else:
    print("Found", available_reas_files)
    reas_input = available_reas_files[0]
    print("Converting reas file", reas_input, "to GRANDroot.")
  print("*****************************************")



  # find inp files
  available_inp_files = glob.glob(path + "RUN??????.inp")

  # inp status messages
  if len(available_inp_files) == 0:
    print("[ERROR] No input file found in this directory. Please check directory and try again.")
    quit()
  elif len(available_inp_files) > 1:
    print("Found", available_inp_files)
    print("[WARNING] More than one input file found in directory. Only input file", available_inp_files[0], "will be used.")
    inp_input = available_inp_files[0]
  else:
    print("Found", available_inp_files)
    inp_input = available_inp_files[0]
    print("Converting input file", inp_input, "to GRANDroot.")


  # ********** load traces **********
  print("Checking subdirectories for *.dat files (traces).")
  available_traces = glob.glob(options.dir + "SIM??????_coreas/*.dat")
  print("Found", len(available_traces), "*.dat files (traces).")
  print("*****************************************")
  # in each dat file:
  # time stamp and the north-, west-, and vertical component of the electric field



  # prepare traces as in DataStoringExample.py
  filename = reas_input.split(".reas")[0] + ".root"
  event_count = 1
  adc_traces = []
  traces = []
  for ev in range(event_count):
      adc_traces.append([])
      traces.append([])
      for i, file in enumerate(available_traces):
          adc_traces[-1].append(
              (
                  np.genfromtxt(file)[:,0].astype(np.int16),
                  np.genfromtxt(file)[:,1].astype(np.int16),
                  np.genfromtxt(file)[:,2].astype(np.int16),
                  np.genfromtxt(file)[:,3].astype(np.int16),
              )
          )
          traces[-1].append(
              (
                  (adc_traces[-1][i][0] * 0.9 / 8192).astype(np.float32),
                  (adc_traces[-1][i][1] * 0.9 / 8192).astype(np.float32),
                  (adc_traces[-1][i][2] * 0.9 / 8192).astype(np.float32),
              )
          )



  #######################
  # Generate ROOT Trees #
  #######################

  
  #############################################################################################################################
  # ShowerSimInfo (deals with the details for the simulation). This might be simulator-dependent (CoREAS has different parameters)
  #############################################################################################################################
  
  # The tree with the Shower information common to ZHAireS and Coreas
  RawShower = RawTrees.RawShowerTree(OutputFileName)
  # The tree with Coreas-only info



  #########################################################################################################################
  # Part I: get the information from Coreas
  #########################################################################################################################   
  # from reas file
  CoreCoordinateNorth = read_params(reas_input, "CoreCoordinateNorth") * 100
  CoreCoordinateWest = read_params(reas_input, "CoreCoordinateWest") * 100
  CoreCoordinateVertical = read_params(reas_input, "CoreCoordinateVertical") * 100

  TimeResolution = read_params(reas_input, "TimeResolution")
  AutomaticTimeBoundaries = read_params(reas_input, "AutomaticTimeBoundaries")
  TimeLowerBoundary = read_params(reas_input, "TimeLowerBoundary")
  TimeUpperBoundary = read_params(reas_input, "TimeUpperBoundary")
  ResolutionReductionScale = read_params(reas_input, "ResolutionReductionScale")

  GroundLevelRefractiveIndex = read_params(reas_input, "GroundLevelRefractiveIndex") # refractive index at 0m asl

  RunID = read_params(reas_input, "RunNumber")
  EventID = read_params(reas_input, "EventNumber")
  GPSSecs = read_params(reas_input, "GPSSecs")
  GPSNanoSecs = read_params(reas_input, "GPSNanoSecs")
  RotationAngleForMagfieldDeclination = read_params(reas_input, "RotationAngleForMagfieldDeclination") # in degrees

  Zenith = read_params(reas_input, "ShowerZenithAngle")
  Azimuth = read_params(reas_input, "ShowerAzimuthAngle")

  Energy = read_params(reas_input, "PrimaryParticleEnergy") # in eV
  Primary = read_params(reas_input, "PrimaryParticleType") # as defined in CORSIKA -> TODO: change to PDG system
  DepthOfShowerMaximum = read_params(reas_input, "DepthOfShowerMaximum") # slant depth in g/cm^2
  DistanceOfShowerMaximum = read_params(reas_input, "DistanceOfShowerMaximum") # geometrical distance of shower maximum from core in cm
  MagneticFieldStrength = read_params(reas_input, "MagneticFieldStrength") # in Gauss
  MagneticFieldInclinationAngle = read_params(reas_input, "MagneticFieldInclinationAngle") # in degrees, >0: in northern hemisphere, <0: in southern hemisphere
  GeomagneticAngle = read_params(reas_input, "GeomagneticAngle") # in degrees


  # from inp file
  nshow = read_params(inp_input, "NSHOW") # number of showers - should be 1 for coreas, so maybe we dont need this parameter at all
  ectmap = read_params(inp_input, "ECTMAP")
  maxprt = read_params(inp_input, "MAXPRT")
  radnkg = read_params(inp_input, "RADNKG")

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # TODO: read and store all seeds as a list
  # RandomSeed = read_params(inp_input, "SEED")
  
  # TODO: read keywords with multiple values as string
  # ecuts = 
  # parallel = 
  # elmflg = 
  # thin = 
  # thinh = 
  
  # TODO: read keywords with T or F flags as string
  # MUMULT  T
  # MUADDI  T
  # PAROUT  T  F
  # LONGI   T   5.     T       T
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  # RawROOT addition
  EventName = "Event_" + str(EventID)

  AtmosphericModel = read_atmos(inp_input)
  Date = "2017-04-01" # from ATM file. TODO: unhardcode this
  t1 = time.strptime(Date.strip(),"%d %B, %Y")
  UnixDate = int(time.mktime(t1))

  HadronicModel = "sibyll" #TODO:Unhardcode this
  print("Warning, hard-coded HadronicModel", HadronicModel) 
  LowEnergyModel = "urqmd" #TODO:Unhardcode this
  print("Warning, hard-coded LowEnergyModel", LowEnergyModel)

  # TODO: add CPU time
  # TODO: add function for reading logs

  ############################################################################################################################
  # Part II: Fill RawShower TTree 
  ############################################################################################################################
  RawShower.run_number = RunID
  RawShower.shower_sim = "Coreas"  #TODO:Unhardcode this, add version, etc
  RawShower.event_number = EventID
  RawShower.prim_energy = [Energy] #TODO: test multiple primaries
  RawShower.shower_azimuth = Azimuth
  RawShower.shower_zenith = Zenith

  RawShower.fill()
  RawShower.write()

  print("### The event written was " + EventName)

  return EventName



if __name__ == "__main__":
    import sys
    path = sys.argv[1]

    # make sure the last character is a slash
    if (path[-1]!="/"):
        path = path + "/"

    CoreasRawToRawROOT(path)
