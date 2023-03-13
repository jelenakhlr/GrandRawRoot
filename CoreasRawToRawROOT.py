import glob
import datetime #to get the unix timestamp
import time #to get the unix timestamp
from grand.io.root_trees import *
import raw_root_trees as RawTrees
from CorsikaInfoFuncs import *



def CoreasRawToRawROOT(path):
  """
  put meaningful comments here

  """
  # output file name hardcoded in the root tree section

  ###############################
  # Part A: load Corsika files  #
  ###############################

  print("Checking", path, "for *.reas and *.inp files (shower info).")
  # TODO: check if the IDs in SIM.reas and RUN.inp match
  # TODO: maybe specify RUN number as input value - or as choice when there is more than one 

  # ********** load SIM.reas **********
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


  # ********** load RUN.inp **********
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
  available_traces = glob.glob(path + "SIM??????_coreas/*.dat")
  print("Found", len(available_traces), "*.dat files (traces).")
  print("*****************************************")
  # in each dat file:
  # time stamp and the north-, west-, and vertical component of the electric field


  ###############################
  # Part B: Generate ROOT Trees #
  ###############################

  #########################################################################################################################
  # Part B.I.i: get the information from Coreas input files
  #########################################################################################################################   
  # from reas file
  CoreCoordinateNorth = read_params(reas_input, "CoreCoordinateNorth") * 100 # convert to m
  CoreCoordinateWest = read_params(reas_input, "CoreCoordinateWest") * 100 # convert to m
  CoreCoordinateVertical = read_params(reas_input, "CoreCoordinateVertical") * 100 # convert to m
  CorePosition = [CoreCoordinateNorth, CoreCoordinateWest, CoreCoordinateVertical]

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
  # fix these values! for now: placeholders
  # TODO: read and store all seeds as a list

  RandomSeed = [1,2,3,4,5,6]
  
  ecuts = [1,2,3,4] 
  # 0: hadrons & nuclei, 1: muons, 2: e-, 3: photons
  GammaEnergyCut    = ecuts[3]
  ElectronEnergyCut = ecuts[2]
  MuonEnergyCut     = ecuts[1] # TODO: check if this exists in Zhaires
  HadronEnergyCut   = ecuts[0]
  NucleonEnergyCut  = ecuts[0]
  MesonEnergyCut    = HadronEnergyCut # mesons are hadronic, so this should be fine

  parallel = [1,2] # COREAS-only
  # PARALLEL = [ECTCUT, ECTMAX, MPIID, FECTOUT]
  # ECTCUT: limit for subshowers
  # ECTMAX: maximum energy for complete shower
  # MPIID: ID for mpi run (ignore for now)
  # T/F flag for extra output file (ignore for now)

  elmflg = [1,2,3] # COREAS-only

  # In Zhaires converter: RelativeThinning, WeightFactor
  # I have:
  thin  = [1,2,3] 
  # THIN = [limit, weight, Rmax]
  thinh = [1,2,3] 
  # THINH = [limit, weight] for hadrons
  
  # !!!! suggestion: 
  ThinLimit  = [thin[0], thinh[0]]
  ThinWeight = [thin[1], thinh[1]]
  ThinRmax   = thin[2]


  # TODO: read keywords with T or F flags as string
  mumult = "T"
  muaddi = "T"
  parout = ["T", "F"]
  longi  = ["T", "5", "T", "T"]
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  # from long file
  # TODO: add params from long file

  EnergyInNeutrinos = 1. # placeholder
  # + energy in all other particles

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

  # TODO: add function for reading logs
  # TODO: add CPU time
  CPUTime = 1.
  # TODO: find injection altitude
  InjectionAltitude = 100.

  ArrayName = "GP13" # TODO: unhardcode this

  # TODO: find xmax
  # SlantXmax
  # XmaxPosition
  # XmaxDistance
  # XmaxAltitude

  ############################################################################################################################
  # Part B.I.ii: Create and fill the RAW Shower Tree
  ############################################################################################################################
  OutputFileName = "Coreas_" + EventName +".root"

  # The tree with the Shower information common to ZHAireS and Coreas
  RawShower = RawTrees.RawShowerTree(OutputFileName)
  # The tree with Coreas-only info
  SimCoreasShower = RawTrees.RawCoreasTree(OutputFileName)
  # TODO: figure out what goes here -> its all leftover info 

  # ********** fill RawShower **********
  RawShower.run_number = RunID
  RawShower.shower_sim = "Coreas"  # TODO:Unhardcode this, add version, etc
  RawShower.event_number = EventID
  RawShower.event_name = EventName
  RawShower.event_date = Date
  RawShower.unix_date = UnixDate

  RawShower.rnd_seed = RandomSeed

  RawShower.energy_in_neutrinos = EnergyInNeutrinos
  RawShower.prim_energy = [Energy]
  RawShower.shower_azimuth = Azimuth
  RawShower.shower_zenith = Zenith
  RawShower.prim_type = [str(Primary)]
  RawShower.prim_inj_alt_shc = [InjectionAltitude]
  RawShower.atmos_model = str(AtmosphericModel)

  RawShower.magnetic_field = np.array([FieldInclination,FieldDeclination,FieldIntensity])
  # RawShower.xmax_grams = SlantXmax
  # RawShower.xmax_pos_shc = XmaxPosition
  # RawShower.xmax_distance = XmaxDistance
  # RawShower.xmax_alt = XmaxAltitude
  RawShower.hadronic_model = HadronicModel
  RawShower.low_energy_model = LowEnergyModel
  RawShower.cpu_time = float(CPUTime)

  # from zhaires converter:
  # RawShower.relative_thinning = RelativeThinning
  # RawShower.weight_factor = WeightFactor  
  # change to:
  RawShower.thinning = ThinLimit # list of floats
  RawShower.weight_factor = ThinWeight # list of floats
  RawShower.thin_Rmax = ThinRmax # float

  RawShower.gamma_energy_cut = GammaEnergyCut
  RawShower.electron_energy_cut = ElectronEnergyCut
  RawShower.muon_energy_cut = MuonEnergyCut
  RawShower.meson_energy_cut = MesonEnergyCut # and hadrons
  RawShower.nucleon_energy_cut = NucleonEnergyCut # same as meson and hadron cut

  RawShower.shower_core_pos = np.array(CorePosition) #do we want np.array(list)?
  #TODO ASAP: ArrayName
  RawShower.array_name = ArrayName
  #TODO ASAP: EventWeight
  #TODO ASAP: TestedCores

####**********************
  # get all info from the long file



  RawShower.fill()
  RawShower.write()


  #########################################################################################################################
  # Part B.II.i: get the information from Coreas output files (i.e. the traces and some extra info)
  #########################################################################################################################   
  
  #****** info from input files: ******

  # placeholders
  # TODO: find these values
  RefractionIndexModel = "model"
  RefractionIndexParameters = [1,2,3] # ? 
  
  TimeWindowMin = TimeLowerBoundary # from reas
  TimeWindowMax = TimeUpperBoundary # from reas
  TimeBinSize   = TimeResolution    # from reas


  #****** load traces ******
  tracefiles = available_traces # from initial file checks

  # get antenna positions from list file (later)
  # get efield for each antenna (later)

  ############################################################################################################################
  # Part B.II.ii: Create and fill the RawEfield Tree
  ############################################################################################################################
 
  #****** fill shower info ******
  RawEfield = RawTrees.RawEfieldTree(OutputFileName)

  RawEfield.run_number = RunID
  RawEfield.event_number = EventID

  RawEfield.efield_sim = "Coreas" # TODO: unhardcode this

  RawEfield.refractivity_model = RefractionIndexModel                                       
  RawEfield.refractivity_model_parameters = RefractionIndexParameters                       
        
  #TODO ASAP: Find how to do this
  #RawEfield.refractivity_profile=Atmostable.T[0].tolist() 
  #RawEfield.refractivity_profile.append(RefrIndex.tolist())
  RawEfield.t_pre = TimeWindowMin
  RawEfield.t_post = TimeWindowMax
  RawEfield.t_bin_size = TimeBinSize

  #****** fill traces ******
 
  RawEfield.du_count = len(tracefiles)


  for antenna in ant_IDs:
    tracefile = glob.glob(path + "SIM??????_coreas/raw_" + str(antenna) + ".dat")

    # load the efield traces for this antenna
    # the files are setup like [timestamp, x polarization, y polarization, z polarization]
    efield = np.loadtxt(tracefile, dtype='f4')
    
    t_0 = np.array(efield[:,0], dtype = np.float32)
    x_polarization = np.array(efield[:,1], dtype = np.float32)
    y_polarization = np.array(efield[:,2], dtype = np.float32)
    z_polarization = np.array(efield[:,3], dtype = np.float32)
    

    # DetectorID = IDs[ant_ID] 
    # is this supposed to be antenna name vs antenna ID?
    
    # TODO: check this
    # in Zhaires converter: AntennaN[ant_ID]
    RawEfield.du_id.append(int(antenna)) 

    RawEfield.t_0.append(t_0)

    # Traces
    RawEfield.trace_x.append(x_polarization)
    RawEfield.trace_y.append(y_polarization)
    RawEfield.trace_z.append(z_polarization)

    # Antenna positions in showers's referential in [m]
    
    pathAntennaList = glob.glob(path + "*.list")[0]
    # TODO: read only the line with the antenna ID

    RawEfield.pos_x.append(ant_position[0])
    RawEfield.pos_y.append(ant_position[1])
    RawEfield.pos_z.append(ant_position[2])


                 
    RawEfield.fill()
    RawEfield.write()



  print("### The event written was " + EventName)

  return EventName



if __name__ == "__main__":
    import sys
    path = sys.argv[1]

    # make sure the last character is a slash
    if (path[-1]!="/"):
        path = path + "/"

    CoreasRawToRawROOT(path)
