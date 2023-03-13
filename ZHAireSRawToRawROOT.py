#!/usr/bin/python
import sys
import os
import glob
import logging
import numpy as np
import datetime #to get the unix timestamp
import time #to get the unix timestamp


logging.basicConfig(level=logging.INFO)
#
#you can get ZHAIRES python from https://github.com/mjtueros/ZHAireS-Python (checkout the Development or DevelopmentLeia branch)
#I use this environment variable to let python know where to find it, but alternatively you just copy the AiresInfoFunctionsGRANDROOT.py file on the same dir you are using this.
#ZHAIRESPYTHON=os.environ["ZHAIRESPYTHON"]
#sys.path.append(ZHAIRESPYTHON)
#sys.path.append("/home/mjtueros/GRAND/ROOTFileFormat/DataChallenge/grandgithub/grand/grand/io")
import AiresInfoFunctionsGRANDROOT as AiresInfo
import ZHAireSInputGenerator as ZHAireSGen #the functions i use from this file should be moved to root_trees_raw, so that we dont need an additional new file. It will be common to Coreas and ZhAireS.


import raw_root_trees as RawTrees

logging.basicConfig(level=logging.INFO)	
logging.getLogger('matplotlib').setLevel(logging.ERROR)



def ZHAiresRawToRawROOT(OutputFileName, RunID, EventID, InputFolder, TaskName="LookForIt", EventName="UseTaskName", NLongitudinal=False, ELongitudinal=False, NlowLongitudinal=False, ElowLongitudinal=False, EdepLongitudinal=False, LateralDistribution=False, EnergyDistribution=False):
    '''
    This routine will read a ZHAireS simulation located in InputFolder and put it in the Desired OutputFileName. 
    RunID is the ID of the run the event is going to be associated with.
    EventID is the ID of the Event is going to be associated with.
    TaskName is what ZHAireS uses to name all the output files. And it generally the same as the "EventName". If you want to change the Name of the event when you store it, you can use the "EventName" optional parameter
    If you dont specify a TaskName, it will look for any .sry file in the directory. There should be only one.
    
    The routine is designed for events simulated with ZHAireS 1.0.30a or later.
    
    It requires to be present in the directory: 
    - 1) A TaskName.sry file, with a "TaskName" inside, that is what ZHAireS uses to name all the output files. If "EventName" is not provided in the function input (maybe you want to override it), it will use TaskName as EventName (recommended).
         Note that Aires will truncate the TaskName if it is too long, and the script will fail. Keep your TaskNames with a reasonable lenght.
         ZHAiresRawToRawROOT will take Energy, Primary, Zenith, Azimuth, etc of the simulation from that file.
         This has the upside that you dont need to keep the idf file (which is some MB) and you dont need to have Aires installed on your system.
         But it has the downside that the values are rounded to 2 decimal places. So if you input Zenith was 78.123 it will be read as 78.12
         
         Since in Aires 19.04.10 there is a python interface that could read idf files, we could get the exact value from there, but for now i will keep it as it is. Just dont produce inputs with more than 2 decimal places!
                  
    - 2) All the a##.trace files produced by ZHAireS with the electric field output (you have to have run your sims with CoREASOutput On)
    - 3) a TaskName.EventParameters file , where the event meta-ZHAireS data is stored (for example the core position used to  generate this particular event, the "ArrayName", the event weight, the event time and nanosecond etc.
    - 4) Optional: the necesary longitudinal tables file. If they dont exist, but the TaskName.idf file is present and  and Aires is installed in the system, AiresInfoFunctions will take care of it.
    -  

    The script output is the shower reference frame: this means a cartesian coordinate system with the shower core at 0,0, GroundAltitude (masl).
    Electric fields are output on the X,Y and Z directions in this coordinate system.

    '''
   
    #TBD: Think about RunID and EventID conventions. (data type, will it just be consecutive numbers or will they code something like the date or type of run -normal,test,calib- etc.)
    #TBD: Do we need a RunID and an EventID at this level? Should it be in the EventParameters file?
    
	#Nice Feature: If EventID is decideded to be just a number, get the next correlative EventID (with respect to the last inside the file) if no EventID is specified
	
    #The function will write two main sections: ShowerSim and EfieldSim. 
    #
    SimShowerInfo=True
    SimEfieldInfo=True 
    
    #We will start by storing the tables Coreas and Zhaires have in common.
    #In the future, i might store additional tables (or other sim info) in a separate tree)
    
    NLongitudinal=False
    
    #########################################################################################################
    #ZHAIRES Sanity Checks
    #########################################################################################################
    #TODO: Handle when InputFolder Does not exist

    #Provide the advertised functionality for sry file
    if TaskName=="LookForIt":
      sryfile=glob.glob(InputFolder+"/*.sry")

      if(len(sryfile)>1):
        logging.critical("there should be one and only one sry file in the input directory!. cannot continue!")
        return -1

      if(len(sryfile)==0):
        logging.critical("there should be one and only one sry file in the input directory!. cannot continue!")
        return -1
    
    else:     
      sryfile=[InputFolder+"/"+TaskName+".sry"]

    #Check Existance of .sry file    
    if  not os.path.isfile(sryfile[0]) :
        logging.critical("Input Summary file not found, {}".format(sryfile))
        return -1
   
    TaskName=AiresInfo.GetTaskNameFromSry(sryfile[0])
   
    EventParametersFile=[InputFolder+"/"+TaskName+".EventParameters"]
    
    if  not os.path.isfile(EventParametersFile[0]) :
        logging.critical("Input EventParametersFile file not found, {}".format(EventParametersFile))
        # return i will not return, in order to be able to handle old sims. I will asign default or dummy values to the required variables
        ArrayName="Unknown"  
    else:        
        ArrayName=ZHAireSGen.GetArrayNameFromParametersFile(EventParametersFile[0])    



    #TODO:idf file is optional in principle, so i dont check for it for now. I should check for the existance of the table files and if anyone fails, check for the existance of the idf file
    idffile=[InputFolder+"/"+EventName+".idf"]

    #provide the advertised functionality for TaskName
    if EventName=="UseTaskName":
      EventName=TaskName
    
    logging.info("###")
    logging.info("###")
    logging.info("### Starting with event "+EventName+" in "+ InputFolder+" to add to "+OutputFileName) 
	
	###########################################################################################################	
    #Root Sanity Checks 
    ###########################################################################################################
    #TODO: Handle when OutputFileName or RootFileHandle is invalid (does not exist, or is not writable, or whatever we decide is invalid (i.e. must already have the RunInfo trees).)
	#TODO: Handle when RunID is invalid: (not a number, not following conventions if any)
	#TODO: Handle when EventID is invalid: it must be unique at least inside the run, (and probabbly unique among all runs in file, depending on conventions)
    	    
    #############################################################################################################################
    # ShowerSimInfo (deals with the details for the simulation). This might be simulator-dependent (CoREAS has different parameters)
    #############################################################################################################################
    if(SimShowerInfo):

        #The tree with the Shower information common to ZHAireS and Coreas       
        RawShower = RawTrees.RawShowerTree(OutputFileName)
        # The tree with ZHAireS only information
        SimZhairesShower = RawTrees.RawZHAireSTree(OutputFileName)

        #########################################################################################################################
        # Part I: get the information from ZHAIRES (for COREAS, its stuff would be here)
        #########################################################################################################################   
        Primary= AiresInfo.GetPrimaryFromSry(sryfile[0],"GRAND")            #Used
        Zenith = AiresInfo.GetZenithAngleFromSry(sryfile[0],"Aires")        #Used
        Azimuth = AiresInfo.GetAzimuthAngleFromSry(sryfile[0],"Aires")      #Used
        Energy = AiresInfo.GetEnergyFromSry(sryfile[0],"Aires")             #Used
        XmaxAltitude, XmaxDistance, XmaxX, XmaxY, XmaxZ = AiresInfo.GetKmXmaxFromSry(sryfile[0])  #Used all
        #Convert to m
        XmaxAltitude= float(XmaxAltitude)*1000.0
        XmaxDistance= float(XmaxDistance)*1000.0
        XmaxPosition= [float(XmaxX)*1000.0, float(XmaxY)*1000.0, float(XmaxZ)*1000.0]
        SlantXmax=AiresInfo.GetSlantXmaxFromSry(sryfile[0])                 #Used        
        InjectionAltitude=AiresInfo.GetInjectionAltitudeFromSry(sryfile[0]) #Used                         
        Date=AiresInfo.GetDateFromSry(sryfile[0])                           #Used
        t1=time.strptime(Date.strip(),"%d/%b/%Y")
        UnixDate = int(time.mktime(t1))                                     #Used
        FieldIntensity,FieldInclination,FieldDeclination=AiresInfo.GetMagneticFieldFromSry(sryfile[0]) #Used
        AtmosphericModel=AiresInfo.GetAtmosphericModelFromSry(sryfile[0])                              #Used
        HadronicModel=AiresInfo.GetHadronicModelFromSry(sryfile[0])                                    #Used
        LowEnergyModel="Aires-Geisha"		                                                           #Used
        print("Warning, hard-coded LowEnergyModel",LowEnergyModel) 
        EnergyInNeutrinos=AiresInfo.GetEnergyFractionInNeutrinosFromSry(sryfile[0])                    #Used
        EnergyInNeutrinos=EnergyInNeutrinos*Energy #to Convert to GeV
        RandomSeed=AiresInfo.GetRandomSeedFromSry(sryfile[0])                                          #Used
        CPUTime=AiresInfo.GetTotalCPUTimeFromSry(sryfile[0],"N/A")                                     #Used
       
        #These might be "run parameters"
        Lat,Long=AiresInfo.GetLatLongFromSry(sryfile[0])                                               # 
        GroundAltitude=AiresInfo.GetGroundAltitudeFromSry(sryfile[0])                                  #
        Site=AiresInfo.GetSiteFromSry(sryfile[0])                                                      #   
        ShowerSimulator=AiresInfo.GetAiresVersionFromSry(sryfile[0])                                   # 
        ShowerSimulator="Aires "+ShowerSimulator                                                       #
        FieldSimulator=AiresInfo.GetZHAireSVersionFromSry(sryfile[0])                                  #
 
  
        RelativeThinning=AiresInfo.GetThinningRelativeEnergyFromSry(sryfile[0])                        #Used
        GammaEnergyCut=AiresInfo.GetGammaEnergyCutFromSry(sryfile[0])                                  #Used
        ElectronEnergyCut=AiresInfo.GetElectronEnergyCutFromSry(sryfile[0])                            #Used
        MuonEnergyCut=AiresInfo.GetMuonEnergyCutFromSry(sryfile[0])                                    #Used
        MesonEnergyCut=AiresInfo.GetMesonEnergyCutFromSry(sryfile[0])                                  #Used
        NucleonEnergyCut=AiresInfo.GetNucleonEnergyCutFromSry(sryfile[0])                              #Used

        #These are ZHAireS specific parameters. Other simulators wont have these parameters, and might have others
        WeightFactor=AiresInfo.GetWeightFactorFromSry(sryfile[0])                                      #Used
 
 
        #MetaZHAires       		        
        #TODO: Document how the core position needs to be stored in the EventParametersFile. 
        #We expet an .EventParameters File that has inside the line:  Core Position: Xcore Ycore Zcore in meters, eg: "Core Position: 2468.927 -4323.117 1998.000" 
 
        if os.path.isfile(EventParametersFile[0]) :
          CorePosition=ZHAireSGen.GetCorePositionFromParametersFile(EventParametersFile[0])
        else:
          CorePosition=(0.0,0.0,0.0)
        #TBD:ASAP should this be removed from the SimData tree? Should the EventParametersFile be included in the root file for future parsing? How are we going to handle this?
   

        ############################################################################################################################# 
        # Part II: Fill RawShower TTree	
        ############################################################################################################################

        RawShower.run_number = RunID
        RawShower.shower_sim = "Aires"  #TODO:Unhardcode this, add version, etc
        RawShower.event_number = EventID
        RawShower.event_name = EventName
        RawShower.event_date = Date
        RawShower.unix_date = UnixDate
        RawShower.rnd_seed = RandomSeed
        RawShower.energy_in_neutrinos = EnergyInNeutrinos    
        RawShower.prim_energy = [Energy] #TODO: test multiple primaries
        RawShower.shower_azimuth = Azimuth
        RawShower.shower_zenith = Zenith
        RawShower.prim_type = [str(Primary)]  #TODO: test multiple primaries
        #TODO ASAP: Add Injection Postion. On neutrino showers this is important (and there can be no Core position)	
        #TODO ASAP prim_injpoint_shc
        RawShower.prim_inj_alt_shc = [InjectionAltitude]
        #TODO ASAP:prim_inj_dir_shc
        RawShower.atmos_model = str(AtmosphericModel) #TODO: Standarize
        #TODO:atmos_model_param  # Atmospheric model parameters: TODO: Think about this. Different models and softwares can have different parameters
        RawShower.magnetic_field = np.array([FieldInclination,FieldDeclination,FieldIntensity])
        RawShower.xmax_grams = SlantXmax
        RawShower.xmax_pos_shc = XmaxPosition
        RawShower.xmax_distance = XmaxDistance                 
        RawShower.xmax_alt = XmaxAltitude
        RawShower.hadronic_model = HadronicModel
        RawShower.low_energy_model = LowEnergyModel
        RawShower.cpu_time = float(CPUTime) 
        
        #ZHAireS/Coreas
        RawShower.relative_thinning = RelativeThinning
        RawShower.weight_factor = WeightFactor  
        RawShower.gamma_energy_cut = GammaEnergyCut
        RawShower.electron_energy_cut = ElectronEnergyCut
        RawShower.muon_energy_cut = MuonEnergyCut
        RawShower.meson_energy_cut = MesonEnergyCut
        RawShower.nucleon_energy_cut = NucleonEnergyCut              
        
        #METAZHAireS
        RawShower.shower_core_pos=np.array(CorePosition) # shower core position TODO: Coordinates?. Undefined for neutrinos.)
        #TODO ASAP: ArrayName
        #TODO ASAP: EventWeight
        #TODO ASAP: TestedCores

        #Fill the tables
        table=AiresInfo.GetLongitudinalTable(InputFolder,1001,Slant=False,Precision="Simple",TaskName=TaskName)               
        RawShower.long_depth.append(np.array(table.T[0], dtype=np.float32))  
        RawShower.long_gammas.append(np.array(table.T[1], dtype=np.float32)) 

        table=AiresInfo.GetLongitudinalTable(InputFolder,1005,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_slantdepth.append(np.array(table.T[0], dtype=np.float32))  
        RawShower.long_eminus.append(np.array(table.T[1], dtype=np.float32))

        table=AiresInfo.GetLongitudinalTable(InputFolder,1006,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_eplus.append(np.array(table.T[1], dtype=np.float32))
        
        table=AiresInfo.GetLongitudinalTable(InputFolder,1008,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_muminus.append(np.array(table.T[1], dtype=np.float32))

        table=AiresInfo.GetLongitudinalTable(InputFolder,1007,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_muplus.append(np.array(table.T[1], dtype=np.float32))
        
        table=AiresInfo.GetLongitudinalTable(InputFolder,1291,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_allch.append(np.array(table.T[1], dtype=np.float32))

        table=AiresInfo.GetLongitudinalTable(InputFolder,1041,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_nuclei.append(np.array(table.T[1], dtype=np.float32))
 
        #TODO: long_hadr is left empty for now, as in zhaires is a combination of several tables...and its rarely used.
        
        table=AiresInfo.GetLongitudinalTable(InputFolder,6796,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_neutrino.append(np.array(table.T[1], dtype=np.float32))
        
        table=AiresInfo.GetLongitudinalTable(InputFolder,7501,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_gamma_cut.append(np.array(table.T[1], dtype=np.float32))        

        table=AiresInfo.GetLongitudinalTable(InputFolder,7705,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_e_cut.append(np.array(table.T[1], dtype=np.float32)) 

        table=AiresInfo.GetLongitudinalTable(InputFolder,7707,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_mu_cut.append(np.array(table.T[1], dtype=np.float32)) 
         
        #TODO: long_hadr_cut is left empty for now, as in zhaires is a combination of several tables...and its rarely used.
                         
        table=AiresInfo.GetLongitudinalTable(InputFolder,7801,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_gamma_ioniz.append(np.array(table.T[1], dtype=np.float32))        

        table=AiresInfo.GetLongitudinalTable(InputFolder,7905,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_e_ioniz.append(np.array(table.T[1], dtype=np.float32)) 

        table=AiresInfo.GetLongitudinalTable(InputFolder,7907,Slant=True,Precision="Simple",TaskName=TaskName)                      
        RawShower.long_mu_ioniz.append(np.array(table.T[1], dtype=np.float32)) 
        
        #TODO: long_hadr_ioniz is left empty for now, as in zhaires is a combination of several tables...and its rarely used.        

        '''
        Atmostable=AiresInfo.GetLongitudinalTable(InputFolder,100,Slant=True,Precision="Simple",TaskName=TaskName)                             
        #RawShower.atmospheric_profile=Atmostable.T[0:3].tolist()
        #TODO:Solve how to store a bidimesional array
        #RawShower.atmospheric_profile.append(np.array(Atmostable.T[0:3],dtype=np.float32))
        '''
        RawShower.fill()
        RawShower.write()


    #############################################################################################################################
    #	SimEfieldInfo
    #############################################################################################################################
    
    #ZHAIRES DEPENDENT
    ending_e = "/a*.trace"
    tracefiles=glob.glob(InputFolder+ending_e)
    #print(tracefiles)
    tracefiles=sorted(tracefiles, key=lambda x:int(x.split(".trace")[0].split("/a")[-1]))
        
    if(SimEfieldInfo and len(tracefiles)>0):

        RawEfield = RawTrees.RawEfieldTree(OutputFileName)                                                                 
        
	    #########################################################################################################################
        # Part I: get the information
        #########################################################################################################################  	
        #Getting all the information i need for	RawEfield
        #
        TimeBinSize=AiresInfo.GetTimeBinFromSry(sryfile[0])
        TimeWindowMin=AiresInfo.GetTimeWindowMinFromSry(sryfile[0])
        TimeWindowMax=AiresInfo.GetTimeWindowMaxFromSry(sryfile[0])

        #make an index of refraction table
        #TODO: Get Refractivity Model parameters from the sry
        RefractionIndexModel="Exponential" #TODO soon: UNHARDCODE THIS
        RefractionIndexParameters=[1.0003250,-0.1218] #TODO soon: UNHARDCODE THIS
        
        print("Warning, hard coded RefractionIndexModel",RefractionIndexModel,RefractionIndexParameters)        
        R0=(RefractionIndexParameters[0]-1.0)*1E6
        #RefrIndex=R0*np.exp(Atmostable.T[0]*RefractionIndexParameters[1]/1000)        
              
        AntennaN,IDs,antx,anty,antz,antt=AiresInfo.GetAntennaInfoFromSry(sryfile[0])


        ############################################################################################################################# 
        # Fill RawEfield part
        ############################################################################################################################ 

        RawEfield.run_number = RunID
        RawEfield.event_number = EventID
        
        
        ############################################################################################################################ 
        # Part II.1: Fill Raw Efield per Event values
        ############################################################################################################################ 
        #Populate what we can
        
        RawEfield.efield_sim="ZHAires"  #TODO: unhardcode this  
        
#       #TODO: Pass this to the event structure 

        RawEfield.refractivity_model = RefractionIndexModel                                       
        RawEfield.refractivity_model_parameters = RefractionIndexParameters                       
         
        #TODO ASAP: Find how to do this
        #RawEfield.refractivity_profile=Atmostable.T[0].tolist() #Atmostable.T[0].tolist()               
        #RawEfield.refractivity_profile.append(RefrIndex.tolist())                                 
        RawEfield.t_pre = TimeWindowMin                                                           
        RawEfield.t_post = TimeWindowMax                                                          
        RawEfield.t_bin_size = TimeBinSize                                                              

        
        ############################################################################################################################ 
        # Part II.2: Fill RawEfield per Antenna
        ############################################################################################################################        

        if(IDs[0]==-1 and antx[0]==-1 and anty[0]==-1 and antz[0]==-1 and antt[0]==-1):
	         logging.critical("hey, no antennas found in event sry "+ str(EventID)+" SimEfield not produced")         #TODO: handle this exeption more elegantly

        else:		

            #convert to 32 bits so it takes less space 
            antx=np.array(antx, dtype=np.float32)
            anty=np.array(anty, dtype=np.float32)
            antz=np.array(antz, dtype=np.float32)
            antt=np.array(antt, dtype=np.float32)
   
            #Important remark. If we need to take into account round earth, then we will need to rotate the electric field components to go to a cartesian frame centered in the array                
            #TODO: check that the number of trace files found is coincidient with the number of antennas found from the sry  
            logging.info("found "+str(len(tracefiles))+" antenna trace files")


            RawEfield.du_count = len(tracefiles)


            for antfile in tracefiles:
                #print("into antenna", antfile)

                ant_number = int(antfile.split('/')[-1].split('.trace')[0].split('a')[-1]) # index in selected antenna list
                                                                                           # TODO soon: Check for this, and handle what hapens if it fails. Maybe there is a more elegant solution                                                                                                                       
                ant_position=(antx[ant_number],anty[ant_number],antz[ant_number])

                efield = np.loadtxt(antfile,dtype='f4') #we read the electric field as a numpy array

                t_0= antt[ant_number]
                
                DetectorID = IDs[ant_number]                                                
                
                print("DetectorID",DetectorID,"AntennaN",AntennaN[ant_number],"ant_number",ant_number,"pos",ant_position,"t0",t_0)

                if(int(AntennaN[ant_number])!=ant_number+1):
                  logging.critical("Warning, check antenna numbers and ids, it seems there is a problem "+str(AntennaN)+" "+str(ant_number+1))
                
                ##########################################################################################################################
                # Part I.1: Convert to GP300 coordinates (here is where customization comes, input specific conventions) (TBD) 
                ##########################################################################################################################
                #TODO: Important remark. If we need to take into account round earth, then we will need to rotate the electric field components to go to a cartesian frame centered in the array
                #this is another reason to have shower coordinates.

                ############################################################################################################################# 
                # Part II: Fill RawEfield	 
                ############################################################################################################################ 
                RawEfield.du_id.append(int(AntennaN[ant_number]))
                #TODO: Save name                
                #RawEfield.du_name.append(DetectorID)
                RawEfield.t_0.append(t_0)            

                # Traces
                RawEfield.trace_x.append(np.array(efield[:,1], dtype=np.float32))
                RawEfield.trace_y.append(np.array(efield[:,2], dtype=np.float32))
                RawEfield.trace_z.append(np.array(efield[:,3], dtype=np.float32))

                # Antenna positions in showers's referential in [m]
                RawEfield.pos_x.append(ant_position[0])
                RawEfield.pos_y.append(ant_position[1])
                RawEfield.pos_z.append(ant_position[2])

                #TODO: Fill p2p and hilbert amplitudes

                         
            #print("Filling RawEfield")
            RawEfield.fill()
            RawEfield.write()
            #print("Wrote RawEfield")

    else:
        logging.critical("no trace files found in "+InputFolder+"Skipping SimEfield") #TODO: handle this exeption more elegantly



	##############################################################################################################################
	# LONGITUDINAL TABLES (not implemented yet, will need to have ZHAIRES installed on your system and the Official version of AiresInfoFunctions).
	##############################################################################################################################

    if(NLongitudinal):
        #the gammas table
        table=AiresInfo.GetLongitudinalTable(InputFolder,1001,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteSlantDepth(HDF5handle, RunID, EventID, table.T[0])
        SimShower.SimShowerWriteNgammas(HDF5handle, RunID, EventID, table.T[1])

        #the eplusminus table, in vertical, to store also the vertical depth
        table=AiresInfo.GetLongitudinalTable(InputFolder,1205,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteVerticalDepth(HDF5handle, RunID, EventID, table.T[0])
        SimShower.SimShowerWriteNeplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the e plus (yes, the positrons)
        table=AiresInfo.GetLongitudinalTable(InputFolder,1006,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNeplus(HDF5handle, RunID, EventID, table.T[1])

        #the mu plus mu minus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1207,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNmuplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the mu plus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1007,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNmuplus(HDF5handle, RunID, EventID, table.T[1])

        #the pi plus pi munus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1211,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNpiplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the pi plus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1011,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNpiplus(HDF5handle, RunID, EventID, table.T[1])

        #and the all charged
        table=AiresInfo.GetLongitudinalTable(InputFolder,1291,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNallcharged(HDF5handle, RunID, EventID, table.T[1])

	##############################################################################################################################
	# Energy LONGITUDINAL TABLES (very important to veryfy the energy balance of the cascade, and to compute the invisible energy)
	##############################################################################################################################
    if(ELongitudinal):
        #the gammas
        table=AiresInfo.GetLongitudinalTable(InputFolder,1501,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEgammas(HDF5handle, RunID, EventID, table.T[1])

        #i call the eplusminus table, in vertical, to store also the vertical depth
        table=AiresInfo.GetLongitudinalTable(InputFolder,1705,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEeplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the mu plus mu minus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1707,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEmuplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the pi plus pi minus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1711,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEpiplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the k plus k minus
        table=AiresInfo.GetLongitudinalTable(InputFolder,1713,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEkplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the neutrons
        table=AiresInfo.GetLongitudinalTable(InputFolder,1521,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEneutrons(HDF5handle, RunID, EventID, table.T[1])

        #the protons
        table=AiresInfo.GetLongitudinalTable(InputFolder,1522,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEprotons(HDF5handle, RunID, EventID, table.T[1])

        #the anti-protons
        table=AiresInfo.GetLongitudinalTable(InputFolder,1523,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEpbar(HDF5handle, RunID, EventID, table.T[1])

        #the nuclei
        table=AiresInfo.GetLongitudinalTable(InputFolder,1541,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEnuclei(HDF5handle, RunID, EventID, table.T[1])

        #the other charged
        table=AiresInfo.GetLongitudinalTable(InputFolder,1591,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEother_charged(HDF5handle, RunID, EventID, table.T[1])

        #the other neutral
        table=AiresInfo.GetLongitudinalTable(InputFolder,1592,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEother_neutral(HDF5handle, RunID, EventID, table.T[1])

        #and the all
        table=AiresInfo.GetLongitudinalTable(InputFolder,1793,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEall(HDF5handle, RunID, EventID, table.T[1])

    ################################################################################################################################
    # NLowEnergy Longitudinal development
    #################################################################################################################################
    if(NlowLongitudinal):
        #the gammas
        table=AiresInfo.GetLongitudinalTable(InputFolder,7001,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNlowgammas(HDF5handle, RunID, EventID, table.T[1])

        #i call the eplusminus table, in vertical, to store also the vertical depth
        table=AiresInfo.GetLongitudinalTable(InputFolder,7005,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNloweplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the positrons (note that they will deposit twice their rest mass!)
        table=AiresInfo.GetLongitudinalTable(InputFolder,7006,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNloweplus(HDF5handle, RunID, EventID, table.T[1])

        #the muons
        table=AiresInfo.GetLongitudinalTable(InputFolder,7207,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNlowmuons(HDF5handle, RunID, EventID, table.T[1])

        #Other Chaged
        table=AiresInfo.GetLongitudinalTable(InputFolder,7091,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNlowother_charged(HDF5handle, RunID, EventID, table.T[1])

        #Other Neutral
        table=AiresInfo.GetLongitudinalTable(InputFolder,7092,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteNlowother_neutral(HDF5handle, RunID, EventID, table.T[1])

    ################################################################################################################################
    # ELowEnergy Longitudinal development
    #################################################################################################################################
    if(ElowLongitudinal):
        #the gammas
        table=AiresInfo.GetLongitudinalTable(InputFolder,7501,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteElowgammas(HDF5handle, RunID, EventID, table.T[1])

        #i call the eplusminus table, in vertical, to store also the vertical depth
        table=AiresInfo.GetLongitudinalTable(InputFolder,7505,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEloweplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the positrons (note that they will deposit twice their rest mass!)
        table=AiresInfo.GetLongitudinalTable(InputFolder,7506,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEloweplus(HDF5handle, RunID, EventID, table.T[1])

        #the muons
        table=AiresInfo.GetLongitudinalTable(InputFolder,7707,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteElowmuons(HDF5handle, RunID, EventID, table.T[1])

        #Other Chaged
        table=AiresInfo.GetLongitudinalTable(InputFolder,7591,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteElowother_charged(HDF5handle, RunID, EventID, table.T[1])

        #Other Neutral
        table=AiresInfo.GetLongitudinalTable(InputFolder,7592,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteElowother_neutral(HDF5handle, RunID, EventID, table.T[1])

    ################################################################################################################################
    # EnergyDeposit Longitudinal development
    #################################################################################################################################
    if(EdepLongitudinal):
        #the gammas
        table=AiresInfo.GetLongitudinalTable(InputFolder,7801,Slant=True,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEdepgammas(HDF5handle, RunID, EventID, table.T[1])

        #i call the eplusminus table, in vertical, to store also the vertical depth
        table=AiresInfo.GetLongitudinalTable(InputFolder,7805,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEdepeplusminus(HDF5handle, RunID, EventID, table.T[1])

        #the positrons (note that they will deposit twice their rest mass!)
        table=AiresInfo.GetLongitudinalTable(InputFolder,7806,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEdepeplus(HDF5handle, RunID, EventID, table.T[1])

        #the muons
        table=AiresInfo.GetLongitudinalTable(InputFolder,7907,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEdepmuons(HDF5handle, RunID, EventID, table.T[1])

        #Other Chaged
        table=AiresInfo.GetLongitudinalTable(InputFolder,7891,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEdepother_charged(HDF5handle, RunID, EventID, table.T[1])

        #Other Neutral
        table=AiresInfo.GetLongitudinalTable(InputFolder,7892,Slant=False,Precision="Simple",TaskName=TaskName)
        SimShower.SimShowerWriteEdepother_neutral(HDF5handle, RunID, EventID, table.T[1])

    ################################################################################################################################
    # Lateral Tables
    #################################################################################################################################
    if(LateralDistribution):
        #the gammas
        table=AiresInfo.GetLateralTable(InputFolder,2001,Density=False,Precision="Simple")
        SimShower.SimShowerWriteLDFradius(HDF5handle, RunID, EventID, table.T[0])
        SimShower.SimShowerWriteLDFgamma(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2205,Density=False,Precision="Simple")
        SimShower.SimShowerWriteLDFeplusminus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2006,Density=False,Precision="Simple")
        SimShower.SimShowerWriteLDFeplus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2207,Density=False,Precision="Simple")
        SimShower.SimShowerWriteLDFmuplusminus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2007,Density=False,Precision="Simple")
        SimShower.SimShowerWriteLDFmuplus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2291,Density=False,Precision="Simple")
        SimShower.SimShowerWriteLDFallcharged(HDF5handle, RunID, EventID, table.T[1])

    ################################################################################################################################
    # Energy Distribution at ground Tables
    #################################################################################################################################
    if(EnergyDistribution):
        #the gammas
        table=AiresInfo.GetLateralTable(InputFolder,2501,Density=False,Precision="Simple")
        SimShower.SimShowerWriteEnergyDist_energy(HDF5handle, RunID, EventID, table.T[0])
        SimShower.SimShowerWriteEnergyDist_gammas(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2705,Density=False,Precision="Simple")
        SimShower.SimShowerWriteEnergyDist_eplusminus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2506,Density=False,Precision="Simple")
        SimShower.SimShowerWriteEnergyDist_eplus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2707,Density=False,Precision="Simple")
        SimShower.SimShowerWriteEnergyDist_muplusminus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2507,Density=False,Precision="Simple")
        SimShower.SimShowerWriteEnergyDist_muplus(HDF5handle, RunID, EventID, table.T[1])

        table=AiresInfo.GetLateralTable(InputFolder,2791,Density=False,Precision="Simple")
        SimShower.SimShowerWriteEnergyDist_allcharged(HDF5handle, RunID, EventID, table.T[1])

    logging.info("### The event written was " + EventName)

    # f.Close()
    # print("****************CLOSED!")
    
    return EventName


#def ZHAiresRawToSimShowerRun(OutputFileName, RunID, EventID, InputFolder):


# TODO: This should probably be part of GRANDRoot.py.?
# Check if the EventID does not already exist in the TTrees
# TODO: Which TTree has all the IDs? Now checking just 2 of them

def CheckIfEventIDIsUnique(EventID, f):
    # Try to get the tree from the file
    try:
        SimShower_tree = f.rawshower
        # This readout should be done with RDataFrame, but it crashes on evt_id :/
        # So doing it the old, ugly way
        # TODO: Ask why it crashes on the ROOT forum and switch to RDataFrame!
        SimShower_tree.Draw("evt_id", "", "goff")
        EventIDs = np.frombuffer(SimShower_tree.GetV1(), dtype=np.float64, count=SimShower_tree.GetSelectedRows()).astype(int)

    # SimShower TTree doesn't exist -> look for SimEfield
    except:
        try:
            SimEfield_tree = f.SimEfield
            SimEfield_tree.Draw("evt_id", "", "goff")
            EventIDs = np.frombuffer(SimEfield_tree.GetV1(), dtype=np.float64, count=SimEfield_tree.GetSelectedRows()).astype(int)

        # No trees - any EventID will do
        except:
            return True

    # If the EventID is already in the trees' EventIDs, return False
    if EventID in EventIDs:
        return False

    return True


if __name__ == '__main__':

	if (len(sys.argv)>6 or len(sys.argv)<6) :
		print("Please point me to a directory with some ZHAires output, and indicate the mode RunID, EventID, EventName and output filename...nothing more, nothing less!")
		print("i.e ZHAiresRawToRawROOT ./MyshowerDir full RunID EventID MyFile.root")
		mode="exit"

	elif len(sys.argv)==6 :
		InputFolder=sys.argv[1]
		mode=sys.argv[2]
		RunID=int(sys.argv[3])
		EventID=int(sys.argv[4])
		OutputFileName=sys.argv[5]

	if(mode=="standard"): 
		ZHAiresRawToRawROOT(OutputFileName,RunID,EventID, InputFolder)

	elif(mode=="full"):

		ZHAiresRawToRawROOT(OutputFileName,RunID,EventID, InputFolder, SimEfieldInfo=True, NLongitudinal=True, ELongitudinal=True, NlowLongitudinal=True, ElowLongitudinal=True, EdepLongitudinal=True, LateralDistribution=True, EnergyDistribution=True)

	elif(mode=="minimal"):
	
		ZHAiresRawToRawROOT(OutputFileName,RunID,EventID, InputFolder, SimEfieldInfo=True, NLongitudinal=False, ELongitudinal=False, NlowLongitudinal=False, ElowLongitudinal=False, EdepLongitudinal=False, LateralDistribution=False, EnergyDistribution=False)

	else:

		print("please enter one of these modes: standard, full or minimal")
		
		
"""
Notes:
	1) Not all SimShower_EventInfo are initialised, for example atmos_model_param
"""
 

