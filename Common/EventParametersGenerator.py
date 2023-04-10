

import os
import numpy as np
    
    
#Author: Matias Tueros, with ChatGP3 help for documentation and error handling. it was Mar 24th 2023 in Barracas, Buenos Aires, Argentina
def GenerateEventParametersFile(EventName, Primary, Energy, Zenith, Azimuth, CorePosition, ArrayName, EventWeight=1, EventUnixTime=0, EventUnixNanosecond=0, OutMode="a", TestedPositions="None"):
    '''
    The function generates an event parameters file in a specific format for use in simulation programs. 
    The file includes information such as the event name, primary particle, energy, zenith and azimuth angles, core position, 
    and array name. 
    Optionally, the Event Weight and Unix Time and Unix Nanosecond can be specified
    Additionally, the file can include a list of tested positions that were tried before generating the event. 

    The idea behind having this file is to make it friendly for other simulation programs, and to avoid putting extra things in the Aires/ZHAireS/CoREAS .inp file
    For now, the only really needed parameters are ArrayName and CorePosition, but if we implement an antenna selection this would be the place to put that information.
    Also is the place for other information external to Aires/ZHAireS regarding event generation (i.e. parameters of the core position generation or the antenna selection)
    
    Inputs:
    - EventName (str): The name of the event.
    - Primary (str): The type of primary particle.
    - Energy (float): The energy of the primary particle in GeV.
    - Zenith (float): The zenith angle of the event in degrees.
    - Azimuth (float): The azimuth angle of the event in degrees.
    - CorePosition (tuple): The (x,y,z) coordinates of the core position in meters.
    - ArrayName (str): The name of the array used for the simulation.
    - EventWeight (float,optional): The statistical weight of the event. Default is 1
    - EventUnixTime (unisgned int, optional): The event time in seconds since EPOCH. Default is 0
    - EventUnixTime (unisgned int, optional): The nanoseconds in the event second since EPOCH. Default is 0    
    - OutMode (str, optional): The mode for opening the output file. Default is "a" for append mode.
    - TestedPositions (list of tuples, optional): A list of (x,y,z) tuples representing tested core positions.
    
    Output:
    - None
    
    '''
    try:
        OutputFile=EventName+'.EventParameters' 
        fout= open(OutputFile, OutMode)
        fout.write('#Event Simulation Parameters##################################################################\n')
        fout.write('EventName: {}\n'.format(EventName))
        fout.write('Primary: {}\n'.format(Primary))
        fout.write('PrimaryEnergy: {0:.5} GeV\n'.format(round(float(Energy),5)))
        fout.write('Zenith: {0:.2f} deg\n'.format(round(float(Zenith),5)))
        fout.write('Azimuth: {0:.2f} deg Magnetic\n'.format(round(float(Azimuth),5)))
        fout.write('Core Position: {0:.3f} {1:.3f} {2:.3f}\n'.format(CorePosition[0],CorePosition[1],CorePosition[2]))
        fout.write('ArrayName: {}\n'.format(ArrayName))  
        fout.write('EventWeight: {0:.3f}\n'.format(EventWeight))
        fout.write('EventUnixTime: {0}\n'.format(EventUnixTime))
        fout.write('EventUnixNanosecond: {0}\n'.format(EventUnixNanosecond))                  
        if(TestedPositions!="None"):
            fout.write('#Core positions tested before generating the event ###########################################\n')
            for a in TestedPositions:
                fout.write(str(a)+"\n")  
      
        fout.write('#End of EventParameters####################################################################\n')     
        fout.close()
    except FileNotFoundError:
        print("Error: The file could not be created because the output directory does not exist.")
        return
    except Exception as e:
        print(f"Error: {str(e)}")
        return
    
#Author: ChatGPT3, with prompts from  Matias Tueros, it was Mar 24th 2023 in Barracas, Buenos Aires, Argentina
def GetCorePositionFromParametersFile(filename):
    """
    This function reads the core position from a file generated by the function GenerateEventParametersFile.
    
    Input:
    filename (string): The name of the file to read from.
    
    Output:
    A tuple containing the x, y, and z coordinates of the core position.
    If the file is not found, or if the core position is not found or the format is incorrect, 
    the function returns (0.0, 0.0, 0.0).
    """
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                if "Core Position:" in line:
                    try:
                        x, y, z = line.strip().split(": ")[1].split()
                        return (float(x), float(y), float(z))
                    except (ValueError, IndexError):
                        print("Error: Incorrect format for core position in file:", filename)
                        sys.exit(1)
            # If the loop completes without finding the core position, return (0,0,0)
            return (0.0,0.0,0.0)
    except FileNotFoundError:
        print("Error: File not found:", filename)
        return (0.0,0.0,0.0)
    

#Author: Matias Tueros, it was Mar 24th 2023 in Barracas, Buenos Aires, Argentina
def GetArrayNameFromParametersFile(inp_file):
  try:
    datafile=open(inp_file,'r')
    with open(inp_file, "r") as datafile:
      for line in datafile:
        if 'ArrayName:' in line:
          line = line.lstrip()
          line = line.rstrip()
          stripedline=line.split(':',-1)
          ArrayName=stripedline[1].lstrip()
          return ArrayName
      try:
        ArrayName
      except NameError:
        #logging.error('warning core position not found, defaulting to (0,0,0)')
        print('warning ArrayName not found')
        return "NOT_FOUND"
  except:
    #logging.error("GetCorePositionFromInp:file not found or invalid:"+inp_file)
    print("GetArrayNameFromParametersFile:file not found or invalid:"+inp_file)
    raise
    return -1

    
#Author: ChatGPT v3 with prompts from Matias Tueros, that is very happy learning that he is soon to be out of work.
def GetTestedPositionsFromParametersFile(file_path):
    """
    Reads the core positions tested before generating the event from the given EventParameters file.

    Args:
    file_path (str): The path to the EventParameters file.

    Returns:
    A list of tuples representing the core positions tested before generating the event. Each tuple has three
    elements (x, y, z) representing the position of the core.

    Raises:
    FileNotFoundError: If the EventParameters file is not found.
    ValueError: If the core positions tested section is found but the expected data cannot be read due to formatting
    errors.
    """
    try:
        # Open the file
        with open(file_path, "r") as f:
            # Read the lines from the file
            lines = f.readlines()

        # Check if the core positions tested section exists
        for i, line in enumerate(lines):
            if "Core positions tested before generating the event" in line:
                positions_start = i + 1
                break
        else:
            # Core positions tested section does not exist, return empty list
            print("tested core positions section not found")
            return []

        # Extract the core positions
        positions = []
        for line in lines[positions_start:]:
            # Remove the parentheses and split the line into x, y, z values
            try:
                x, y, z = line.strip()[1:-1].split(" ")
                x = float(x.strip())
                y = float(y.strip())
                z = float(z.strip())
            except (ValueError, IndexError):
                # Skip any lines that do not match the expected format
                print("tested core positions do not match format, expecting (x y z)",line)
                continue
            # Add the core position to the list
            positions.append((x, y, z))

        return positions

    except FileNotFoundError:
        print(f"Error: Could not find the file {file_path}")
        return []
    except:
        print(f"Error: An unexpected error occurred while reading the file {file_path}")
        return []


def GetEventWeightFromParametersFile(filename):
    """
    This function reads the EventWeight from a given Event Parameters file.

    Args:
    filename (str): The name of the Event Parameters file.

    Returns:
    float: The value of the EventWeight parameter.

    Raises:
    FileNotFoundError: If the file cannot be found.
    """
    try:
        with open(filename, 'r') as f:
            for line in f:
                if 'EventWeight:' in line:
                    return float(line.strip().split(':')[1])
            print('EventWeight not found in EventParameters file, defaulting to 1')        
            return 1        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        return None
        

def GetEventUnixTimeFromParametersFile(filename):
    """
    This function reads the EventUnixTime and EventUnixNanosecond from a given Event Parameters file.

    Args:
    filename (str): The name of the Event Parameters file.

    Returns:
    tuple: A tuple containing the Unix time in seconds and nanoseconds.

    Raises:
    FileNotFoundError: If the file cannot be found.
    """
    try:
        with open(filename, 'r') as f:
            unix_time=0
            unix_ns =0
            for line in f:
                if 'EventUnixTime:' in line:
                    unix_time = int(line.strip().split(':')[1])
                elif 'EventUnixNanosecond:' in line:
                    unix_ns = int(line.strip().split(':')[1])
            if( unix_time==0 and unix_ns==0):
              print('Unix time not found in file,or was 0,0...defaulting to 0,0')            
            return (unix_time, unix_ns)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        return None



    
if __name__ == '__main__':

  import numpy as np
  import sys
  
  print (np.size(sys.argv))
  
  if np.size(sys.argv)!=7:
    print ("Arguments = EventName, Primary, Energy, Zenith, Azimuth, ArrayName")
    print ("\n")
    print ("For example: python3 Test Proton 10 20 30 JuanDomingo")
  #
  else:
     
    #Get the parameters  
    EventName = sys.argv[1]
    Primary = sys.argv[2] 
    Energy = float(sys.argv[3]) #in deg
    Zenith = float(sys.argv[4]) #in deg
    Azimuth = float(sys.argv[5]) #in deg
    ArrayName = sys.argv[6]
    #
    print("Going to try creating an EventParameters file")
    print("EventName", EventName)
    print("Primary", Primary)
    print("Energy [GeV]", Energy)
    print("Zenith [Deg]" , Zenith)
    print("Azimuth [Deg]", Azimuth)
    print("ArrayName",ArrayName)
  
    TestedPositions=[]  
    for i in range(0,10):
      TestedPositions.append((np.random.rand(),np.random.rand(),np.random.rand()))
     
    
    EventWeight=3.2
    
    UnixTime=123456789
    UnixNano=987654321
    
    CorePosition=(np.random.rand(),np.random.rand(),np.random.rand())
    print("Going in with Core Position", CorePosition)       
    GenerateEventParametersFile(EventName, Primary, Energy, Zenith, Azimuth, CorePosition, ArrayName, EventWeight, UnixTime, UnixNano, OutMode="w", TestedPositions=TestedPositions )
    
    print("Found Tested positions")
    positions=GetTestedPositionsFromParametersFile(EventName+".EventParameters")
    print(positions)
    
    CorePosition=GetCorePositionFromParametersFile(EventName+".EventParameters")
    print("Found Core Position")
    print(CorePosition)
    
    
    
    EventWeight=GetEventWeightFromParametersFile(EventName+".EventParameters")
    print("Found Weight")
    print(EventWeight)
    
    
    UnixTime, UnixNano = GetEventUnixTimeFromParametersFile(EventName+".EventParameters")
    print("Found Time")
    print(UnixTime,UnixNano)
    
    
    
    