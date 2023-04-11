# Welcome to RawRoot
The RawRoot file format is a step towards common output between CoREAS and ZHAireS. Its based on the GRAND root file format.

The idea is for it to be the starting point for the generation of GRANDRoot files.

## Common
Here we have the files that are common to both simulators

## CoREASRawRoot
Here we have the scripts to produce RawRoot files from CoREAS simulations

To run the script on the provided example event just go to this directory and do

python3 CoreasRawToRawROOT.py 000004/

## ZHAireSRawRoot
Here we have the scripts to produce RawRoot files from ZHAireS simulations

To run the script on the provided example event just go to this directory and do

python3 ZHAireSRawToRawROOT.py ./GP10_192745211400_SD075V standard 0 1  GP10_192745211400_SD075V.root
