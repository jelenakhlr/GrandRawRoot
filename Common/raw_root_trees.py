"""
Read/Write python interface to GRAND data (real and simulated) stored in Cern ROOT TTrees.

This is the interface for accessing GRAND ROOT TTrees that do not require the user (reader/writer of the TTrees) to have any knowledge of ROOT. It also hides the internals from the data generator, so that the changes in the format are not concerning the user.
"""

from logging import getLogger
import sys
import datetime
import os

import ROOT
import numpy as np
import glob

from collections import defaultdict

# This import changes in Python 3.10
if sys.version_info.major >= 3 and sys.version_info.minor < 10:
    from collections import MutableSequence
else:
    from collections.abc import MutableSequence
        
from dataclasses import dataclass, field    

thismodule = sys.modules[__name__]

from root_trees import *

###########################################################################################################################################################################################################
#
# RawShowerTree
#
##########################################################################################################################################################################################################


@dataclass
## The class for storing a shower simulation-only data for each event
class RawShowerTree(MotherEventTree):
    """The class for storing a shower simulation-only data for each event"""

    _type: str = "rawshower"

    _tree_name: str = "trawshower"
    
    ## Name and version of the shower simulator
    _shower_sim: StdString = StdString("")

    ### Event name (the task name, can be usefull to track the original simulation)
    _event_name: StdString = StdString("")

    ### Event Date  (used to define the atmosphere and/or the magnetic field)
    _event_date: StdString = StdString("")
    
    ### Unix time of this event date
    _unix_date: np.ndarray = field(default_factory=lambda: np.zeros(1, np.uint32))
   
    ### Random seed
    _rnd_seed: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))
    
    ### Energy in neutrinos generated in the shower (GeV). Useful for invisible energy computation
    _energy_in_neutrinos: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))
    
    ### Primary energy (GeV) 
    _energy_primary: StdVectorList = field(default_factory=lambda: StdVectorList("float"))
    
    ### Shower azimuth (deg, CR convention)
    _azimuth: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))

    ### Shower zenith  (deg, CR convention)
    _zenith: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))
    
    ### Primary particle type (PDG)
    _primary_type: StdVectorList = field(default_factory=lambda: StdVectorList("string"))

    # Primary injection point [m] in Shower coordinates
    _prim_injpoint_shc: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))

    ### Primary injection altitude [m] in Shower Coordinates
    _prim_inj_alt_shc: StdVectorList = field(default_factory=lambda: StdVectorList("float"))

    # primary injection direction in Shower Coordinates
    _prim_inj_dir_shc: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))

    ### Atmospheric model name TODO:standardize
    _atmos_model: StdString = StdString("")

    # Atmospheric model parameters: TODO: Think about this. Different models and softwares can have different parameters
    _atmos_model_param: np.ndarray = field(default_factory=lambda: np.zeros(3, np.float32))
    
    # Table of air density [g/cm3] and vertical depth [g/cm2] versus altitude [m]
    _atmos_altitude: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    _atmos_density: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    _atmos_depth: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))

        
    ### Magnetic field parameters: Inclination, Declination, Fmodulus.: In shower coordinates. Declination
    #The Earth’s magnetic field, B, is described by its strength, Fmodulus = ∥B∥; its inclination, I, defined
    # as the angle between the local horizontal plane and the field vector; and its declination, D, defined
    # as the angle between the horizontal component of B, H, and the geographical North (direction of
    # the local meridian). The angle I is positive when B points downwards and D is positive when H is 
    # inclined towards the East.    
    _magnetic_field: np.ndarray = field(default_factory=lambda: np.zeros(3, np.float32))

    ### Shower Xmax depth  (g/cm2 along the shower axis)
    _xmax_grams: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))
    
    ### Shower Xmax position in shower coordinates [m]
    _xmax_pos_shc: np.ndarray = field(default_factory=lambda: np.zeros(3, np.float64))
    
    ### Distance of Xmax  [m] to the ground
    _xmax_distance: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))
    
    ### Altitude of Xmax  [m]. Its important for the computation of the index of refraction at maximum, and of the cherenkov cone
    _xmax_alt: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    ### high energy hadronic model (and version) used TODO: standarize
    _hadronic_model: StdString = StdString("")
    
    ### low energy model (and version) used TODO: standarize
    _low_energy_model: StdString = StdString("")
    
    ### Time it took for the simulation of the cascade (s). In the case shower and radio are simulated together, use TotalTime/(nant-1) as an approximation
    _cpu_time: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))


    #### ZHAireS/Coreas
    # * THINNING *
    # Thinning energy, relative to primary energy
    # this is EFRCTHN in Coreas (the 0th THIN value)
    _relative_thinning: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    # Main weight factor parameter. This has different meaning for Coreas and Zhaires
    # this is WMAX in Coreas (the 1st THIN value) - Weight limit for thinning
    _maximum_weight: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    # this is THINRAT in Coreas (the 0th THINH value) - hadrons
    _hadronic_thinning: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    # this is THINRAT in Coreas (the 1st THINH value)
    _hadronic_thinning_weight: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    # this is RMAX in Coreas (the 2nd THIN value)
    # Maximum radius (in cm) at observation level within which all particles are subject to inner radius thinning
    _rmax: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    # * CUTS *
    #gamma energy cut (GeV)
    _gamma_energy_cut: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    #electron/positron energy cut (GeV)
    _electron_energy_cut: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    #muons energy cut (GeV)
    _muon_energy_cut: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    #mesons energy cut (GeV)
    _meson_energy_cut: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))

    #nucleons energy cut (GeV)
    _nucleon_energy_cut: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float64))


    ###META ZHAireS/Coreas

    ### Core position with respect to the antenna array (undefined for neutrinos)
    _shower_core_pos: np.ndarray = field(default_factory=lambda: np.zeros(3, np.float32))
    
 
    ### Longitudinal Pofiles (those compatible between Coreas/ZHAires)
    
    ## Longitudinal Profile of vertical depth (g/cm2)
    _long_depth: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Longitudinal Profile of slant depth (g/cm2)
    _long_slantdepth: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Longitudinal Profile of Number of Gammas      
    _long_pd_gammas: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Longitudinal Profile of Number of e+
    _long_pd_eplus: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Longitudinal Profile of Number of e-
    _long_pd_eminus: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>")) 
    ## Longitudinal Profile of Number of mu+
    _long_pd_muplus: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Longitudinal Profile of Number of mu-
    _long_pd_muminus: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))      
    ## Longitudinal Profile of Number of All charged particles
    _long_pd_allch: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))   
    ## Longitudinal Profile of Number of Nuclei
    _long_pd_nuclei: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Longitudinal Profile of Number of Hadrons
    _long_pd_hadr:StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))

    ## Longitudinal Profile of Energy of created neutrinos (GeV)
    _long_ed_neutrino: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))           


    ## Longitudinal Profile of low energy gammas (GeV)
    _long_ed_gamma_cut: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))           
    ## Longitudinal Profile of low energy e+/e- (GeV)
    _long_ed_e_cut: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))           
    ## Longitudinal Profile of low energy mu+/mu- (GeV)
    _long_ed_mu_cut: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))           
    ## Longitudinal Profile of low energy hadrons (GeV)
    _long_ed_hadr_cut: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    
    ## Longitudinal Profile of energy deposit by gammas (GeV)
    _long_ed_gamma_ioniz: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))                          
    ## Longitudinal Profile of energy deposit by e+/e-  (GeV)
    _long_ed_e_ioniz: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))           
    ## Longitudinal Profile of energy deposit by muons  (GeV)
    _long_ed_mu_ioniz: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))           
    ## Longitudinal Profile of energy deposit by hadrons (GeV)
    _long_ed_hadr_ioniz: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))     
 

    @property
    def shower_sim(self):
         return str(self._shower_sim)
    
    @shower_sim.setter
    def shower_sim(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(f"Incorrect type for site {type(value)}. Either a string or a ROOT.std.string is required.")
    
        self._shower_sim.string.assign(value)

    @property
    def long_depth(self):
        """Longitudinal profile depth (g/cm2)"""
        return self._long_depth

    @long_depth.setter
    def long_depth(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_depth.clear()
            self._long_depth += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_depth._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_depth {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )            
        

    @property
    def long_slantdepth(self):
        """Longitudinal profile of slant depth (g/cm2)"""
        return self._long_slantdepth

    @long_slantdepth.setter
    def long_slantdepth(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_slantdepth.clear()
            self._long_slantdepth += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_slantdepth._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_slantdepth {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )
 
    @property
    def long_pd_gammas(self):
        """Longitudinal profile of gammas"""
        return self._long_pd_gammas

    @long_pd_gammas.setter
    def long_pd_gammas(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_gammas.clear()
            self._long_pd_gammas += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_gammas._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_gammas {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )    

    @property
    def long_pd_eplus(self):
        """Longitudinal profile of positrons"""
        return self._long_pd_eplus

    @long_pd_eplus.setter
    def long_pd_eplus(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_eplus.clear()
            self._long_pd_eplus += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_eplus._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_eplus {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )    

    @property
    def long_pd_eminus(self):
        """Longitudinal profile of electrons"""
        return self._long_pd_eminus

    @long_pd_eminus.setter
    def long_pd_eminus(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_eminus.clear()
            self._long_pd_eminus += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_eminus._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_eminus {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )    

    @property
    def long_pd_muplus(self):
        """Longitudinal profile of positrons"""
        return self._long_pd_muplus

    @long_pd_muplus.setter
    def long_pd_muplus(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_muplus.clear()
            self._long_pd_muplus += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_muplus._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_muplus {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )    

    @property
    def long_pd_muminus(self):
        """Longitudinal profile of electrons"""
        return self._long_pd_muminus

    @long_pd_muminus.setter
    def long_pd_muminus(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_muminus.clear()
            self._long_pd_muminus += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_muminus._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_muminus {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )  

    @property
    def long_pd_allch(self):
        """Longitudinal profile of all charged particles"""
        return self._long_pd_allch

    @long_pd_allch.setter
    def long_pd_allch(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_allch.clear()
            self._long_pd_allch += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_allch._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_allch {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )  

    @property
    def long_pd_nuclei(self):
        """Longitudinal profile of nuclei"""
        return self._long_pd_nuclei
        
    @long_pd_nuclei.setter
    def long_pd_nuclei(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_nuclei.clear()
            self._long_pd_nuclei += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_nuclei._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_nuclei {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )  


    @property
    def long_pd_hadr(self):
        """Longitudinal profile of hadrons"""
        return self._long_pd_hadr
        
    @long_pd_hadr.setter
    def long_pd_hadr(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_pd_hadr.clear()
            self._long_pd_hadr += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_pd_hadr._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_pd_hadr {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )  

    @property
    def long_ed_neutrino(self):
        """Longitudinal profile of created neutrinos"""
        return self._long_ed_neutrino
        
    @long_ed_neutrino.setter
    def long_ed_neutrino(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_neutrino.clear()
            self._long_ed_neutrino += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_neutrino._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_neutrino {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )  

    @property
    def long_ed_gamma_cut(self):
        """Longitudinal profile of low energy gammas"""
        return self._long_ed_gamma_cut

    @long_ed_gamma_cut.setter
    def long_ed_gamma_cut(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_gamma_cut.clear()
            self._long_ed_gamma_cut += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_gamma_cut._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_gamma_cut {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )
                
    @property
    def long_ed_gamma_ioniz(self):
        """Longitudinal profile of gamma energy deposit"""
        return self._long_ed_gamma_ioniz

    @long_ed_gamma_ioniz.setter
    def long_ed_gamma_ioniz(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_gamma_ioniz.clear()
            self._long_ed_gamma_ioniz += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_gamma_ioniz._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_gamma_ioniz {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )    
            
    @property
    def long_ed_e_cut(self):
        """Longitudinal profile of low energy e+/e-"""
        return self._long_ed_e_cut

    @long_ed_e_cut.setter
    def long_ed_e_cut(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_e_cut.clear()
            self._long_ed_e_cut += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_e_cut._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_e_cut {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )
                
    @property
    def long_ed_e_ioniz(self):
        """Longitudinal profile of energy deposit by e+/e-"""
        return self._long_ed_e_ioniz

    @long_ed_e_ioniz.setter
    def long_ed_e_ioniz(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_e_ioniz.clear()
            self._long_ed_e_ioniz += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_e_ioniz._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_e_ioniz {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )  

    @property
    def long_ed_mu_cut(self):
        """Longitudinal profile of low energy muons"""
        return self._long_ed_mu_cut

    @long_ed_mu_cut.setter
    def long_ed_mu_cut(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_mu_cut.clear()
            self._long_ed_mu_cut += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_mu_cut._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_mu_cut {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )
                
    @property
    def long_ed_mu_ioniz(self):
        """Longitudinal profile of muon energy deposit"""
        return self._long_ed_mu_ioniz

    @long_ed_mu_ioniz.setter
    def long_ed_mu_ioniz(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_mu_ioniz.clear()
            self._long_ed_mu_ioniz += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_mu_ioniz._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_mu_ioniz {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )
            
    @property
    def long_ed_hadr_cut(self):
        """Longitudinal profile of low energy hadrons"""
        return self._long_ed_hadr_cut

    @long_ed_hadr_cut.setter
    def long_ed_hadr_cut(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_hadr_cut.clear()
            self._long_ed_hadr_cut += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_hadr_cut._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_hadr_cut {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )
                
    @property
    def long_ed_hadr_ioniz(self):
        """Longitudinal profile of hadrons energy deposit"""
        return self._long_ed_hadr_ioniz

    @long_ed_hadr_ioniz.setter
    def long_ed_hadr_ioniz(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._long_ed_hadr_ioniz.clear()
            self._long_ed_hadr_ioniz += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._long_ed_hadr_ioniz._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _long_ed_hadr_ioniz {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )         
  
    @property
    def relative_thinning(self):
        """Thinning energy, relative to primary energy"""
        return self._relative_thinning[0]

    @relative_thinning.setter
    def relative_thinning(self, value: np.float64) -> None:
        self._relative_thinning[0] = value 
 

    @property
    def maximum_weight(self):
        """Weight limit for thinning"""
        return self._maximum_weight[0]

    @maximum_weight.setter
    def maximum_weight(self, value: np.float64) -> None:
        self._maximum_weight[0] = value
 

    @property
    def hadronic_thinning(self):
        """hadronic thinning ratio"""
        return self._hadronic_thinning[0]

    @hadronic_thinning.setter
    def hadronic_thinning(self, value: np.float64) -> None:
        self._hadronic_thinning[0] = value


    @property
    def hadronic_thinning_weight(self):
        """hadronic thinning weight ratio"""
        return self._hadronic_thinning_weight[0]

    @hadronic_thinning_weight.setter
    def hadronic_thinning_weight(self, value: np.float64) -> None:
        self._hadronic_thinning_weight[0] = value

    
    @property
    def rmax(self):
        """Maximum radius (in cm) at observation level within which all particles are subject to inner radius thinning"""
        return self._rmax[0]

    @rmax.setter
    def rmax(self, value: np.float64) -> None:
        self._rmax[0] = value



    @property
    def gamma_energy_cut(self):
        """gamma energy cut (GeV)"""
        return self._gamma_energy_cut[0]

    @gamma_energy_cut.setter
    def gamma_energy_cut(self, value: np.float64) -> None:
        self._gamma_energy_cut[0] = value  
      

    @property
    def electron_energy_cut(self):
        """electron energy cut (GeV)"""
        return self._electron_energy_cut[0]

    @electron_energy_cut.setter
    def electron_energy_cut(self, value: np.float64) -> None:
        self._electron_energy_cut[0] = value 

    @property
    def muon_energy_cut(self):
        """muon energy cut (GeV)"""
        return self._muon_energy_cut[0]

    @muon_energy_cut.setter
    def muon_energy_cut(self, value: np.float64) -> None:
        self._muon_energy_cut[0] = value 

    @property
    def meson_energy_cut(self):
        """meson energy cut (GeV)"""
        return self._meson_energy_cut[0]

    @meson_energy_cut.setter
    def meson_energy_cut(self, value: np.float64) -> None:
        self._meson_energy_cut[0] = value 

    @property
    def nucleon_energy_cut(self):
        """nucleon energy cut (GeV)"""
        return self._nucleon_energy_cut[0]

    @nucleon_energy_cut.setter
    def nucleon_energy_cut(self, value: np.float64) -> None:
        self._nucleon_energy_cut[0] = value 


    @property
    def event_name(self):
        """Event name"""
        return str(self._event_name)

    @event_name.setter
    def event_name(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for event_name {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._event_name.string.assign(value)

    @property
    def event_date(self):
        """Event Date"""
        return str(self._date)

    @event_date.setter
    def event_date(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for date {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._event_date.string.assign(value)

    @property
    def rnd_seed(self):
        """Random seed"""
        return self._rnd_seed[0]

    @rnd_seed.setter
    def rnd_seed(self, value):
        self._rnd_seed[0] = value

    @property
    def energy_in_neutrinos(self):
        """Energy in neutrinos generated in the shower (GeV). Usefull for invisible energy"""
        return self._energy_in_neutrinos[0]

    @energy_in_neutrinos.setter
    def energy_in_neutrinos(self, value):
        self._energy_in_neutrinos[0] = value

    @property
    def energy_primary(self):
        """Primary energy (GeV) TODO: Check unit conventions. # LWP: Multiple primaries? I guess, variable count. Thus variable size array or a std::vector"""
        return self._energy_primary

    @energy_primary.setter
    def energy_primary(self, value):
        # A list of strings was given
        if isinstance(value, list):
            # Clear the vector before setting
            self._energy_primary.clear()
            self._energy_primary += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("float")):
            self._energy_primary._vector = value
        else:
            raise ValueError(
                f"Incorrect type for energy_primary {type(value)}. Either a list or a ROOT.vector of floats required."
            )

    @property
    def azimuth(self):
        """Shower azimuth TODO: Discuss coordinates Cosmic ray convention is bad for neutrinos, but neurtino convention is problematic for round earth. Also, geoid vs sphere problem"""
        return self._azimuth[0]

    @azimuth.setter
    def azimuth(self, value):
        self._azimuth[0] = value

    @property
    def zenith(self):
        """Shower zenith TODO: Discuss coordinates Cosmic ray convention is bad for neutrinos, but neurtino convention is problematic for round earth"""
        return self._zenith[0]

    @zenith.setter
    def zenith(self, value):
        self._zenith[0] = value

    @property
    def primary_type(self):
        """Primary particle type TODO: standarize (PDG?)"""
        return self._primary_type

    @primary_type.setter
    def primary_type(self, value):
        # A list of strings was given
        if isinstance(value, list):
            # Clear the vector before setting
            self._primary_type.clear()
            self._primary_type += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("string")):
            self._primary_type._vector = value
        else:
            raise ValueError(
                f"Incorrect type for primary_type {type(value)}. Either a list or a ROOT.vector of strings required."
            )

    @property
    def prim_injpoint_shc(self):
        """Primary injection point in Shower coordinates"""
        return np.array(self._prim_injpoint_shc)

    @prim_injpoint_shc.setter
    def prim_injpoint_shc(self, value):
        set_vector_of_vectors(value, "vector<float>", self._prim_injpoint_shc, "prim_injpoint_shc")

    @property
    def prim_inj_alt_shc(self):
        """Primary injection altitude in Shower Coordinates"""
        return self._prim_inj_alt_shc

    @prim_inj_alt_shc.setter
    def prim_inj_alt_shc(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._prim_inj_alt_shc.clear()
            self._prim_inj_alt_shc += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("float")):
            self._prim_inj_alt_shc._vector = value
        else:
            raise ValueError(
                f"Incorrect type for prim_inj_alt_shc {type(value)}. Either a list, an array or a ROOT.vector of floats required."
            )

    @property
    def prim_inj_dir_shc(self):
        """primary injection direction in Shower Coordinates"""
        return np.array(self._prim_inj_dir_shc)

    @prim_inj_dir_shc.setter
    def prim_inj_dir_shc(self, value):
        set_vector_of_vectors(value, "vector<float>", self._prim_inj_dir_shc, "prim_inj_dir_shc")

    @property
    def atmos_model(self):
        """Atmospheric model name TODO:standarize"""
        return str(self._atmos_model)

    @atmos_model.setter
    def atmos_model(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for site {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._atmos_model.string.assign(value)

    @property
    def atmos_model_param(self):
        """Atmospheric model parameters: TODO: Think about this. Different models and softwares can have different parameters"""
        return np.array(self._atmos_model_param)

    @atmos_model_param.setter
    def atmos_model_param(self, value):
        self._atmos_model_param = np.array(value).astype(np.float32)
        self._tree.SetBranchAddress("atmos_model_param", self._atmos_model_param)

    @property
    def atmos_altitude(self):
        """height above sea level in meters, for the atmos_density and atmos_depth table"""
        return self._atmos_altitude


    @atmos_altitude.setter
    def atmos_altitude(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._atmos_altitude.clear()
            self._atmos_altitude += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._atmos_altitude._vector = value
        else:
            raise ValueError(
                f"Incorrect type for atmos_altitude {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )            
        
    @property
    def atmos_density(self):
        """Table of air density [g/cm3]"""
        return self._atmos_density
 
    @atmos_density.setter
    def atmos_density(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._atmos_density.clear()
            self._atmos_density += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._atmos_density._vector = value
        else:
            raise ValueError(
                f"Incorrect type for atmos_density {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            ) 
 
        
    @property
    def atmos_depth(self):
        """Table of vertical depth [g/cm2]"""
        return self._atmos_depth
        
    @atmos_depth.setter
    def atmos_depth(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._atmos_depth.clear()
            self._atmos_depth += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._atmos_depth._vector = value
        else:
            raise ValueError(
                f"Incorrect type for atmos_depth {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )




    @property
    def magnetic_field(self):
        """Magnetic field parameters: Inclination, Declination, modulus. TODO: Standarize. Check units. Think about coordinates. Shower coordinates make sense."""
        return np.array(self._magnetic_field)

    @magnetic_field.setter
    def magnetic_field(self, value):
        self._magnetic_field = np.array(value).astype(np.float32)
        self._tree.SetBranchAddress("magnetic_field", self._magnetic_field)

    @property
    def xmax_grams(self):
        """Shower Xmax depth (g/cm2 along the shower axis)"""
        return self._xmax_grams[0]

    @xmax_grams.setter
    def xmax_grams(self, value):
        self._xmax_grams[0] = value

    @property
    def xmax_pos_shc(self):
        """Shower Xmax position in shower coordinates"""
        return np.array(self._xmax_pos_shc)

    @xmax_pos_shc.setter
    def xmax_pos_shc(self, value):
        self._xmax_pos_shc = np.array(value).astype(np.float64)
        self._tree.SetBranchAddress("xmax_pos_shc", self._xmax_pos_shc)

    @property
    def xmax_distance(self):
        """Distance of Xmax [m]"""
        return self._xmax_distance[0]

    @xmax_distance.setter
    def xmax_distance(self, value):
        self._xmax_distance[0] = value

    @property
    def xmax_alt(self):
        """Altitude of Xmax (m, in the shower simulation earth. Its important for the index of refraction )"""
        return self._xmax_alt[0]

    @xmax_alt.setter
    def xmax_alt(self, value):
        self._xmax_alt[0] = value

    @property
    def hadronic_model(self):
        """High energy hadronic model (and version) used TODO: standarize"""
        return str(self._hadronic_model)

    @hadronic_model.setter
    def hadronic_model(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for site {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._hadronic_model.string.assign(value)

    @property
    def low_energy_model(self):
        """High energy model (and version) used TODO: standarize"""
        return str(self._low_energy_model)

    @low_energy_model.setter
    def low_energy_model(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for site {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._low_energy_model.string.assign(value)

    @property
    def cpu_time(self):
        """Time it took for the shower + efield simulation."""
        return np.array(self._cpu_time)

    @cpu_time.setter
    def cpu_time(self, value):
        self._cpu_time = np.array(value).astype(np.float32)
        self._tree.SetBranchAddress("cpu_time", self._cpu_time)
        
    @property
    def shower_core_pos(self):
        """Shower core position"""
        return np.array(self._shower_core_pos)

    @shower_core_pos.setter
    def shower_core_pos(self, value):
        self._shower_core_pos = np.array(value).astype(np.float32)
        self._tree.SetBranchAddress("shower_core_pos", self._shower_core_pos)        
        

    @property
    def unix_date(self):
        """The date of the event in seconds since epoch"""
        return self._unix_date[0]

    @unix_date.setter
    def unix_date(self, val: np.uint32) -> None:
        self._unix_date[0] = val



#####################################################################################################################################################################################################
#
# RawEfieldTree
# 
#####################################################################################################################################################################################################
        
@dataclass
## The class for storing Efield simulation-only data common for each event
class RawEfieldTree(MotherEventTree):
    """The class for storing Efield simulation-only data common for each event"""

    _type: str = "rawefield"

    _tree_name: str = "trawefield"

    #Per Event Things
    ## Name and version of the electric field simulator
    _efield_sim: StdString = StdString("")

    ## Name of the atmospheric index of refraction model
    _refractivity_model: StdString = StdString("")
    _refractivity_model_parameters: StdVectorList = field(default_factory=lambda: StdVectorList("double"))    
    _atmos_refractivity: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    
    
    ## The antenna time window is defined around a t0 that changes with the antenna, starts on t0+t_pre (thus t_pre is usually negative) and ends on t0+post
    _t_pre: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))
    _t_post: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))
    _t_bin_size: np.ndarray = field(default_factory=lambda: np.zeros(1, np.float32))    
    
    #Per antenna things
    _du_id: StdVectorList = field(default_factory=lambda: StdVectorList("int"))  # Detector ID
    _du_name: StdVectorList = field(default_factory=lambda: StdVectorList("string"))  # Detector Name
    ## Number of detector units in the event - basically the antennas count
    _du_count: np.ndarray = field(default_factory=lambda: np.zeros(1, np.uint32))


        
    _t_0: StdVectorList = field(default_factory=lambda: StdVectorList("float"))  # Time window t0
    _p2p: StdVectorList = field(default_factory=lambda: StdVectorList("float"))  # peak 2 peak amplitudes (x,y,z,modulus)

    ## X position in shower referential
    _du_x: StdVectorList = field(default_factory=lambda: StdVectorList("float"))
    ## Y position in shower referential
    _du_y: StdVectorList = field(default_factory=lambda: StdVectorList("float"))
    ## Z position in shower referential
    _du_z: StdVectorList = field(default_factory=lambda: StdVectorList("float"))    
    
    ## Efield trace in X direction
    _trace_x: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Efield trace in Y direction
    _trace_y: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))
    ## Efield trace in Z direction
    _trace_z: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))

    ## Efield trace in X,Y,Z direction
    # value = np.concatenate((self._du_id, self._trace_x, self._trace_y, self._trace_z, self.t_bin_size))
    _trace: StdVectorList = field(default_factory=lambda: StdVectorList("vector<float>"))

    @property
    def du_count(self):
        """Number of detector units in the event - basically the antennas count"""
        return self._du_count[0]

    @du_count.setter
    def du_count(self, value: np.uint32) -> None:
        self._du_count[0] = value

    @property
    def efield_sim(self):
         return str(self._efield_sim)
    
    @efield_sim.setter
    def efield_sim(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(f"Incorrect type for site {type(value)}. Either a string or a ROOT.std.string is required.")
    
        self._efield_sim.string.assign(value)

    @property
    def refractivity_model(self):
        """Name of the atmospheric index of refraction model"""
        return str(self._refractivity_model)

    @refractivity_model.setter
    def refractivity_model(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for site {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._refractivity_model.string.assign(value)

    @property
    def refractivity_model_parameters(self):
        """Refractivity model parameters"""
        return self._refractivity_model_parameters

    @refractivity_model_parameters.setter
    def refractivity_model_parameters(self, value) -> None:
        # A list of strings was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._refractivity_model_parameters.clear()
            self._refractivity_model_parameters += value
        # A vector was given
        elif isinstance(value, ROOT.vector("double")):
            self._refractivity_model_parameters._vector = value
        else:
            raise ValueError(
                f"Incorrect type for refractivity_model_parameters {type(value)}. Either a list, an array or a ROOT.vector of unsigned shorts required."
            )


    @property
    def atmos_refractivity(self):
        """refractivity for each altitude at atmos_altiude table"""
        return self._atmos_refractivity


    @atmos_refractivity.setter
    def atmos_refractivity(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._atmos_refractivity.clear()
            self._atmos_refractivity += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._atmos_refractivity._vector = value
        else:
            raise ValueError(
                f"Incorrect type for atmos_refractivity {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )            
        

    @property
    def t_pre(self):
        """Starting time of antenna data collection time window. The window starts at t0+t_pre, thus t_pre is usually negative."""
        return self._t_pre[0]

    @t_pre.setter
    def t_pre(self, value):
        self._t_pre[0] = value

    @property
    def t_post(self):
        """Finishing time of antenna data collection time window. The window ends at t0+t_post."""
        return self._t_post[0]

    @t_post.setter
    def t_post(self, value):
        self._t_post[0] = value

    @property
    def t_bin_size(self):
        """Time bin size"""
        return self._t_bin_size[0]

    @t_bin_size.setter
    def t_bin_size(self, value):
        self._t_bin_size[0] = value




    @property
    def du_id(self):
        """Detector ID"""
        return self._du_id

    @du_id.setter
    def du_id(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._du_id.clear()
            self._du_id += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("int")):
            self._du_id._vector = value
        else:
            raise ValueError(
                f"Incorrect type for du_id {type(value)}. Either a list, an array or a ROOT.vector of float required."
            )

    @property
    def du_name(self):
        """Detector Name"""
        return self._du_name    


    @du_name.setter
    def du_name(self, value):
        # A list of strings was given
        if isinstance(value, list):
            # Clear the vector before setting
            self._du_name.clear()
            self._du_name += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("string")):
            self._du_name._vector = value
        else:
            raise ValueError(
                f"Incorrect type for du_name {type(value)}. Either a list or a ROOT.vector of strings required."
            )



    @property
    def t_0(self):
        """Time window t0"""
        return self._t_0

    @t_0.setter
    def t_0(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._t_0.clear()
            self._t_0 += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("float")):
            self._t_0._vector = value
        else:
            raise ValueError(
                f"Incorrect type for t_0 {type(value)}. Either a list, an array or a ROOT.vector of float required."
            )

    @property
    def p2p(self):
        """Peak 2 peak amplitudes (x,y,z,modulus)"""
        return self._p2p

    @p2p.setter
    def p2p(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._p2p.clear()
            self._p2p += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("float")):
            self._p2p._vector = value
        else:
            raise ValueError(
                f"Incorrect type for p2p {type(value)}. Either a list, an array or a ROOT.vector of float required."
            )        

    @property
    def trace_x(self):
        """Efield trace in X direction"""
        return self._trace_x

    @trace_x.setter
    def trace_x(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._trace_x.clear()
            self._trace_x += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._trace_x._vector = value
        else:
            raise ValueError(
                f"Incorrect type for trace_x {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )

    @property
    def trace_y(self):
        """Efield trace in Y direction"""
        return self._trace_y

    @trace_y.setter
    def trace_y(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._trace_y.clear()
            self._trace_y += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._trace_y._vector = value
        else:
            raise ValueError(
                f"Incorrect type for trace_y {type(value)}. Either a list, an array or a ROOT.vector of float required."
            )

    @property
    def trace_z(self):
        """Efield trace in Z direction"""
        return self._trace_z

    @trace_z.setter
    def trace_z(self, value):
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._trace_z.clear()
            self._trace_z += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._trace_z._vector = value
        else:
            raise ValueError(
                f"Incorrect type for trace_z {type(value)}. Either a list, an array or a ROOT.vector of float required."
            )
        
    @property
    def trace(self):
        return self._trace
    
    @trace.setter
    def trace(self, value):
        # value = np.concatenate((self._du_id, self._trace_x, self._trace_y, self._trace_z, self.t_bin_size))
        # A list was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._trace.clear()
            self._trace += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._trace._vector = value
        else:
            raise ValueError(
                f"Incorrect type for trace {type(value)}. Either a list, an array or a ROOT.vector of float required."
            )
                

    @property
    def du_x(self):
        """X position in site's referential"""
        return self._du_x

    @du_x.setter
    def du_x(self, value) -> None:
        # A list of strings was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._du_x.clear()
            self._du_x += value
        # A vector was given
        elif isinstance(value, ROOT.vector("float")):
            self._du_x._vector = value
        else:
            raise ValueError(
                f"Incorrect type for du_x {type(value)}. Either a list, an array or a ROOT.vector of floats required."
            )

    @property
    def du_y(self):
        """Y position in site's referential"""
        return self._du_y

    @du_y.setter
    def du_y(self, value) -> None:
        # A list of strings was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._du_y.clear()
            self._du_y += value
        # A vector was given
        elif isinstance(value, ROOT.vector("float")):
            self._du_y._vector = value
        else:
            raise ValueError(
                f"Incorrect type for du_y {type(value)}. Either a list, an array or a ROOT.vector of floats required."
            )

    @property
    def du_z(self):
        """Z position in site's referential"""
        return self._du_z

    @du_z.setter
    def du_z(self, value) -> None:
        # A list of strings was given
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._du_z.clear()
            self._du_z += value
        # A vector was given
        elif isinstance(value, ROOT.vector("float")):
            self._du_z._vector = value
        else:
            raise ValueError(
                f"Incorrect type for du_z {type(value)}. Either a list, an array or a ROOT.vector of floats required."
            )


#############################################################################################################################################################################################################################
#
#   RawMetaTree
#
#############################################################################################################################################################################################################################

@dataclass
## The class for storing meta-data for each event, that is meta to the shower and efield simulation (like, coreposition, array name, antenna selection, etc)
class RawMetaTree(MotherEventTree):
    """The class for storing data about the event generation that is meta to the shower/efield simulation itself"""


    _type: str = "rawmeta"

    _tree_name: str = "trawmeta"
    
    ### Array over wich the event was simulated (use "starshape" for...starshapes)
    _array_name: StdString = StdString("")
    
    #In the simulation, the coordinates are in "shower coordinates" whose origin is at the core position. So core position is always 0,0,0. The core position this represents in your array is meta the simulator 
    ### Core position with respect to the antenna array (undefined for neutrinos)
    _shower_core_pos: np.ndarray = field(default_factory=lambda: np.zeros(3, np.float32)) 
        
    #In the simulation, the origin of time is when the shower hits the ground.
    #There is a date (no time of the day, just the date) that can be used to get the magnetic field from the magnetic field model, but nothing else.
    #If the event you are simulating represents an event that happened at a specific time, you set here its second and nanosecond.     
    ### Unix second of the shower t0 (when the core traveling at c arrives at the ground?)
    _unix_second: np.ndarray = field(default_factory=lambda: np.zeros(1, np.uint32))    
    
    ### Unix nanosecond of the shower t0 (when the core traveling at c arrives at the ground?)
    _unix_nanosecond: np.ndarray = field(default_factory=lambda: np.zeros(1, np.uint32))    
    
    #The event you are simulating might be the result of several trials at input generation until you found one that has some chance of triggering. You need to store this info for efficiency studies.
    ### statistical weight given to the event
    _event_weight: np.ndarray = field(default_factory=lambda: np.zeros(1, np.uint32))    
    ### tested core positions
    _tested_cores: StdVectorList = field(default_factory=lambda:StdVectorList("vector<float>"))    

    @property
    def array_name(self):
        """array name"""
        return str(self._array_name)

    @array_name.setter
    def array_name(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for array_name {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._array_name.string.assign(value)

    @property
    def shower_core_pos(self):
        """Shower core position"""
        return np.array(self._shower_core_pos)

    @shower_core_pos.setter
    def shower_core_pos(self, value):
        self._shower_core_pos = np.array(value).astype(np.float32)
        self._tree.SetBranchAddress("shower_core_pos", self._shower_core_pos)         

    @property
    def unix_second(self):
        """Unix second of the shower t0 (when the core traveling at c arrives at the ground?)"""
        return self._unix_second[0]

    @unix_second.setter
    def unix_second(self, val: np.uint32) -> None:
        self._unix_second[0] = val


    @property
    def unix_nanosecond(self):
        """Unix nanosecond of the shower t0 (when the core traveling at c arrives at the ground?)"""
        return self._unix_nanosecond[0]

    @unix_nanosecond.setter
    def unix_nanosecond(self, val: np.uint32) -> None:
        self._unix_nanosecond[0] = val
        
        
    @property
    def event_weight(self):
        """The event statistical weight"""
        return self._event_weight[0]

    @event_weight.setter
    def event_weight(self, val: np.uint32) -> None:
        self._event_weight[0] = val
        
    @property
    def tested_cores(self):
        """tested cores"""
        return self._tested_cores
        
    @tested_cores.setter
    def long_pd_hadr(self, value):
        # A list was given  
        if (
            isinstance(value, list)
            or isinstance(value, np.ndarray)
            or isinstance(value, StdVectorList)
        ):
            # Clear the vector before setting
            self._tested_cores.clear()
            self._tested_cores += value
        # A vector of strings was given
        elif isinstance(value, ROOT.vector("vector<float>")):
            self._tested_cores._vector = value
        else:
            raise ValueError(
                f"Incorrect type for _tested_cores {type(value)}. Either a list, an array or a ROOT.vector of vector<float> required."
            )                

#############################################################################################################################################################################################################################
#
#   RawZHAiresTree
#
#############################################################################################################################################################################################################################
       
@dataclass
## The class for storing shower data for each event specific to ZHAireS only
class RawZHAireSTree(MotherEventTree):
    """The class for storing shower data for each event specific to ZHAireS only"""

    _type: str = "eventshowerzhaires"

    _tree_name: str = "teventshowerzhaires"

    # ToDo: we need explanations of these parameters
      

    _other_parameters: StdString = StdString("")



    @property
    def other_parameters(self):
        """Other parameters"""
        return str(self._other_parameters)

    @other_parameters.setter
    def other_parameters(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for other_parameters {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._other_parameters.string.assign(value)        

#############################################################################################################################################################################################################################
#
#   RawCoreasTree
#
#############################################################################################################################################################################################################################
     
@dataclass
## The class for storing shower data for each event specific to Coreas only
class RawCoreasTree(MotherEventTree):
    """The class for storing shower data for each event specific to Coreas only"""

    _type: str = "eventshowercoreas"

    _tree_name: str = "teventshowercoreas"

    _AutomaticTimeBoundaries: StdString = StdString("")
    _ResolutionReductionScale: StdString = StdString("")
    _GroundLevelRefractiveIndex: StdString = StdString("")
    _GPSSecs: StdString = StdString("")
    _GPSNanosecs: StdString = StdString("")
    _DepthOfShowerMaximum: StdString = StdString("")
    _DistanceOfShowerMaximum: StdString = StdString("")
    _GeomagneticAngle: StdString = StdString("")
    _nshow: StdString = StdString("")
    _ectmap: StdString = StdString("")
    _maxprt: StdString = StdString("")
    _radnkg: StdString = StdString("")
    _parallel: StdString = StdString("")
    _mumult: StdString = StdString("")
    _muaddi: StdString = StdString("")
    _parout: StdString = StdString("")
    _longi: StdString = StdString("")

    @property
    def AutomaticTimeBoundaries(self):
        """Automatic time boundaries"""
        return str(self._AutomaticTimeBoundaries)
    
    @AutomaticTimeBoundaries.setter
    def AutomaticTimeBoundaries(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for AutomaticTimeBoundaries {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._AutomaticTimeBoundaries.string.assign(value)
    
    @property
    def ResolutionReductionScale(self):
        """Resolution reduction scale"""
        return str(self._ResolutionReductionScale)
    
    @ResolutionReductionScale.setter
    def ResolutionReductionScale(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for ResolutionReductionScale {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._ResolutionReductionScale.string.assign(value)

    @property
    def GroundLevelRefractiveIndex(self):
        """Ground level refractive index"""
        return str(self._GroundLevelRefractiveIndex)
    
    @GroundLevelRefractiveIndex.setter
    def GroundLevelRefractiveIndex(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for GroundLevelRefractiveIndex {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._GroundLevelRefractiveIndex.string.assign(value)

    @property
    def GPSSecs(self):
        """GPS seconds"""
        return str(self._GPSSecs)
    
    @GPSSecs.setter
    def GPSSecs(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for GPSSecs {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._GPSSecs.string.assign(value)

    @property
    def GPSNanosecs(self):
        """GPS nanoseconds"""
        return str(self._GPSNanosecs)
    
    @GPSNanosecs.setter
    def GPSNanosecs(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for GPSNanosecs {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._GPSNanosecs.string.assign(value)

    @property
    def DepthOfShowerMaximum(self):
        """Depth of shower maximum"""
        return str(self._DepthOfShowerMaximum)
    
    @DepthOfShowerMaximum.setter
    def DepthOfShowerMaximum(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for DepthOfShowerMaximum {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._DepthOfShowerMaximum.string.assign(value)

    @property
    def DistanceOfShowerMaximum(self):
        """Distance of shower maximum"""
        return str(self._DistanceOfShowerMaximum)
    
    @DistanceOfShowerMaximum.setter
    def DistanceOfShowerMaximum(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for DistanceOfShowerMaximum {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._DistanceOfShowerMaximum.string.assign(value)

    @property
    def GeomagneticAngle(self):
        """Geomagnetic angle"""
        return str(self._GeomagneticAngle)
    
    @GeomagneticAngle.setter
    def GeomagneticAngle(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for GeomagneticAngle {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._GeomagneticAngle.string.assign(value)

    @property
    def nshow(self):
        """nshow"""
        return str(self._nshow)
    
    @nshow.setter
    def nshow(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for nshow {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._nshow.string.assign(value)

    @property
    def ectmap(self):
        """ectmap"""
        return str(self._ectmap)
    
    @ectmap.setter
    def ectmap(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for ectmap {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._ectmap.string.assign(value)

    @property
    def maxprt(self):
        """maxprt"""
        return str(self._maxprt)
    
    @maxprt.setter
    def maxprt(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for maxprt {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._maxprt.string.assign(value)

    @property
    def radnkg(self):
        """radnkg"""
        return str(self._radnkg)
    
    @radnkg.setter
    def radnkg(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for radnkg {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._radnkg.string.assign(value)

    @property
    def parallel(self):
        """parallel"""
        return str(self._parallel)
    
    @parallel.setter
    def parallel(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for parallel {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._parallel.string.assign(value)

    @property
    def mumult(self):
        """mumult"""
        return str(self._mumult)
    
    @mumult.setter
    def mumult(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for mumult {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._mumult.string.assign(value)

    @property
    def muaddi(self):
        """muaddi"""
        return str(self._muaddi)
    
    @muaddi.setter
    def muaddi(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for muaddi {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._muaddi.string.assign(value)

    @property
    def parout(self):
        """parout"""
        return str(self._parout)
    
    @parout.setter
    def parout(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for parout {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._parout.string.assign(value)

    @property
    def longi(self):
        """longi"""
        return str(self._longi)
    
    @longi.setter
    def longi(self, value):
        # Not a string was given
        if not (isinstance(value, str) or isinstance(value, ROOT.std.string)):
            raise ValueError(
                f"Incorrect type for longi {type(value)}. Either a string or a ROOT.std.string is required."
            )

        self._longi.string.assign(value)
