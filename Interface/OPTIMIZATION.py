RCP = "RCP0" # no climate change
filename = "DATA/DATA.h5" # Test data from Central Finland with 3579 forest stands
scenario = "NFS"
extension = "test" # some additional info to the saved output 

import wget
import os
import pandas as pd
import sys
import MultiFunctionalOptimization as MFO

from importlib import reload
reload(MFO)
module_path = os.path.abspath(os.path.join(''))
#module_path+"\\Documents"
if module_path not in sys.path:
    sys.path.append(module_path+"/py_class")

class OptGUI:
    path_to_zip_file = module_path + "/DATA.zip"
    directory_to_extract_to = module_path+"/DATA"

    try:
        import zipfile
        with zipfile.ZipFile(path_to_zip_file, 'r') as zip_ref:
            zip_ref.extractall(directory_to_extract_to)
        print("Zip Extracted")
    except:
        print("Zip already existed")
    
    def __init__(self,scenario):
        self.scenario = scenario
        
        self.mfo = MFO.MultiFunctionalOptimization() 
        self.mfo.readData(module_path+"\\"+filename,
                     sampleRatio=1 #If no sample ratio given, the ratio is assumed to be 1.
                    )
        self.columnTypes = {
            'i_Vm3':(float,"Relative to Area"),
            'Harvested_V':(float,"Relative to Area"),
            'Harvested_V_log_under_bark':(float,"Relative to Area"), 
            'Harvested_V_pulp_under_bark':(float,"Relative to Area"),
            'Harvested_V_under_bark':(float,"Relative to Area"), 
            'Biomass':(float,"Relative to Area"),
            'ALL_MARKETED_MUSHROOMS':(float,"Relative to Area"), 
            'BILBERRY':(float,"Relative to Area"), 
            'COWBERRY':(float,"Relative to Area"),
            'HSI_MOOSE':(float,"Relative to Area"),
            'CAPERCAILLIE':(float,"Relative to Area"), 
            'HAZEL_GROUSE':(float,"Relative to Area"), 
            'V_total_deadwood':(float,"Relative to Area"), 
            'N_where_D_gt_40':(float,"Relative to Area"),
            'prc_V_deciduous':(float,"Relative to Area"),
            'CARBON_SINK':(float,"Relative to Area"), 
            'Recreation':(float,"Relative to Area"),
            'Scenic':(float,"Relative to Area")
        }

        self.mfo.calculateTotalValuesFromRelativeValues(columnTypes=self.columnTypes)

        self.regimeClassNames = {"regimeClass0name":"CCF","regimeClass1name":"SA","regimeClass2name":"Broadleave"}
        self.regimeClassregimes = {"regimeClass0regimes":["CCF_3","CCF_4","BAUwGTR"],"regimeClass1regimes":["SA"],"regimeClass2regimes":["BAUwT_B", "BAUwT_5_B", "BAUwT_15_B", "BAUwT_30_B", "BAUwT_GTR_B"]}

        self.mfo.addRegimeClassifications(regimeClassNames = self.regimeClassNames,regimeClassregimes=self.regimeClassregimes)

        self.mfo.finalizeData(initialRegime="initial_state")

        if self.scenario == 'NFS':
            
            self.wood_production_bioenergy = { 
            # Increment - target 2025
            "Total_i_Vm3_2025": ["Total annual timber volume increment by 2025 (m3)",
                                 "Total_i_Vm3",
                                 "max","targetYearWithSlope","sum",2025], 
            # Increment - target 2050
            "Total_i_Vm3_2050": ["Total annual timber volume increment by 2050 (m3)",
                                 "Total_i_Vm3",
                                 "max","targetYearWithSlope","sum",2050], 
            # Harvested roundwood - target 2025
            "Total_Harvested_V_2025" :["Total annual harvested timber volume by 2025 (log & pulp) (m3)",
                                       "Total_Harvested_V",
                                       "max","targetYearWithSlope","sum",2025], 
            # Harvested biomass - target 2025
            "Total_Biomass_2025": ["Total annual harvested biomass volume by 2025 (m3)",
                                   "Total_Biomass",
                                   "max","targetYearWithSlope","sum",2025]
            }
            
            self.nonwood = { 
            # Bilberry - no decline, maximise it
            "Relative_BILBERRY": ["Bilberry yield (relative to 2016, max minimum over yrs)",
                                  "Relative_Total_BILBERRY",
                                  "max","min","sum"],
            # Cowberry - no decline, maximise it
            "Relative_COWBERRY": ["Cowberry yield (relative to 2016, max minimum over yrs)",
                                  "Relative_Total_COWBERRY",
                                  "max","min","sum"],
            # Mushrooms - no decline, maximise it
            "Relative_ALL_MARKETED_MUSHROOMS": ["All marketed mushroom yield (relative to 2016, max minimum over yrs)",
                                 "Relative_Total_ALL_MARKETED_MUSHROOMS",
                                 "max","min","sum"]    
            }
            
            self.game = {
            # HSI moose - maximise  
            "Sum_Total_HSI_MOOSE": ["Total habitat index for MOOSE (max average over all years)",
                                   "Total_HSI_MOOSE",
                                   "max","average","sum"],
            # HSI hazel grouse - maximise
            "Sum_Total_HAZEL_GROUSE": ["Total habitat index for HAZEL_GROUSE (max average over yrs)",
                                   "Total_HAZEL_GROUSE",
                                   "max","average","sum"],
            # HSI capercaillie - maximise
            "Sum_Total_CAPERCAILLIE": ["Total habitat index for CAPERCAILLIE (max average over yrs)",
                                   "Total_CAPERCAILLIE",
                                   "max","average","sum"]    
            }
            
            self.biodiversity = {
            # Deadwood - target 2025
            "Average_Deadwood_V_2025": ["Average Deadwood volume by 2025 (m3/ha)", 
                                        "V_total_deadwood",
                                        "max", "targetYear", "areaWeightedAverage", 2025], 
            # Large trees (>40 cm) - maximise
            "Total_N_where_D_gt_40": ["Total No. of trees diameter >= 40 cm  (max end value)",
                                      "Total_N_where_D_gt_40",
                                      "max","lastYear","sum"],    
            # Deciduous tree volume - maximise
            "Total_prc_V_deciduous":  ["Total %-share of deciduous trees (related to V) (max end value)", 
                                       "Total_prc_V_deciduous",
                                       "max", "lastYear","sum"],
            # Conservation regime - target
            "Ratio_CCF_forests": ["Ratio of BC sites in commercial forests (%, CCF_3, CCF_4 and BAUwGTR)",
                                  "CCF_forests",
                                  "max","firstYear","areaWeightedAverage"]
            }
            
            self.climate_regulation = {
            # Carbon sink - target 2025
            "Total_CARBON_SINK_2025": ["Total sequestration in carbon dioxide by 2025 (t CO2)",
                                       "Total_CARBON_SINK",
                                       "max","targetYearWithSlope","sum",2025] 
            }
            
            self.recreation = {
            # Recreation index - maximise
            "Sum_Total_Recreation" : ["Total Recreation index (max minimum over yrs)",
                                      "Total_Recreation",
                                      "max","min","sum"],
            
            # Scenic index - maximise
            "Sum_Total_Scenic" : ["Total Scenic index (max minimum over yrs)",
                                  "Total_Scenic",
                                  "max","min","sum"]
            }
            
            self.resilience = {
            # CC adaption regimes - maximise
            "Ratio_Broadleave_forests": ["Ratio of adaptive management regimes (%, increasing broadleave share)",
                                         "Broadleave_forests",
                                         "max","firstYear","areaWeightedAverage"]
            }
            
            self.objectives = {
                      **self.wood_production_bioenergy,
                      **self.nonwood,
                      **self.game,
                      **self.biodiversity,
                      **self.climate_regulation,
                      **self.recreation,
                      **self.resilience
            }
            
            print("objectives for NFS loaded")
            
        if self.scenario == 'BDS':
            
            self.wood_production_bioenergy = { 
            # Harvested roundwood - maximise (even flow)
            "Average_Harvested_V" : ["Average harvested timber volume (log & pulp) (m3/ha, evenflow)",
                                     "Harvested_V",
                                     "max","min","areaWeightedAverage"]
            }
            
            self.game = {
            # HSI moose - maximise       
            "Sum_Total_HSI_MOOSE": ["Total habitat index for MOOSE (max average over all years)",
                                   "Total_HSI_MOOSE",
                                   "max","average","sum"],
            # HSI hazel grouse - maximise
            "Sum_Total_HAZEL_GROUSE": ["Total habitat index for HAZEL_GROUSE (max average over yrs)",
                                   "Total_HAZEL_GROUSE",
                                   "max","average","sum"],
            # HSI carpercaillie - maximise
            "Sum_Total_CAPERCAILLIE": ["Total habitat index for CAPERCAILLIE (max average over yrs)",
                                   "Total_CAPERCAILLIE",
                                   "max","average","sum"]    
            }
            
            self.biodiversity = {
            # Deadwood - target 2050, increase by XX%
            "relative_Amount_Deadwood_2050" : ["Total Deadwood volume by 2050 (%, relative to 2016 values)",
                                               "Relative_Total_V_total_deadwood",
                                               "max","targetYearWithSlope","sum",2050],
            # Large trees - target 2050, increase by XX%
            #"relative_N_where_D_gt_40_2050": ["Total No. of trees diameter >= 40 cm  by 2050 (%, relative to 2016 values)",
            #                                  "Relative_Total_N_where_D_gt_40",
            #                                  "max","targetYear","sum",2050],
            # Deciduous tree volume - target 2050, increase by XX% 
            "relative_prc_V_deciduous_2050": ["Total share of deciduous trees by 2050 (related to V) (%, relative to 2016 values)",
                                              "Relative_Total_prc_V_deciduous",
                                              "max","targetYearWithSlope","sum",2050],
            # Regime SA - target
            "Ratio_CCF_forests": ["Ratio of BC sites in commercial forests (%, CCF_3, CCF_4 and BAUwGTR)",
                                  "CCF_forests",
                                  "max","firstYear","areaWeightedAverage"],
            # Conservation regimes - target
            "Ratio_SA_forests": ["Ratio of protected areas (%, SA forests)",
                                 "SA_forests",
                                 "max","firstYear","areaWeightedAverage"]    
            
            }
            
            self.recreation = {
            # Recreation index - maximise
            "Sum_Total_Recreation" : ["Total Recreation index (max minimum over yrs)",
                                      "Total_Recreation",
                                      "max","min","sum"],
            # Scenic index - maximise
            "Sum_Total_Scenic" : ["Total Scenic index (max minimum over yrs)",
                                  "Total_Scenic",
                                  "max","min","sum"]
            }
            
            self.objectives = {
                      **self.wood_production_bioenergy,
                      **self.game,
                      **self.biodiversity,
                      **self.recreation
            }
            
            print("objectives for BDS loaded")
            
        if self.scenario == 'BES':
            
            self.wood_production_bioenergy = { 
            # Harvested roundwood - maximise even flow
            "Average_Harvested_V" : ["Average harvested timber volume (log & pulp) (m3/ha, evenflow)",
                                     "Harvested_V",
                                     "max","min","areaWeightedAverage"],
            # Harvested biomass - maximise even flow
            "Biomass_Evenflow": ["Average harvested biomass volume (m3/ha, evenflow)",
                                 "Biomass",
                                 "max","min","areaWeightedAverage"]
            }
            
            self.biodiversity = {
            # Deadwood - no decline (no target value)
            "relative_Amount_Deadwood_2050" : ["Total Deadwood volume by 2050 (%, relative to 2016 values)",
                                               "Relative_Total_V_total_deadwood",
                                               "max","targetYearWithSlope","sum",2050],
            # Large trees - no decline (no target value)
            #"relative_N_where_D_gt_40_2050": ["Total No. of trees diameter >= 40 cm  by 2050 (%, relative to 2016 values)",
            #                                  "Relative_Total_N_where_D_gt_40",
            #                                  "max","targetYear","sum",2050],
            # Deciduous tree volume - no decline (no target value)
            "relative_prc_V_deciduous_2050": ["Total share of deciduous trees by 2050 (related to V) (%, relative to 2016 values)",
                                              "Relative_Total_prc_V_deciduous",
                                              "max","targetYearWithSlope","sum",2050]
            }
            
            self.recreation = {
            # Recreation index - maximise
            "Sum_Total_Recreation" : ["Total Recreation index (max minimum over all years)",
                                      "Total_Recreation",
                                      "max","min","sum"],
            # Scenic index - maximise
            "Sum_Total_Scenic" : ["Total Scenic index (max minimum over all years)",
                                  "Total_Scenic",
                                  "max","min","sum"]
            }
            
            self.objectives = {
                      **self.wood_production_bioenergy,
                      **self.biodiversity,
                      **self.recreation,
            }
            
            print("objectives for BES loaded")
            
        self.initialValues = {"Total_i_Vm3":107*10**6 / 19,               # from National Forest Policy            
                         "Total_Harvested_V": 72.3*10**6 / 19,       # from National Forest Policy 
                         "Total_Biomass": 2.9*10**6 / 19,            # from National Forest Policy  
                         "Total_CARBON_SINK" : 34.1*10**6 / 19,      # from National Forest Policy  
                                                    
                         "SA_forests" : 0.106,     # from ForestStatistics 2018
                         "CCF_forests" : 0.015,    # from ForestStatistics 2018
                         "BAUwGTR_forests":0.015}  # from ForestStatistics 2018    

        self.mfo.defineObjectives(self.objectives,initialValues = self.initialValues)

        self.CCFregimes = [regime for regime in self.mfo.regimes if "CCF" in regime] + ["SA"]

        self.constraintTypes = {"CCFonPeat":["Allowed regimes","Only CCF on peat lands",self.CCFregimes,"PEAT"],"SPATIAL":["SPATIAL"]}

        self.mfo.defineConstraints(self.constraintTypes)

        self.mfo.calculateObjectiveRanges(debug=False)

        self.mfo.showGUI(debug=False)
