import numpy as np
import os
import sys
import argparse
from ProjectUtils.config_editor import *
from ProjectUtils.utils_optim import *
from ProjectUtils.runSimulator import *

def GenerateRandomDesign(designParams, data):
    GeneratedDesign = dict()
    for k, v in designParams.items():
        tmp = np.random.uniform(v["min"], v["max"])
        GeneratedDesign[k] = tmp
    return GeneratedDesign

parser = argparse.ArgumentParser("Run Random Design Point")
parser.add_argument("-c", "--config", help="Path to the config json file", required=True)
parser.add_argument("-r", "--HowToRun", help="insideSig = 0, Visualize = 1, outsideSig = 2", type=int, default=1)

args = parser.parse_args()
if not os.path.exists(args.config):
    print("Config file does not exist")
    sys.exit(1)
    
data = ReadJsonFile(args.config)
designParams = GetDesignParamNames(data["Signatures_NParams"], data["DesignParams_Ranges"])
Valid = False
while(not Valid):
    GeneratedDesign = GenerateRandomDesign(designParams, data)
    Valid = CheckConstraints(GeneratedDesign, data["Constraints"], data["Signatures_NParams"])
RunSim(GeneratedDesign, data, args.HowToRun)
