# common utilities for analyzers
import os
import sys
import glob

import ROOT

from DevTools.Utilities.utilities import ZMASS, getCMSSWVersion
from DevTools.Utilities.hdfsUtils import get_hdfs_root_files

def deltaPhi(phi0,phi1):
    result = phi0-phi1
    while result>ROOT.TMath.Pi():
        result -= 2*ROOT.TMath.Pi()
    while result<=-ROOT.TMath.Pi():
        result += 2*ROOT.TMath.Pi()
    return result

def deltaR(eta0,phi0,eta1,phi1):
    deta = eta0-eta1
    dphi = deltaPhi(phi0,phi1)
    return ROOT.TMath.Sqrt(deta**2+dphi**2)

def getTestFiles(sample):
    sampleMap = {
        'hzz': '/store/mc/RunIIFall17NanoAOD/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/22014A1F-BD44-E811-A374-0025905B85B8.root',
    }

    if sample in sampleMap:
        return sampleMap[sample]

    return []
