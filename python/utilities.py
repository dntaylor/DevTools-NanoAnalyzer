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
        'hzz': '/store/mc/RunIIFall17NanoAODv4/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/40000/4AAF9451-2C48-8247-867A-3BB8980C4893.root',
        'data': '/store/data/Run2017B/DoubleMuon/NANOAOD/Nano14Dec2018-v1/280000/FB901F01-98AA-214F-A2C2-D67630861952.root',
    }

    if sample in sampleMap:
        return sampleMap[sample]

    return []
