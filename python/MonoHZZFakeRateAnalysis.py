#!/usr/bin/env python
import argparse
import logging
import sys

from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR, getTestFiles

from Candidates import *
from leptonIds import *

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("MonoHZZFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class MonoHZZFakeRateAnalysis(AnalysisBase):
    '''
    MonoHZZFakeRate analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','monoHZZFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MonoHZZFakeRateTree')
        self.preselection = 'muons_count+electrons_count>2'
        super(MonoHZZFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.threeLoose,'threeLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')
        self.cutTree.add(self.vertex,'vertex')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',4])

        # trigger

        # z1 leptons
        self.addDiLepton('z')
        self.addLepton('z1')
        self.tree.add(lambda cands: passHZZLooseNoIso(cands['z1']), 'z1_passLooseNoIso','I')
        self.tree.add(lambda cands: passHZZLoose(     cands['z1']), 'z1_passLoose','I')
        self.tree.add(lambda cands: passHZZTight(     cands['z1']), 'z1_passTight','I')
        self.addLepton('z2')
        self.tree.add(lambda cands: passHZZLooseNoIso(cands['z2']), 'z2_passLooseNoIso','I')
        self.tree.add(lambda cands: passHZZLoose(     cands['z2']), 'z2_passLoose','I')
        self.tree.add(lambda cands: passHZZTight(     cands['z2']), 'z2_passTight','I')

        # z2 leptons
        self.addLepton('l')
        self.tree.add(lambda cands: passHZZLooseNoIso(cands['l']), 'l_passLooseNoIso','I')
        self.tree.add(lambda cands: passHZZLoose(     cands['l']), 'l_passLoose','I')
        self.tree.add(lambda cands: passHZZTight(     cands['l']), 'l_passTight','I')

        # met
        self.addMet('met')


    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            'l'  : None,
            'z1' : None,
            'met': self.pfmet,
        }


        # get leptons
        muons = self.getCands(self.muons,passHZZLooseMuonNoIso)
        electrons = self.getCands(self.electrons,passHZZLooseElectronNoIso)

        tightmuons = self.getCands(self.muons,passHZZTightMuon)
        tightelectrons = self.getCands(self.electrons,passHZZTightElectron)

        # cross clean electrons
        electrons = self.cleanCands(electrons,tightmuons,0.05)

        leps = electrons + muons
        if len(leps)!=3: return candidate # need at least 4 leptons

        # TODO FSR not added

        # get the candidates
        zlCands = []
        for tri in itertools.permutations(leps,3):
            # require +-+-
            if tri[0].charge()+tri[1].charge()!=0: continue
            # require same type
            if tri[0].collName!=tri[1].collName: continue
            # require hzz cuts
            z1 = DiCandidate(tri[0],tri[1])
            m1 = z1.M()
            if m1<60 or m1>120: continue
            # combinatoric cuts
            keep = True
            for i,j in itertools.combinations(range(3),2):
                dicand = DiCandidate(tri[i],tri[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if tri[i].charge()+tri[j].charge()==0 and m<4.: keep = False
            if not keep: continue
            # pt cuts
            pts = sorted([tri[i].pt() for i in range(3)],reverse=True)
            if pts[0]<20: continue
            if pts[1]<10: continue
            # its a good candidate
            zlCands += [tri]
        if not zlCands: return candidate

        # sort by best z
        bestZ = 0
        bestCand = []
        for tri in zlCands:
            z = DiCandidate(tri[0],tri[1])
            zmass = z.M()
            if abs(ZMASS-zmass)<(ZMASS-bestZ):
                bestZ = zmass
                bestCand = tri

        z1 = bestCand[0] if bestCand[0].pt()>bestCand[1].pt() else bestCand[1]
        z2 = bestCand[1] if bestCand[0].pt()>bestCand[1].pt() else bestCand[0]
        l = bestCand[2]

        candidate['z1'] = z1
        candidate['z2'] = z2
        candidate['l'] = l
        candidate['z'] = DiCandidate(z1,z2)

        return candidate


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2','l']:
            chanString += self.getCollectionString(cands[c])
        return chanString


    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self):
        ne = len(self.electrons)
        nm = len(self.muons)
        return ne+nm>=3

    def vertex(self):
        # !isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2
        return self.event.PV_npvsGood()

    def trigger(self):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
        # 2018 commented
        triggerNames = {
            'DoubleMuon'     : [
                'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                #'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8',
                #'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8',
                'TripleMu_10_5_5_DZ',
                'TripleMu_12_10_5',
            ],
            'DoubleEG'       : [
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'DoubleEle33_CaloIdL_MW',
                'Ele16_Ele12_Ele8_CaloIdL_TrackIdL',
            ],
            'MuonEG'         : [
                #'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
                'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                'Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'DiMu9_Ele9_CaloIdL_TrackIdL_DZ',
                'Mu8_DiEle12_CaloIdL_TrackIdL',
                'Mu8_DiEle12_CaloIdL_TrackIdL_DZ',
            ],
            'SingleMuon'     : [
                'IsoMu27',
            ],
            'SingleElectron' : [
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf',
            ],
        }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleMuon', 
            'DoubleEG', 
            'MuonEG',
            'SingleMuon',
            'SingleElectron',
        ]
        # reject triggers if they are in another dataset
        # looks for the dataset name in the filename
        # for MC it accepts all
        reject = True if self.event.isData()>0.5 else False
        for dataset in datasets:
            # if we match to the dataset, start accepting triggers
            if dataset in self.fileNames[0]: reject = False
            for trigger in triggerNames[dataset]:
                var = 'HLT_{0}'.format(trigger)
                passTrigger = getattr(self.event,var)()
                if passTrigger>0.5:
                    # it passed the trigger
                    # in data: reject if it corresponds to a higher dataset
                    return False if reject else True
            # dont check the rest of data
            if dataset in self.fileNames[0]: break
        return False


    




def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('data'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='monoHZZFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    zzAnalysis = MonoHZZFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MonoHZZFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        shift = args.shift,
    )

    try:
       zzAnalysis.analyze()
       zzAnalysis.finish()
    except KeyboardInterrupt:
       zzAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
