#!/usr/bin/env python
import argparse
import logging
import sys

from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR, getTestFiles

from Candidates import *

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("MonoHZZAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class MonoHZZAnalysis(AnalysisBase):
    '''
    MonoHZZ analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','monoHZZTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MonoHZZTree')
        self.preselection = 'muons_count+electrons_count>3'
        super(MonoHZZAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.fourLoose,'fourLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',5])

        # trigger

        # 4 lepton
        self.addComposite('4l')
        self.addCompositeMet('4lmet')

        # z1 leptons
        self.addDiLepton('z1')
        self.addLepton('z11')
        self.addLepton('z12')

        # z2 leptons
        self.addDiLepton('z2')
        self.addLepton('z21')
        self.addLepton('z22')

        # met
        self.addMet('met')


    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z11' : None,
            'z12' : None,
            'z21' : None,
            'z22' : None,
            'z1' : None,
            'z2' : None,
            '4l' : None,
            '4lmet' : None,
            'met': self.pfmet,
        }

        # get leptons
        leps = self.electrons+self.muons
        if len(leps)<4: return candidate # need at least 4 leptons


        # get the candidates
        zzCands = []
        for quad in itertools.permutations(leps,4):
            # require +-+-
            if quad[0].charge()+quad[1].charge()!=0: continue
            if quad[2].charge()+quad[3].charge()!=0: continue
            # require same type
            if quad[0].collName!=quad[1].collName: continue
            if quad[2].collName!=quad[3].collName: continue
            # require deltaR seperation of 0.02 m(ll)>12
            keep = True
            for i,j in itertools.combinations(range(4),2):
                dicand = DiCandidate(quad[i],quad[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
            if not keep: continue
            # require lead e/m pt > 25
            ems = [cand for cand in quad if cand.collName in ['electrons','muons']]
            # its a good candidate
            zzCands += [quad]
        if not zzCands: return candidate

        # sort by best z then st
        bestZ = 0
        highestSt = 0
        bestCand = []
        for quad in zzCands:
            z = DiCandidate(quad[0],quad[1])
            zmass = z.M()
            st = sum([cand.pt() for cand in quad[2:]])
            if abs(ZMASS-zmass)<(ZMASS-bestZ):
                bestZ = zmass
                highestSt = st
                bestCand = quad
            elif abs(ZMASS-zmass)==(ZMASS-bestZ):
                bestZ = zmass
                highestSt = st
                bestCand = quad

        z11 = bestCand[0] if bestCand[0].pt()>bestCand[1].pt() else bestCand[1]
        z12 = bestCand[1] if bestCand[0].pt()>bestCand[1].pt() else bestCand[0]
        z21 = bestCand[2] if bestCand[2].pt()>bestCand[3].pt() else bestCand[3]
        z22 = bestCand[3] if bestCand[2].pt()>bestCand[3].pt() else bestCand[2]

        candidate['z11'] = z11
        candidate['z12'] = z12
        candidate['z21'] = z21
        candidate['z22'] = z22
        candidate['z1'] = DiCandidate(z11,z12)
        candidate['z2'] = DiCandidate(z21,z22)
        candidate['4l'] = CompositeCandidate(z11,z12,z21,z22)
        candidate['4lmet'] = MetCompositeCandidate(self.pfmet,z11,z12,z21,z22)

        return candidate


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z11','z12','z21','z22']:
            chanString += self.getCollectionString(cands[c])
        return chanString


    ###########################
    ### analysis selections ###
    ###########################
    def fourLoose(self,cands):
        return len(self.electrons+self.muons)>=4

    def trigger(self,cands):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
        triggerNames = {
            'DoubleMuon'     : [
                'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL',
                'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL',
                'DiMu9_Ele9_CaloIdL_TrackIdL',
                'TripleMu_12_10_5',
            ],
            'DoubleEG'       : [
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'Ele16_Ele12_Ele8_CaloIdL_TrackIdL',
                'Mu8_DiEle12_CaloIdL_TrackIdL',
            ],
            'MuonEG'         : [
                'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
            ],
            'SingleMuon'     : [
                'IsoMu24',
                'IsoTkMu24',
                'Mu45_eta2p1',
                'Mu50',
            ],
            'SingleElectron' : [
                'Ele27_WPTight_Gsf',
                'Ele27_eta2p1_WPLoose_Gsf',
                'Ele45_WPLoose_Gsf',
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
                var = '{0}Pass'.format(trigger)
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

    #parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('hzz'), help='Input files')
    parser.add_argument('--inputFiles', type=str, nargs='*', default='hzz.root', help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='monoHZZTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    zzAnalysis = MonoHZZAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MonoHZZTree',
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
