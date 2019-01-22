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
        self.cutTree.add(self.vertex,'vertex')

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
    def passMuon(self,cand):
        cuts = {
            'pt' : lambda cand: cand.pt()>5,
            'eta': lambda cand: abs(cand.eta())<2.4,
            'dxy': lambda cand: abs(cand.dxy())<0.5,
            'dz' : lambda cand: abs(cand.dz())<1,
            # muon best track type not in nanoaod?
            #'id' : lambda cand: (cand.isGlobal() or (cand.isTracker() and cand.nStations()>0)), # and cand.muonBestTrackType!=2,
            # global muon not in current version of nanoaod (but is in cmssw master)
            'id' : lambda cand: True, # isLooseMuon part of cut to include
            # note, this needs to be done after FSR, so probably move
            'iso': lambda cand: cand.pfRelIso03_all()<0.35,
            'sip': lambda cand: abs(cand.sip3d())<4,
        }
        # TODO ghost cleaning
        return all([cuts[cut](cand) for cut in cuts])

    def passTightMuon(self,cand):
        cuts = {
            'loose': lambda cand: self.passMuon(cand),
            'tight': lambda cand: cand.isPFcand() or (cand.highPtId()>0 and cand.pt()>200),
        }
        return all([cuts[cut](cand) for cut in cuts])

    def passElectron(self,cand):
        pt = cand.pt()
        eta = abs(cand.eta() + cand.deltaEtaSC())
        # not included in nanoaod yet (but is on master)
        #mva = cand.mvaFall17V1Iso() # note, should be V2
        cuts = {
            'pt' : lambda cand: cand.pt()>7,
            'eta': lambda cand: abs(cand.eta())<2.5,
            'dxy': lambda cand: abs(cand.dxy())<0.5,
            'dz' : lambda cand: abs(cand.dz())<1,
            'id' : lambda cand: cand.mvaFall17Iso_WPL()>0.5, # tmp
            #'id' : lambda cand: (pt<=10 and ((eta<0.8                and mva>0.5739521065342641)\
            #                            or   (eta>=0.8 and eta<1.479 and mva>0.5504628790992929)\
            #                            or   (eta>=1.479             and mva>0.5924627534389098)))\
            #                    or\
            #                    (pt>10  and ((eta<0.8                and mva>-0.03391387993354392)\
            #                            or   (eta>=0.8 and eta<1.479 and mva>-0.018451958064666783)\
            #                            or   (eta>=1.479             and mva>-0.38565459150737535))),
            'iso': lambda cand: cand.pfRelIso03_all()<0.35,
            'sip': lambda cand: abs(cand.sip3d())<4,
        }
        return all([cuts[cut](cand) for cut in cuts])

    def passPhoton(self,cand):
        # note, not possible given the pt cuts in nanoaod (5 vs 2 needed)
        return True

    def passJet(self,cand):
        return True

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
        muons = self.getCands(self.muons,self.passMuon)
        electrons = self.getCands(self.electrons,self.passElectron)

        # cross clean electrons
        tightmuons = self.getCands(self.muons,self.passTightMuon)
        electrons = self.cleanCands(electrons,tightmuons,0.05)

        leps = electrons + muons
        if len(leps)<4: return candidate # need at least 4 leptons

        # TODO FSR not added

        # get the candidates
        zzCands = []
        for quad in itertools.permutations(leps,4):
            # require +-+-
            if quad[0].charge()+quad[1].charge()!=0: continue
            if quad[2].charge()+quad[3].charge()!=0: continue
            # require same type
            if quad[0].collName!=quad[1].collName: continue
            if quad[2].collName!=quad[3].collName: continue
            # require hzz cuts
            z1 = DiCandidate(quad[0],quad[1])
            m1 = z1.M()
            z2 = DiCandidate(quad[2],quad[3])
            m2 = z2.M()
            if m1<40 or m1>120 or m2<12 or m2>120: continue
            # combinatoric cuts
            keep = True
            for i,j in itertools.combinations(range(4),2):
                dicand = DiCandidate(quad[i],quad[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if quad[i].charge()+quad[j].charge()==0 and m<4.: keep = False
            if not keep: continue
            # pt cuts
            pts = sorted([quad[i].pt() for i in range(4)],reverse=True)
            if pts[0]<20: continue
            if pts[1]<10: continue
            # smart cut
            if quad[0].collName==quad[2].collName:
                if quad[0].charge()!=quad[2].charge():
                   pairs = [[0,2],[1,3]]
                else:
                   pairs = [[0,3],[1,2]]
                zX = DiCandidate(*[quad[i] for i in pairs[0]])
                zY = DiCandidate(*[quad[i] for i in pairs[1]])
                za = zX if abs(zX.M()-ZMASS) < abs(zY.M()-ZMASS) else zY
                zb = zY if abs(zX.M()-ZMASS) < abs(zY.M()-ZMASS) else zX
                if abs(za.M()-ZMASS)<abs(z1.M()-ZMASS) and zb.M()<12: continue
            # 4 body cut
            zz = CompositeCandidate(*quad)
            if zz.M()<70: continue
            # its a good candidate
            zzCands += [quad]
        if not zzCands: return candidate

        # sort by best z then st
        # note: should be selected by P_sig/P_bkg (but dont know what this is)
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

    def vertex(self,cands):
        # !isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2
        return self.event.PV_npvsGood()

    def trigger(self,cands):
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

    #parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('hzz'), help='Input files')
    #parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('data'), help='Input files')
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
