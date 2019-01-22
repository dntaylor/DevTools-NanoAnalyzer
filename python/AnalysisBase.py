# AnalysisBase.py

import logging
import os
import sys
import math
import time

sys.argv.append('-b')
import ROOT
sys.argv.pop()
from array import array

from CutTree import CutTree
from AnalysisTree import AnalysisTree

from utilities import deltaR, deltaPhi
from DevTools.Utilities.utilities import getCMSSWVersion
from Candidates import *

class AnalysisBase(object):
    '''
    Analysis Tree
    '''

    def __init__(self,**kwargs):
        inputFileNames = kwargs.pop('inputFileNames',[])
        inputTreeDirectory = kwargs.pop('inputTreeDirectory','')
        inputTreeName = kwargs.pop('inputTreeName','Events')
        outputFileName = kwargs.pop('outputFileName','analysisTree.root')
        outputTreeName = kwargs.pop('outputTreeName','AnalysisTree')
        self.shift = kwargs.pop('shift','')
        self.events = []
        self.outputTreeName = outputTreeName
        # preselection
        if not hasattr(self,'preselection'): self.preselection = '1'
        # input files
        self.fileNames = []
        if os.path.isfile('PSet.py'):                # grab input files from crab pset
            import PSet
            self.fileNames = list(PSet.process.source.fileNames)
        elif isinstance(inputFileNames, basestring): # inputFiles is a file name
            if inputFileNames.startswith('/store'):  # xrootd access
                self.fileNames += [inputFileNames]
            if os.path.isfile(inputFileNames):       # single file
                if inputFileNames[-4:] == 'root':    # file is a root file
                    self.fileNames += [inputFileNames]
                else:                                # file is list of files
                    with open(inputFileNames,'r') as f:
                        for line in f:
                            self.fileNames += [line.strip()]
        else:
            self.fileNames = inputFileNames          # already a python list or a cms.untracked.vstring()
        if not isinstance(outputFileName, basestring): # its a cms.string(), get value
            outputFileName = outputFileName.value()
        if not self.fileNames:
            logging.error('No files found')
            raise ValueError
        # test for hdfs
        #self.hasHDFS = os.path.exists('/hdfs/store/user')
        self.hasHDFS = False
        # input tchain
        self.treename = inputTreeName
        self.totalEntries = 0
        self.numLumis = 0
        self.numEvents = 0
        self.summedWeights = 0
        self.summedWeightsLHEScale = [0]*9
        logging.info('Getting Lumi information')
        #self.skims = {}
        #self.tfiles = []
        for f,fName in enumerate(self.fileNames):
            if fName.startswith('/store/user'):
                fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
            elif fName.startswith('/store'):
                fName = '{0}/{1}'.format('root://cmsxrootd.fnal.gov/',fName)
            tfile = ROOT.TFile.Open(fName)
            tree = tfile.Get(self.treename)
            self.totalEntries += tree.GetEntries()
            if not hasattr(self,'version'):
                tree.GetEntry(1)
                if hasattr(tree,'provenance'):
                    ver = tree.provenance[0].split('_')
                    self.version = ''.join([ver[1],ver[2],'X'])
                else:
                    self.version = getCMSSWVersion()
            lumitree = tfile.Get('LuminosityBlocks')
            for entry in lumitree:
                self.numLumis += 1
            runtree = tfile.Get('Runs')
            for entry in runtree:
                if hasattr(runtree,'genEventCount'):
                    self.numEvents += runtree.genEventCount
                    self.summedWeights += runtree.genEventSumw
                    for i in range(9):
                        self.summedWeightsLHEScale[i] += runtree.LHEScaleSumw[i]
            #self.tfiles += [tfile]
            tfile.Close('R')
        logging.info('Analysis is running with version {0}'.format(self.version))
        logging.info("Will process {0} lumi sections with {1} events ({2}).".format(self.numLumis,self.numEvents,self.summedWeights))
        self.flush()
        if not len(self.fileNames): raise Exception
        # tfile
        self.outfile = ROOT.TFile(outputFileName,"recreate")
        # cut tree
        self.cutTree = CutTree()
        # analysis tree
        self.tree = AnalysisTree(outputTreeName)
        self.eventsStored = 0

        # some things we always need:

        # pileup
        self.tree.add(lambda cands: self.event.nOtherPV()+1, 'numVertices', 'I')
        self.tree.add(lambda cands: self.event.fixedGridRhoFastjetAll(), 'rho', 'F')

        # gen
        self.tree.add(lambda cands: 0 if self.event.isData() else self.event.Pileup_nPU(), 'numTrueVertices', 'I')
        self.tree.add(lambda cands: self.event.isData(), 'isData', 'I')
        self.tree.add(lambda cands: 0 if self.event.isData() else self.event.genWeight(), 'genWeight', 'F')
        # scale shifts
        weightMap = {
            0: {'muR':1.0, 'muF':1.0},
            1: {'muR':1.0, 'muF':2.0},
            2: {'muR':1.0, 'muF':0.5},
            3: {'muR':2.0, 'muF':1.0},
            4: {'muR':2.0, 'muF':2.0},
            5: {'muR':2.0, 'muF':0.5},
            6: {'muR':0.5, 'muF':1.0},
            7: {'muR':0.5, 'muF':2.0},
            8: {'muR':0.5, 'muF':0.5},
        }
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[0] if len(self.event.LHEScaleWeight())>0 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[0]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[1] if len(self.event.LHEScaleWeight())>1 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[1]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[2] if len(self.event.LHEScaleWeight())>2 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[2]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[3] if len(self.event.LHEScaleWeight())>3 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[3]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[4] if len(self.event.LHEScaleWeight())>4 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[4]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[5] if len(self.event.LHEScaleWeight())>5 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[5]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[6] if len(self.event.LHEScaleWeight())>6 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[6]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[7] if len(self.event.LHEScaleWeight())>7 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[7]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.LHEScaleWeight()[8] if len(self.event.LHEScaleWeight())>8 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[8]), 'F')

    def __exit__(self, type, value, traceback):
        self.finish()

    def __del__(self):
        self.finish()

    def finish(self):
        print ''
        logging.info('Finishing')
        logging.info('Writing {0} events'.format(self.eventsStored))
        self.outfile.cd()
        cutflowHist = ROOT.TH1D('summedWeights','summedWeights',1,0,1)
        cutflowHist.SetBinContent(1,self.summedWeights)
        self.outfile.Write()
        self.outfile.Close()

    def flush(self):
        sys.stdout.flush()
        sys.stderr.flush()

    def get_event(self):
        return '{run}:{lumi}:{event}'.format(run=self.event.run(),lumi=self.event.lumi(),event=self.event.event())

    #############################
    ### primary analysis loop ###
    #############################
    def analyze(self):
        '''
        The primary analyzer loop.
        '''
        logging.info('Beginning Analysis')
        start = time.time()
        total = 0
        for f, fName in enumerate(self.fileNames):
            if fName.startswith('/store/user'): 
                fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
            elif fName.startswith('/store'):
                fName = '{0}/{1}'.format('root://cmsxrootd.fnal.gov/',fName)
            logging.info('Processing file {0} of {1}: {2}'.format(f+1, len(self.fileNames), fName))
            tfile = ROOT.TFile.Open(fName,'READ')
            tree = tfile.Get(self.treename)
            for row in tree:
                total += 1
                if total==2: start = time.time() # just ignore first event for timing
                if total % 1000 == 1:
                    cur = time.time()
                    elapsed = cur-start
                    remaining = float(elapsed)/total * float(self.totalEntries) - float(elapsed)
                    mins, secs = divmod(int(remaining),60)
                    hours, mins = divmod(mins,60)
                    logging.info('{}: Processing event {}/{} ({} selected) - {}:{:02d}:{:02d} remaining'.format(self.outputTreeName,total,self.totalEntries,self.eventsStored,hours,mins,secs))
                    self.flush()
                self.setupEvent(tree)
                self.perRowAction()
            tfile.Close('R')

    def setupEvent(self,tree):
        '''Setup the event objects'''
        # load objects
        self.event     = Event(tree)
        if self.event.isData(): self.shift = ''
        if not self.event.isData(): self.gen = [GenParticle(tree,entry=i) for i in range(tree.nGenPart)]
        self.electrons = [Electron(tree,entry=i,shift=self.shift) for i in range(tree.nElectron)]
        self.muons     = [Muon(tree,entry=i,shift=self.shift) for i in range(tree.nMuon)]
        self.taus      = [Tau(tree,entry=i,shift=self.shift) for i in range(tree.nTau)]
        self.photons   = [Photon(tree,entry=i,shift=self.shift) for i in range(tree.nPhoton)]
        self.jets      = [Jet(tree,entry=i,shift=self.shift) for i in range(tree.nJet)]
        self.pfmet     = Met(tree,shift=self.shift)

    def perRowAction(self):
        '''Per row action, can be overridden'''
        goodToCheck = self.cutTree.evaluate(self.event)
        if not goodToCheck: return

        # select candidates
        cands = self.selectCandidates()
        cands['event'] = self.event

        # store event?
        goodToStore = self.cutTree.checkCands(cands)
        if not goodToStore: return

        self.tree.fill(cands)
        self.eventsStored += 1

    def selectCandidates(self):
        '''
        Select candidates
            format should be:
            candidates = {
                "objectName" : ("collectionName", position),
                ...
            }
        '''
        logging.warning("You must override selectCandidates.")
        return {}

    #################
    ### utilities ###
    #################
    def getCands(self,coll,func):
        cands = []
        for cand in coll:
            if func(cand): cands += [cand]
        return cands

    def cleanCands(self,src,other,dr):
        cleaned = []
        for s in src:
            keep = True
            for o in other:
                if deltaR(s.eta(),s.phi(),o.eta(),o.phi())<dr:
                    keep = False
                    break
            if keep:
                cleaned += [s]
        return cleaned

    def getCollectionString(self,cand):
        if isinstance(cand,Electron): return 'e'
        elif isinstance(cand,Muon):   return 'm'
        elif isinstance(cand,Tau):    return 't'
        elif isinstance(cand,Photon): return 'g'
        elif isinstance(cand,Jet):    return 'j'
        else:                         return 'a'

    def checkTrigger(self,*datasets,**triggerNames):
        '''Check trigger using trigger map'''
        isData = self.event.isData()>0.5
        # reject triggers if they are in another dataset
        # looks for the dataset name in the filename
        # for MC it accepts any
        reject = True if isData else False
        for dataset in datasets:
            # if we match to the dataset, start accepting triggers
            if dataset in self.fileNames[0] and isData: reject = False
            for trigger in triggerNames[dataset]:
                var = 'HLT_{0}'.format(trigger)
                passTrigger = getattr(self.event,var)()
                if passTrigger>0.5:
                    # it passed the trigger
                    # in data: reject if it corresponds to a higher dataset
                    return False if reject else True
            # dont check the rest of data
            if dataset in self.fileNames[0] and isData: break
        return False

    def metFilter(self):
        if not self.event.isData(): return True

        return self.event.Flag_METFilters()

        # legacy, keep in case we revert
        filterList = [
            'HBHENoiseFilter',
            'HBHENoiseIsoFilter',
            'globalTightHalo2016Filter',
            'EcalDeadCellTriggerPrimitiveFilter',
            'goodVertices',
            'eeBadScFilter',
            'noBadMuons',
            'BadChargedCandidateFilter',
        ]
        notFilterList = [
            # Dont use, "noBadMuons" contains this
            #'duplicateMuons',
            #'badMuons',
        ]
        for f in filterList:
            if getattr(self.event,'Flag_{}'.format(f))()==0:
                logging.info('Rejecting event {0}:{1}:{2} for {3}={4}'.format(self.event.run(), self.event.lumi(), self.event.event(), f, getattr(self.event,f)()))
                self.report_failure('fails {}'.format(f))
                return False
        for f in notFilterList:
            if getattr(self.event,'Flag_{}'.format(f))()>0:
                logging.info('Rejecting event {0}:{1}:{2} for {3}={4}'.format(self.event.run(), self.event.lumi(), self.event.event(), f, getattr(self.event,f)()))
                self.report_failure('fails {}'.format(f))
                return False
        return True

    def report_failure(self,message):
        if self.events:
            event = self.get_event()
            if event in self.events:
                 print event, message

    ##########################
    ### add object to tree ###
    ##########################
    def addTriggers(self):
        pass

    def addCandVar(self,label,varLabel,var,rootType,permissive=False):
        '''Add a variable for a cand'''
        if permissive:
            self.tree.add(lambda cands: getattr(cands[label],var)() if cands[label].has(var) else 0., '{0}_{1}'.format(label,varLabel), rootType)
        else:
            self.tree.add(lambda cands: getattr(cands[label],var)(), '{0}_{1}'.format(label,varLabel), rootType)

    def addFlavorDependentCandVar(self,label,varLabel,varMap,rootType,default=0):
        '''Add a variable for a cand based on flavor'''
        self.tree.add(lambda cands: getattr(cands[label],varMap[cands[label].collName])() if cands[label].collName in varMap else default, '{0}_{1}'.format(label,varLabel), rootType)

    def addMet(self,label):
        '''Add Met variables'''
        self.addCandVar(label,'pt','pt','F')
        self.addCandVar(label,'phi','phi','F')
        #self.addCandVar(label,'cov00','covXX','F')
        #self.addCandVar(label,'cov01','covXY','F')
        #self.addCandVar(label,'cov10','covXY','F')
        #self.addCandVar(label,'cov11','covYY','F')

    def addCandidate(self,label):
        '''Add variables relevant for all objects'''
        self.addCandVar(label,'pt','pt','F')
        self.addCandVar(label,'eta','eta','F')
        self.tree.add(lambda cands: abs(getattr(cands[label],'eta')()), '{}_abseta'.format(label), 'F')
        self.addCandVar(label,'phi','phi','F')
        self.addCandVar(label,'mass','mass','F')

    def addJet(self,label):
        '''Add variables relevant for jets'''
        self.addCandidate(label)
        self.addCandVar(label,'CSVv2','btagCSVV2','F')
        self.addCandVar(label,'hadronFlavour','hadronFlavour','I')

    def addLepton(self,label):
        '''Add variables relevant for leptons'''
        self.addCandidate(label)
        self.addCandVar(label,'charge','charge','I')
        self.addCandVar(label,'dz','dz','F')
        self.addCandVar(label,'dxy','dxy','F')

    def addDetailedMuon(self,label):
        '''Add detailed  variables'''
        pass

    def addDetailedTau(self,label):
        '''Add detailed variables'''
        pass

    def addPhoton(self,label):
        '''Add variables relevant for photons'''
        self.addCandidate(label)
        self.addCandVar(label,'r9','r9','F')

    def genDeltaR(self,cand):
        '''Get the gen level deltaR'''
        if cand.genMatch()==0: return 0.
        eta = cand.eta()
        genEta = cand.genEta()
        phi = cand.phi()
        genPhi = cand.genPhi()
        return deltaR(eta,phi,genEta,genPhi)

    def genJetDeltaR(self,cand):
        '''Get the gen level deltaR'''
        if cand.genJetMatch()==0: return 0.
        eta = cand.eta()
        genEta = cand.genJetEta()
        phi = cand.phi()
        genPhi = cand.genJetPhi()
        return deltaR(eta,phi,genEta,genPhi)

    def addDiCandidate(self,label):
        '''Add variables relevant for a two object candidate'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'M2','M2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaR','deltaR','F')
        self.addCandVar(label,'deltaEta','deltaEta','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')

    def addDiJet(self,label):
        '''Add variables relevant for a dijet candidate'''
        self.addDiCandidate(label)

    def addDiLepton(self,label):
        '''Add variables relevant for a dilepton candidate'''
        self.addDiCandidate(label)

    def addLeptonMet(self,label):
        '''Add variables related to a lepton + met'''
        self.addCandVar(label,'mt','Mt','F')
        self.addCandVar(label,'mt2','Mt2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')

    def addComposite(self,label):
        '''Add variables related to multi object variables'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'M2','M2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')

    def addCompositeMet(self,label):
        '''Add variables related to multi object variables'''
        self.addCandVar(label,'mt','Mt','F')
        self.addCandVar(label,'mt2','Mt2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
