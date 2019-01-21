# AnalysisTree.py
import sys
import logging
sys.argv.append('-b')
import ROOT
sys.argv.pop()
from array import array


typeMap = {
    'I': int,
    'l': long,
    'F': float,
    'C': str,
}
arrayMap = {
    'I': 'i',
    'l': 'L',
    'F': 'f',
}

class AnalysisTree(object):
    '''
    Wrapper for the analysis tree
    '''
    def __init__(self,name):
        self.tree = ROOT.TTree(name,name)
        self.filled = set()
        self.branches = {}
        self.branchSet = set()
        self.branches['run']   = {'var': array('i',[0]), 'rootType': 'I', 'function': 'run', 'branchName': 'run'}
        self.branches['lumi']  = {'var': array('i',[0]), 'rootType': 'I', 'function': 'luminosityBlock', 'branchName': 'lumi'}
        self.branches['event'] = {'var': array('L',[0]), 'rootType': 'l', 'function': 'event', 'branchName': 'event'}
        self.__addBranch('run')
        self.__addBranch('lumi')
        self.__addBranch('event')

    def __addBranch(self,label):
        if label not in self.branchSet:
            self.tree.Branch(label, self.branches[label]['var'], '{0}/{1}'.format(self.branches[label]['branchName'],self.branches[label]['rootType']))
            self.branchSet.add(label)
        else:
            logging.error('Branch with label "{0}" already exists.'.format(label))

    def add(self, fun, label, rootType):
        if label not in self.branches:
            if rootType[0]=='C': # special handling of string
                self.branches[label] = {'var': bytearray(rootType[1]), 'rootType': rootType[0], 'function': fun, 'branchName': '{0}[{1}]'.format(label,rootType[1]), 'size': rootType[1]}
            else:
                self.branches[label] = {'var': array(arrayMap[rootType],[0]), 'rootType': rootType, 'function': fun, 'branchName': label}
            self.__addBranch(label)
        else:
            logging.error("{0} already in AnalysisTree.".format(label))
            raise ValueError("{0} already in AnalysisTree.".format(label))

    def __evaluate(self,label,cands):
        pyType = typeMap[self.branches[label]['rootType']]
        if isinstance(self.branches[label]['function'], basestring): # its just a branch in the tree
            self.branches[label]['var'][0] = pyType(getattr(cands['event'],self.branches[label]['function'])())
        else:
            if self.branches[label]['rootType']=='C': # special handling of string
                strSize = self.branches[label]['size']
                funcVal = pyType(self.branches[label]['function'](cands))
                #if len(funcVal)==strSize-1:
                if len(funcVal)<strSize:
                    self.branches[label]['var'][:strSize] = funcVal
                else:
                    logging.error('Size mismatch function with label {0}.'.format(label))
            else:
                self.branches[label]['var'][0] = pyType(self.branches[label]['function'](cands))

    def fill(self,cands,**kwargs):
        allowDuplicates = kwargs.pop('allowDuplicates',False)
        eventkey = '{0}:{1}:{2}'.format(cands['event'].run(), cands['event'].lumi(), cands['event'].event())
        if eventkey in self.filled and not allowDuplicates:
            logging.warning("Event {0} already filled.".format(eventkey))
        else:
            for label in self.branches:
                try:
                    self.__evaluate(label,cands)
                except:
                    logging.error('Error processing branch {0}'.format(label))
                    raise
            self.tree.Fill()
            self.filled.add(eventkey)
