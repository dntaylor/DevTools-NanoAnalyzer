# CutTree.py
# Hold the pass fail decision of events
import sys
import logging
sys.argv.append('-b')
import ROOT
sys.argv.pop()
from array import array

class CutTree(object):
    '''
    Stores pass-fail decisions of events
    '''
    def __init__(self):
        self.labels = []
        self.selections = {}
        self.results = {}
        #self.tree = ROOT.TTree("CutTree","CutTree")
        #self.results['run'] = array('i',[0])
        #self.results['lumi'] = array('i',[0])
        #self.results['event'] = array('L',[0])
        #self.tree.Branch('run',self.results['run'],'run/I')
        #self.tree.Branch('lumi',self.results['lumi'],'lumi/I')
        #self.tree.Branch('event',self.results['event'],'event/l')
        self.filled = set()

    def add(self, fun, label):
        if label not in self.results:
            self.labels += [label]
            self.selections[label] = fun
            #self.results[label] = array('i',[0])
            self.results[label] = 0
            #self.tree.Branch(label,self.results[label],'{0}/I'.format(label))
        else:
            logging.warning("{0} already in CutTree.".format(label))

    def getLabels(self):
        return self.labels

    def evaluate(self,cands):
        eventkey = '{0}:{1}:{2}'.format(cands['event'].run(), cands['event'].lumi(), cands['event'].event())
        if eventkey in self.filled:
            logging.warning("Event {0} already filled.".format(eventkey))
            return False
        else:
            passAll = True
            # verify each cut
            for label in self.selections:
                cut = self.selections[label]
                self.results[label] = int(cut(cands))
                if not self.results[label]: passAll = False
            # verify we have a candidate
            if passAll:
                for cname,cand in cands.iteritems():
                    if cand==None:
                        passAll = False
            return passAll

    def getResults(self):
        return self.results
