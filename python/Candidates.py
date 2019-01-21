import math

import ROOT

from utilities import deltaR, deltaPhi

##########################
### Basic tree wrapper ###
##########################
class Event(object):
    '''
    Wrap a tree.
    '''
    def __init__(self,tree):
        self.tree = tree

    def __getattr__(self,name):
        return lambda: self.get(name) # returns attribute as a function

    def get(self,var):
        return getattr(self.tree,var)

    def has(self,var):
        return hasattr(self.tree,var)

    def isData(self,):
        return not hasattr(self.tree,'genWeight') # hack

    def lumi(self,):
        return self.luminosityBlock()

##############################
### Basic candidate access ###
##############################
class Candidate(object):
    '''
    Encapsulate access to an object in a TTree.
    '''
    def __init__(self,tree,entry=-1,collName=''):
        self.tree = tree
        self.collName = collName
        self.entry = entry

    def __getattr__(self,name):
        return lambda: self.get(name) # returns the attribute as a function

    def __repr__(self):
        return '<{}[{}] {:.2f}:{:.2f}:{:.2f}:{:.2f}>'.format(
            self.collName,
            self.entry,
            self.pt(),
            self.eta(),
            self.phi(),
            self.mass(),
        )

    def get(self,var):
        '''Default variable access from tree.'''
        if not self.tree: return 0
        varName = '{0}_{1}'.format(self.collName,var)
        if self.entry<0: # for flat trees
            return getattr(self.tree,varName)
        else: # for vector branches
            return getattr(self.tree,varName)[self.entry]

    def has(self,var):
        if not self.tree: return False
        varName = '{0}_{1}'.format(self.collName,var)
        return hasattr(self.tree,varName)

    def p4(self):
        p4 = ROOT.TLorentzVector()
        p4.SetPtEtaPhiM(self.pt(), self.eta(), self.phi(), self.mass())
        return p4

############
### Muon ###
############
class Muon(Candidate):
    '''
    Muon object access.
    '''
    def __init__(self,tree,entry=-1,collName='Muon',shift=None):
        super(Muon, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

################
### Electron ###
################
class Electron(Candidate):
    '''
    Electron object access.
    '''
    def __init__(self,tree,entry=-1,collName='Electron',shift=None):
        super(Electron, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

###########
### Tau ###
###########
class Tau(Candidate):
    '''
    Tau object access.
    '''
    def __init__(self,tree,entry=-1,collName='Tau',shift=None):
        super(Tau, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

###########
### Jet ###
###########
class Jet(Candidate):
    '''
    Jet object access.
    '''
    def __init__(self,tree,entry=-1,collName='Jet',shift=None):
        super(Jet, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

##############
### Photon ###
##############
class Photon(Candidate):
    '''
    Photon object access.
    '''
    def __init__(self,tree,entry=-1,collName='Photon',shift=None):
        super(Photon, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

###################
### GenParticle ###
###################
class GenParticle(Candidate):
    '''
    Gen particle object access.
    '''
    def __init__(self,tree,entry=-1,collName='GenPart',shift=None):
        super(GenParticle, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

###########
### MET ###
###########
class Met(Candidate):
    '''
    Met object access.
    '''
    def __init__(self,tree,entry=-1,collName='MET',shift=None):
        super(Met, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

    def __repr__(self):
        return '<{0} {1} {2:.2f}:{3:.2f}>'.format(
            self.collName,
            self.entry,
            self.pt(),
            self.phi(),
        )

    def p4(self):
        metP4 = ROOT.TLorentzVector()
        metP4.SetPtEtaPhiM(self.pt(),0.,self.phi(),0)
        return metP4

############################
### Composite candidates ###
############################
class CompositeCandidate(object):
    '''
    Primary object for access to composite variables.
    '''
    def __init__(self,*objects):
        self.objects = objects

    def __getitem__(self,key):
        if isinstance(key,int): return self.objects[key]

    def __setitem__(self,key,value):
        if isinstance(key,int): self.objects[key] = value

    def __getattr__(self,name):
        try:
            return self.get(name)
        except:
            return self.get(name.capitalize()) # catch things like 'pt' instead of 'Pt'

    def __repr__(self):
        return '<CompositeCandidate {0}>'.format(
            ' '.join([obj.__repr__() for obj in self.objects])
        )

    def get(self,var):
        '''Default variable access from TLorentzVector'''
        vec = self.p4()
        return getattr(vec,var)

    def p4(self):
        p4 = ROOT.TLorentzVector()
        for obj in self.objects: p4 += obj.p4()
        return p4

    def st(self):
        return sum([obj.pt() for obj in self.objects])

###################
### Dicandidate ###
###################
class DiCandidate(CompositeCandidate):
    '''
    Dicandidate variable access.
    '''
    def __init__(self,obj0,obj1):
        super(DiCandidate, self).__init__(obj0,obj1)

    def deltaR(self):
        return deltaR(self.objects[0].eta(),
                      self.objects[0].phi(),
                      self.objects[1].eta(),
                      self.objects[1].phi())

    def deltaPhi(self):
        return deltaPhi(self.objects[0].phi(),
                       self.objects[1].phi())

    def deltaEta(self):
        return abs(self.objects[0].eta()-self.objects[1].eta())

##########################
### Candidate plus met ###
##########################
class MetCompositeCandidate(CompositeCandidate):
    '''
    Met + candidate variable specialization.
    '''
    def __init__(self,met,*objects):
        if not isinstance(met,Met):
            raise TypeError('First argument must be Met object')
        super(MetCompositeCandidate, self).__init__(met,*objects)
        self.met = met
        self.cands = objects

    def metP4(self):
        return self.met.p4()

    def candP4(self):
        p4 = ROOT.TLorentzVector()
        for cand in self.cands: p4 += cand.p4()
        return p4

    def mt(self):
        metP4 = self.metP4()
        candP4 = self.candP4()
        return math.sqrt(abs((candP4.Et()+metP4.Et())**2 - ((candP4+metP4).Pt())**2))

    def mt2(self):
        return self.mt()**2

    def Mt(self):
        return self.mt()

    def Mt2(self):
        return self.mt2()

    def deltaPhi(self):
        candP4 = self.candP4()
        return deltaPhi(self.met.phi(),candP4.Phi())

    def _x(self,c,o):
        dphi = deltaPhi(c.phi(),self.met.phi())
        codphi = deltaPhi(c.phi(),o.phi())
        if abs(codphi)==math.pi or codphi==0:
            x = 0
        else:
            x = c.pt()/(c.pt() + self.met.pt()*(math.cos(dphi) - math.sin(dphi)/math.tan(codphi)))
        # note: large met and small deltaphi between c/o results in large negative values for denom
        # this doesnt work for our boosted topology in a->tt since met resolution is poor
        #if x<0:
        #    print x, c.pt(), self.met.pt(), dphi, codphi
        return x

    def mcat(self,i1=0,i2=1):
        c1 = self.cands[i1]
        c2 = self.cands[i2]
        x1 = self._x(c1,c2)
        x2 = self._x(c2,c1)
        if x1*x2<=0: return 0.
        cp4 = ROOT.TLorentzVector()
        cp4 += self.cands[i1].p4()
        cp4 += self.cands[i2].p4()
        if len(self.cands)==2: # both candidates had neutrinos
            mcoll = cp4.M() / math.sqrt(x1*x2)
        else: # add in the other cands (assuming no contribution from met)
            op4 = ROOT.TLorentzVector()
            for i,c in enumerate(self.cands):
                if i in [i1,i2]: continue
                op4 += c.p4()
            ncp4 = ROOT.TLorentzVector()
            ncp4.SetPtEtaPhiM(cp4.Pt()/(x1*x2), cp4.Eta(), cp4.Phi(), cp4.M() / math.sqrt(x1*x2))
            tcp4 = ncp4+op4
            mcoll = tcp4.M()
        return mcoll

    def Mcat(self,*args):
        return self.mcat(*args)
