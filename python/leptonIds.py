

def passHZZLooseMuonNoIso(cand):
    cuts = {
        'pt' : lambda cand: cand.pt()>5,
        'eta': lambda cand: abs(cand.eta())<2.4,
        'dxy': lambda cand: abs(cand.dxy())<0.5,
        'dz' : lambda cand: abs(cand.dz())<1,
        # muon best track type not in nanoaod?
        'id' : lambda cand: (cand.isGlobal() or (cand.isTracker() and cand.nStations()>0)), # and cand.muonBestTrackType!=2,
        'sip': lambda cand: abs(cand.sip3d())<4,
    }
    # TODO ghost cleaning
    return all([cuts[cut](cand) for cut in cuts])

def passHZZLooseMuon(cand):
    cuts = {
        'loose': lambda cand: passHZZLooseMuonNoIso(cand),
        'iso'  : lambda cand: cand.pfRelIso03_all()<0.35,
    }
    return all([cuts[cut](cand) for cut in cuts])

def passHZZTightMuon(cand):
    cuts = {
        'loose': lambda cand: passHZZLooseMuon(cand),
        'tight': lambda cand: cand.isPFcand() or (cand.highPtId()>0 and cand.pt()>200),
    }
    return all([cuts[cut](cand) for cut in cuts])

def passHZZLooseElectronNoIso(cand):
    cuts = {
        'pt' : lambda cand: cand.pt()>7,
        'eta': lambda cand: abs(cand.eta())<2.5,
        'dxy': lambda cand: abs(cand.dxy())<0.5,
        'dz' : lambda cand: abs(cand.dz())<1,
        'sip': lambda cand: abs(cand.sip3d())<4,
    }
    return all([cuts[cut](cand) for cut in cuts])

def passHZZLooseElectron(cand):
    cuts = {
        'loose': lambda cand: passHZZLooseElectronNoIso(cand),
        'iso'  : lambda cand: cand.pfRelIso03_all()<0.35,
    }
    return all([cuts[cut](cand) for cut in cuts])

def passHZZTightElectron(cand):
    pt = cand.pt()
    eta = abs(cand.eta() + cand.deltaEtaSC())
    mva = cand.mvaFall17V2Iso()
    cuts = {
        'loose' : lambda cand: passHZZLooseElectron(cand),
        'id' : lambda cand: (pt<=10 and ((eta<0.8                and mva>0.5739521065342641)\
                                    or   (eta>=0.8 and eta<1.479 and mva>0.5504628790992929)\
                                    or   (eta>=1.479             and mva>0.5924627534389098)))\
                            or\
                            (pt>10  and ((eta<0.8                and mva>-0.03391387993354392)\
                                    or   (eta>=0.8 and eta<1.479 and mva>-0.018451958064666783)\
                                    or   (eta>=1.479             and mva>-0.38565459150737535))),
    }
    return all([cuts[cut](cand) for cut in cuts])

def passHZZLooseNoIso(cand):
    if cand.collName=='Electron': return passHZZLooseElectronNoIso(cand)
    if cand.collName=='Muon'    : return passHZZLooseMuonNoIso(cand)
    return False

def passHZZLoose(cand):
    if cand.collName=='Electron': return passHZZLooseElectron(cand)
    if cand.collName=='Muon'    : return passHZZLooseMuon(cand)
    return False

def passHZZTight(cand):
    if cand.collName=='Electron': return passHZZTightElectron(cand)
    if cand.collName=='Muon'    : return passHZZTightMuon(cand)
    return False
