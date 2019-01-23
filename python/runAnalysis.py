#!/usr/bin/env python

# import run script
from DevTools.NanoAnalyzer.MonoHZZAnalysis import main as runMonoHZZ
from DevTools.NanoAnalyzer.MonoHZZFakeRateAnalysis import main as runMonoHZZFakeRate

analyses = {
    'MonoHZZ': runMonoHZZ,
    'MonoHZZFakeRate': runMonoHZZFakeRate,
}

def runAnalysis(analysis,argv):
    '''Return analysis function'''
    func = analyses[analysis]

    return func(argv)

