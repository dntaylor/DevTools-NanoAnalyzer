#!/usr/bin/env python

# import run script
from DevTools.NanoAnalyzer.MonoHZZAnalysis import main as runMonoHZZ

analyses = {
    'MonoHZZ': runMonoHZZ,
}

def runAnalysis(analysis,argv):
    '''Return analysis function'''
    func = analyses[analysis]

    return func(argv)

