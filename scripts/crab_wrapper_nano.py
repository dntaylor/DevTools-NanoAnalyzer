#!/usr/bin/env python
import argparse
import sys
import subprocess
import os

# import run script
from DevTools.NanoAnalyzer.runAnalysis import runAnalysis

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit analyzers')

    # hack to deal with crab format
    argv = argv[1:] if argv[0].isdigit() else argv # first one is a job num
    argv = ['--'+a if not a.startswith('-') and '=' in a else a for a in argv]

    parser.add_argument('--analysis', type=str, default='', help='Analysis to submit')
    parser.add_argument('--shift', type=str, nargs='?', default='', help='Energy shift')
    parser.add_argument('--outputFile', type=str, nargs='?', default='crab.root', help='Output filename')

    return parser.parse_args(argv)



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # add on the input and output
    # import the crab generated PSet.py
    argv = []
    if os.path.isfile('PSet.py'):
        import PSet
        fileList = list(PSet.process.source.fileNames)

        # add on the inputfiles
        argv += ['--inputFiles'] + fileList
    else:
        print 'Warning: Failed to find PSet.py'

    argv += ['--outputFile',args.outputFile,'--shift',args.shift]

    status = runAnalysis(args.analysis,argv)

    # write xml file
    with open('FrameworkJobReport.xml','w') as f:
        f.write(
'''<FrameworkJobReport>
<ReadBranches>
</ReadBranches>
<PerformanceReport>
  <PerformanceSummary Metric="StorageStatistics">
    <Metric Name="Parameter-untracked-bool-enabled" Value="true"/>
    <Metric Name="Parameter-untracked-bool-stats" Value="true"/>
    <Metric Name="Parameter-untracked-string-cacheHint" Value="application-only"/>
    <Metric Name="Parameter-untracked-string-readHint" Value="auto-detect"/>
    <Metric Name="ROOT-tfile-read-totalMegabytes" Value="0"/>
    <Metric Name="ROOT-tfile-write-totalMegabytes" Value="0"/>
  </PerformanceSummary>
</PerformanceReport>
<GeneratorInfo>
</GeneratorInfo>
</FrameworkJobReport>'''
        )

    return status

if __name__ == "__main__":
    func = main()
    sys.exit(func)
