# here we get the current path of this script so that we can specify things relative to
# this script, rather than the execution folder. And the modules path...
import sys, os
from os.path import dirname, join
ScriptPath = dirname(__file__)
RunScript = os.path.abspath(os.path.join(ScriptPath, "run.py"))
Root = os.getcwd()

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import matplotlib.pyplot as plt
import json
from subprocess import run, DEVNULL
from pprint import pprint

colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colourIndex = 0

parameterSets = [
                 ["1.0", "5.0", "10.0", "5.0"], #defaults
                 ["0.25", "5.0", "10.0", "5.0"],
                 ["0.5", "5.0", "10.0", "5.0"],
                 ["0.75", "5.0", "10.0", "5.0"],
                 ["1.25", "5.0", "10.0", "5.0"],
                 ["1.5", "5.0", "10.0", "5.0"],
                 ["1.75", "5.0", "10.0", "5.0"],
                 ["2.0", "5.0", "10.0", "5.0"]
                 ]

for i, parameters in enumerate(parameterSets, 0):
    resultsFileName = os.path.abspath(os.path.join(Root, "results-{:03d}.json".format(i+1)))
    run(["python", RunScript] + parameters + [resultsFileName], stdout=DEVNULL)
    with open(resultsFileName) as resultsFile:
        results = json.load(resultsFile)
    labelText = str(results["materialParameters"])
    plt.plot(results["fibre"]["strain"], results["fibre"]["stress"], colours[colourIndex] + "-", label=labelText)
    plt.plot(results["cross"]["strain"], results["cross"]["stress"], colours[colourIndex] + "1--")
    colourIndex = colourIndex + 1
    if colourIndex >= len(colours):
        colourIndex = 0

plt.ylim(ymax=50)
plt.ylabel('axial force')
plt.xlabel('axial strain')
plt.legend(loc=0)
plt.show()
