
from pprint import pprint
from simulation import simulate
import json
import sys

materialParameters = [1.0, 5.0, 10.0, 5.0]
resultsFileName = "results.json"
if len(sys.argv) >= 5:
    materialParameters = [float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])]
    
if len(sys.argv) == 6:
    resultsFileName = sys.argv[5]

resultsFibre = simulate(0.0, materialParameters)

# by re-orienting the fibre direction in the unit cube we can simulate a cross fibre extension
resultsCross = simulate(90.0, materialParameters)

results = {}
results["materialParameters"] = materialParameters
results["fibre"] = resultsFibre
results["cross"] = resultsCross

pprint(resultsFibre)
pprint(resultsCross)
pprint(results)

with open(resultsFileName, 'w') as rjson:
    json.dump(results, rjson)
