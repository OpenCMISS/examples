import os
from time import strftime
from datetime import date
import socket

examplesList = ["ClassicalField/Laplace/Laplace/Fortran","FiniteElasticity/Cantilever/Fortran"]
exampleDir = os.getcwd()
hostname = socket.gethostname()
mastersite = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_nesi"

def open_log(logPath) :
  f1 = open(logPath,"w")
  f1.write("<pre>")
  f1.close()

def close_log(logPath) :
  f1 = open(logPath,"a")
  f1.write("</pre>")
  f1.close()

def add_history(historyPath,err) :
  global hostname
  if os.path.exists(historyPath) :
    history = open(historyPath,"a")
  else :
    history = open(historyPath,"w")
    history.write("<pre>Completed Time\t\tStatus\tHostname\n")
  if err==0 :
    history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tsuccess\t'+hostname+'\n')
  else :
    history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tfail\t'+hostname+'\n')
  history.close()

for example in examplesList :
  os.chdir(example)
  logDir = "%s/../logs/examples/%s" %(exampleDir,example)
  logPath = "%s/nesi_build_%s.log" %(logDir,str(date.today()))
  open_log(logPath)
  command = "make MPI_DIR=/usr/mpi/gcc/openmpi-1.6 >> %s 2>&1" %(logPath)
  err = os.system(command)
  close_log(logPath)
  add_history("%s/nesi_build_history.log" %(logDir),err)
  if err==0:
    print "%s example builds successfully." %(example)
    command = "llsubmit nesi.ll"
    err = os.system(command)
    if err==0 :
       print "Testing job for %s example has been successfully submitted. Please refer to the number above for the job id. You may go to %s/examples/%s/joblogs and select the relevant job log. Please note, since the job execution is queued, it does not guarantee that the job is executed when the logs are uploaded. It is suggested to view the job logs one day after the job is submitted." %(example, mastersite,example)
    else :
      print "Unable to submit the job for %s example" %(example)
  else :
    print "Building %s example failed. Please refer to %s/examples/%s/nesi_build_%s.log for detail." %(mastersite,example,str(date.today()))
  os.chdir(exampleDir)
