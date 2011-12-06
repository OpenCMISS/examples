import os, subprocess,sys,commands
from time import strftime
from datetime import date
import socket
import json


### Utility Methods and Global Variables ###

sys.path.append(os.environ['OPENCMISSEXAMPLES_ROOT'])
hostname = socket.gethostname()

def append_path(prefix, suffix) :
  if (not prefix.endswith("/")) :
    prefix = prefix + "/"
  return prefix + suffix

globalExamplesDir = os.environ['OPENCMISSEXAMPLES_ROOT'] 
examplesDir = globalExamplesDir if (not 'DIR' in os.environ) else append_path(globalExamplesDir,os.environ['DIR'])
logsDir = append_path(os.environ['OPENCMISS_ROOT'],"build/logs") 

compiler = 'gnu' if (not 'COMPILER' in os.environ) else os.environ['COMPILER']
def get_version(compiler) :
  operating_system = commands.getoutput('uname')
  if operating_system == 'Linux':
    if compiler == 'gnu':
      version_info = commands.getoutput('gfortran -v')
      version_location = version_info.find('gcc version ')
      if version_location != -1:
        return '_' + version_info[version_location+12:version_location+15]
  return ''

compiler_version = compiler + get_version(compiler) if (not 'OPENCMISS_COMPILER_PATH' in os.environ) else os.environ['OPENCMISS_COMPILER_PATH']

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

class Test:
  def __init__(self, dct, example,globalExamplesDir):
    self.args = "" 
    self.processors = 1 
    if ("expectedPath" in dct):
      self.outputPath = "." 
      self.tolerance = 1e-7 
    self.__dict__.update(dct)
    self.path = example.path if ("path" not in dct) else append_path(append_path(globalExamplesDir,example.globalTestDir),dct["path"])

  def __repr__(self):
    return "processors: %d, args: %s" %(self.processors,self.args)

class Example:
  
  def __init__(self,dct,compiler_version,globalExamplesDir,examplesDir,logsDir):
    self.path = os.getcwd()
    self.system = os.uname()[0].lower()
    self.arch = os.uname()[4]
    self.compilerVersion = compiler_version
    self.logPath = self.path.replace(examplesDir,logsDir)
    self.globalTestDir = "" if ("globalTestDir" not in dct) else dct["globalTestDir"]
    test_dct = dct["test"]
    self.tests = []
    for test_entry in test_dct :
      self.addTest(Test(test_entry,self,globalExamplesDir))

  def addTest(self, test):
    self.tests.append(test)

  def __repr__(self):
    return self.path

        
def object_encode(dct) :
  global compiler_version,globalExamplesDir,examplesDir,logsDir
  example = Example(dct["example"],compiler_version,globalExamplesDir,examplesDir,logsDir)
  return example

def load_log_dir(examplePath) :
  global logsDir
  logDir = examplePath.replace(os.environ['OPENCMISSEXAMPLES_ROOT'],logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  return logDir
  


size = 'small' if (not 'SIZE' in os.environ) else os.environ['SIZE']

  
############ Beginning of tests #######################

def test_build_library():
  yield check_build_library,compiler_version,os.uname()[4],os.uname()[0].lower()
   
def check_build_library(compiler_version,arch,system):
  global logsDir, compiler
  rootDir=os.environ['OPENCMISS_ROOT']+"/cm"
  os.chdir(rootDir)
  logDir = os.getcwd().replace(rootDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_library_build_" + compiler_version + str(date.today())
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=%s >> %s 2>&1" %(compiler,logPath)
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_library_build_history_" + compiler_version,err)
  assert err==0
  
def test_example():
  global examplesDir
  for examplePath, subFolders, files in os.walk(examplesDir) :
    if examplePath.find(".svn")==-1 :
      for f in files :
        if (size=='small' and f=='nightlytest.json') or (size=='large' and f in ('nightlytest.json','weeklytest.json')) :
          os.chdir(examplePath)
          json_data=open(f).read()
          example = object_encode(json.loads(json_data))
          yield check_build, 'build',example
          for test in example.tests : 
            yield check_run, 'run', example, test
            if (hasattr(test, 'expectedPath')):
              yield check_output,'check', example, test
  
def check_build(status,example):
  global compiler
  logDir = load_log_dir(example.path)
  logPath = append_path(logDir,"nose_build_" + example.compilerVersion + str(date.today()))
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=%s >> %s 2>&1" %(compiler,logPath)
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_build_history_" + example.compilerVersion,err)
  assert err==0

def check_run(status,example,test):
  logDir = load_log_dir(test.path)
  os.chdir(test.path)
  logPath = append_path(logDir,"nose_run_" + example.compilerVersion + str(date.today()))
  open_log(logPath)
  exampleName = example.path.rpartition("/")[2]
  command = "%s/bin/%s-%s/mpich2/%s/%sExample-debug %s >> %s 2>&1" %(example.path,example.arch,example.system,example.compilerVersion,exampleName,test.args,logPath)
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_run_history_" + example.compilerVersion,err)
  assert err==0

def check_output(status,example,test):
  logDir = load_log_dir(test.path)
  os.chdir(test.path)
  logPath = append_path(logDir,"nose_check_" + example.compilerVersion + str(date.today()))
  open_log(logPath)
  ndiff = os.environ['OPENCMISS_ROOT']+"/cm/utils/ndiff"
  errall =0
  for outputFile in os.listdir(test.expectedPath) :
    if outputFile!='.svn' :
      command = "%s --tolerance=%f %s/%s %s/%s >> %s 2>&1" %(ndiff,test.tolerance,test.expectedPath,outputFile,test.outputPath,outputFile,logPath)
      err = os.system(command)
      if err!=0 :
        errall = -1
  if errall==0 :
    f1 = open(logPath,"a")
    f1.write("The output values are identical with the expected ones. However, the outputs may be generated from previous build. You need also check if the latest test run passes.")
    f1.close()
  close_log(logPath)
  add_history(logDir+"/nose_check_history_" + example.compilerVersion,errall)
  assert errall==0 
 
if __name__ == '__main__':
  test_example()   
