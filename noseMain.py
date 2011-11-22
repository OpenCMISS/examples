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

examplesDir = os.environ['OPENCMISSEXAMPLES_ROOT'] if (not 'DIR' in os.environ) else append_path(os.environ['OPENCMISSEXAMPLES_ROOT'],os.environ['DIR'])
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

class Struct:
    def __init__(self, dct):
      self.__dict__.update(dct)
        
          

def object_encode(dct) :
  return Struct(dct)




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
  yield check_build_library,compiler_version
   
def check_build_library(compiler_version):
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
  command = "make USEFIELDML=true COMPILER=" + compiler + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_library_build_history_" + compiler_version,err)
  assert err==0
  
def test_example():
  global compiler_version,examplesDir,logsDir
  for examplePath, subFolders, files in os.walk(examplesDir) :
    if examplePath.find(".svn")==-1 :
      for f in files :
        if (size=='small' and f=='nightlytest.json') or (size=='large' and (f=='nightlytest.json' or f=='weeklytest.json')) :
          json_data=open(examplePath + "/" + f).read()
          example = json.loads(json_data,object_hook=object_encode).example
          example.path = examplePath
          example.system = os.uname()[0].lower()
          example.arch = os.uname()[4]
          example.compilerVersion = compiler_version
          example.logPath = example.path.replace(examplesDir,logsDir)
          yield check_build, 'build',example
          for test in example.test : 
            if (not hasattr(test,"path")) :
              test.path = example.path        
            yield check_run, 'run', example, test
            if (hasattr(test, 'expectedPath')):
              yield check_output,'check', example, test
  
def check_build(status,example):
  global compiler
  os.chdir(example.path)
  logDir = load_log_dir(example.path)
  logPath = append_path(logDir,"nose_build_" + example.compilerVersion + str(date.today()))
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=" + compiler + " >> "  + logPath + " 2>&1"
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
  command = example.path+"/bin/"+example.arch+"-"+example.system+"/mpich2/"+example.compilerVersion+"/"+exampleName+"Example-debug "+test.args + " >> "  + logPath + " 2>&1"
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
      command = ndiff+" --tolerance=" + test.tolerance +" "+test.expectedPath+"/"+outputFile+" "+test.outputPath+"/"+outputFile + ' >> '  + logPath + " 2>&1"
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
