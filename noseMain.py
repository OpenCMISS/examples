import os, subprocess,sys,commands
from time import strftime
from datetime import date
import socket
sys.path.append(os.environ['OPENCMISS_ROOT']+"/cm/examples")
examplesDir = os.environ['OPENCMISS_ROOT']+"/cm/examples"
logsDir = os.environ['OPENCMISS_ROOT']+"/build/logs"
hostname = socket.gethostname()

def load_prop(propFile, properties) :
  properties['TestingPoint']=[]
  for propLine in propFile:
    propDef= propLine.strip()
    if len(propDef) == 0:
        continue
    if propDef[0] in ( '!', '#' ):
        continue
    punctuation= [ propDef.find(c) for c in '= ' ] + [ len(propDef) ]
    found= min( [ pos for pos in punctuation if pos != -1 ] )
    name= propDef[:found].rstrip()
    value= propDef[found:].lstrip("= ").rstrip()
    if value.find(',')==-1 :
      properties[name]= value
    else :
      if (name.startswith('TestingPoint.')) :
        properties['TestingPoint'].append(value.split(','))
      else :
        properties[name]=value.split(',')
  propFile.close()

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

def getVersion(compiler) :
  operating_system = commands.getoutput('uname') 
  if operating_system == 'Linux':
    if compiler == 'gnu':
      version_info = commands.getoutput('gfortran -v')
      if version_info.find('gcc version 4.5') != -1:
        return '_4.5'
      elif version_info.find('gcc version 4.4') != -1 :
        return '_4.4'
  return ''
  

def test_build_library():
  global logsDir
  rootDir=os.environ['OPENCMISS_ROOT']+"/cm"
  os.chdir(rootDir)
  logDir = os.getcwd().replace(rootDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_library_build"
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=gnu" + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_library_build_history",err)
  assert err==0
  
  

def test_example():
  rootdir = os.getcwd()+"/examples/FiniteElasticity/testingPoints"
  for root, subFolders, files in os.walk(rootdir) :
    if root.find(".svn")==-1 :
      for f in files :
        if f=="nightlytest.prop" or (f=="weeklytest.prop" and date.today().weekday() == 6) :
          os.chdir(root)
          yield check_build, 'build',root,"gnu"
          system = os.uname()[0].lower()
          arch = os.uname()[4]
          propFile= file(f, "r")
          properties = dict()
          load_prop(propFile,properties)
          if '42TestingPointsPATH' in properties :
            testingPointsPath = os.environ['OPENCMISS_ROOT']+"/cm/examples/"+properties['42TestingPointsPATH']
          else:
            testingPointsPath = ""
          testpoints = properties['TestingPoint']
          for testpoint in testpoints :
            if (testpoint[0].find("${")!=-1) :
              testingPointPath = properties[testpoint[0][2:testpoint[0].find("}")]]+testpoint[0][testpoint[0].find("}")+1:]
            else:
              testingPointPath = testpoint[0]
            os.chdir(testingPointsPath + testingPointPath)
            if len(testpoint)<=2 :            
              yield check_run, 'run', os.getcwd(), system, arch, "gnu", root, testpoint[1], testingPointPath, False
            else :
              yield check_run, 'run', os.getcwd(), system, arch, "gnu", root, testpoint[1], testingPointPath
              yield check_output,'check',os.getcwd(), testpoint[2], testpoint[3]
  
def check_build(status,root,compiler):
  global examplesDir, logsDir
  logDir = os.getcwd().replace(examplesDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_build"
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=" + compiler + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_build_history",err)
  assert err==0

def check_run(status,cwd, system,arch,compiler,masterPath,testArgs,testPath,noCheck=None):
  global examplesDir, logsDir
  logDir = os.getcwd().replace(examplesDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_run"
  open_log(logPath)
  exampleName = masterPath.rpartition("/")[2]
  compiler_version = getVersion(compiler)
  command = masterPath+"/bin/"+arch+"-"+system+"/mpich2/"+compiler+compiler_version+"/"+exampleName+"Example-debug "+testArgs + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_run_history",err)
  assert err==0

def check_output(status, cwd, ndiffDir, outputDir):
  global examplesDir, logsDir
  logDir = os.getcwd().replace(examplesDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_check"
  open_log(logPath)
  ndiff = os.environ['OPENCMISS_ROOT']+"/cm/utils/ndiff"
  errall =0
  for outputFile in os.listdir(ndiffDir) :
    if outputFile!='.svn' :
      command = ndiff+" --tolerance=1e-10 "+ndiffDir+"/"+outputFile+" "+outputDir+"/"+outputFile + ' >> '  + logPath + " 2>&1"
      err = os.system(command)
      if err!=0 :
        errall = -1
  if errall==0 :
    f1 = open(logPath,"a")
    f1.write("The output values are identical with the expected ones. However, the outputs may be generated from previous build. You need also check if the latest test run passes.")
    f1.close()
  close_log(logPath)
  add_history(logDir+"/nose_check_history",errall)
  assert errall==0 
 
if __name__ == '__main__':
  test_example()   
