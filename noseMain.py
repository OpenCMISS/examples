import os, subprocess,sys,commands
from time import strftime
from datetime import date
import socket
sys.path.append(os.environ['OPENCMISSEXAMPLES_ROOT'])
examplesDir = os.environ['OPENCMISSEXAMPLES_ROOT']
logsDir = os.environ['OPENCMISS_ROOT']+"/build/logs"
hostname = socket.gethostname()
compiler = os.environ['COMPILER']
size = os.environ['SIZE']
parentdir = os.environ['DIR']

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
      version_location = version_info.find('gcc version ')
      if version_location != -1:
        return '_' + version_info[version_location+12:version_location+15]
  return ''


def test_build_library():
  global logsDir, compiler
  rootDir=os.environ['OPENCMISS_ROOT']+"/cm"
  os.chdir(rootDir)
  logDir = os.getcwd().replace(rootDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_library_build_" + compiler + str(date.today())
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=" + compiler + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_library_build_history_" + compiler,err)
  assert err==0
  
def test_example():
  global compiler,examplesDir,parentdir
  if (not examplesDir.endswith("/")) :
     rootdir = examplesDir + "/"
  rootdir = rootdir+parentdir
  for root, subFolders, files in os.walk(rootdir) :
    if root.find(".svn")==-1 :
      for f in files :
        if (size=='small' and f=='nightlytest.prop') or (size=='large' and (f=='nightlytest.prop' or f=='weeklytest.prop')) :
          os.chdir(root)
          yield check_build, 'build',root,compiler
          system = os.uname()[0].lower()
          arch = os.uname()[4]
          propFile= file(f, "r")
          properties = dict()
          load_prop(propFile,properties)
          if '42TestingPointsPATH' in properties :
            testingPointsPath = examplesDir+"/"+properties['42TestingPointsPATH']
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
              yield check_run, 'run', os.getcwd(), system, arch, compiler, root, testpoint[1], testingPointPath, False
            else :
              yield check_run, 'run', os.getcwd(), system, arch, compiler, root, testpoint[1], testingPointPath
              if len(testpoint)<=4 :
                yield check_output,'check',os.getcwd(), testpoint[2], testpoint[3]
              else :
                yield check_output,'check',os.getcwd(), testpoint[2], testpoint[3], testpoint[4]
  
def check_build(status,root,compiler):
  global examplesDir, logsDir
  logDir = os.getcwd().replace(examplesDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_build_" + compiler + str(date.today())
  open_log(logPath)
  command = "make USEFIELDML=true COMPILER=" + compiler + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_build_history_" + compiler,err)
  assert err==0

def check_run(status,cwd, system,arch,compiler,masterPath,testArgs,testPath,noCheck=None):
  global examplesDir, logsDir
  logDir = os.getcwd().replace(examplesDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_run_" + compiler + str(date.today())
  open_log(logPath)
  exampleName = masterPath.rpartition("/")[2]
  compiler_version = getVersion(compiler)
  command = masterPath+"/bin/"+arch+"-"+system+"/mpich2/"+compiler+compiler_version+"/"+exampleName+"Example-debug "+testArgs + " >> "  + logPath + " 2>&1"
  err = os.system(command)
  close_log(logPath)
  add_history(logDir+"/nose_run_history_" + compiler,err)
  assert err==0

def check_output(status, cwd, ndiffDir, outputDir,tolerance=None):
  global examplesDir, logsDir
  logDir = os.getcwd().replace(examplesDir,logsDir)
  newDir = ''
  for folder in logDir.split('/') :
    newDir = newDir + '/' + folder
    if not os.path.isdir(newDir):
      os.mkdir(newDir)
  logPath = logDir+"/nose_check_" + compiler + str(date.today())
  open_log(logPath)
  ndiff = os.environ['OPENCMISS_ROOT']+"/cm/utils/ndiff"
  errall =0
  for outputFile in os.listdir(ndiffDir) :
    if outputFile!='.svn' :
      if tolerance != None :
        command = ndiff+" --tolerance=" + tolerance +" "+ndiffDir+"/"+outputFile+" "+outputDir+"/"+outputFile + ' >> '  + logPath + " 2>&1"
      else :
        command = ndiff+" "+ndiffDir+"/"+outputFile+" "+outputDir+"/"+outputFile + ' >> '  + logPath + " 2>&1" 
      err = os.system(command)
      if err!=0 :
        errall = -1
  if errall==0 :
    f1 = open(logPath,"a")
    f1.write("The output values are identical with the expected ones. However, the outputs may be generated from previous build. You need also check if the latest test run passes.")
    f1.close()
  close_log(logPath)
  add_history(logDir+"/nose_check_history_" + compiler,errall)
  assert errall==0 
 
if __name__ == '__main__':
  test_example()   
