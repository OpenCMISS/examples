import os,subprocess
  
def check_build(system,compiler,masterPath):
  f1 = open("temp0","w")
  temp = os.getcwd()
  os.chdir(masterPath)
  if system=="linux" :
    command = "make COMPILER=" + compiler
  execCommand = subprocess.Popen(args=command.split(), stdout=f1,stderr=f1)
  f1.close()
  err = os.waitpid(execCommand.pid, 0)[1]
  assert err==0
  os.chdir(temp)

def check_run(system,arch,compiler,masterPath,testArgs):
  f1 = open("temp1","w")
  exampleName = masterPath.rpartition("/")[2]
  command = masterPath+"/bin/"+arch+"-"+system+"/mpich2/"+compiler+"/"+exampleName+"Example-debug "+testArgs
  execCommand = subprocess.Popen(args=command.split(), stdout=f1,stderr=f1)
  f1.close()
  err = os.waitpid(execCommand.pid, 0)[1]
  assert err==0

def check_output(ndiffDir, outputDir):
  ndiff = os.environ['OPENCMISS_ROOT']+"/cm/utils/ndiff"
  for outputFile in os.listdir(ndiffDir) :
    if outputFile!='.svn' :
      command = ndiff+" --tolerance=1e-10 "+ndiffDir+"/"+outputFile+" "+outputDir+"/"+outputFile
      err = os.system(command)
      assert err==0 
