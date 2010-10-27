import os, subprocess,sys
sys.path.append(os.environ['OPENCMISS_ROOT']+"/cm/examples")

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

def test_example():
  rootdir = os.getcwd()
  for root, subFolders, files in os.walk(rootdir) :
    if root.find(".svn")==-1 :
      for f in files :
        if f=="nightlytest.prop" :
          os.chdir(root)
          yield check_build, 'build',root,"gnu"
          system = os.uname()[0].lower()
          arch = os.uname()[4]
          propFile= file(f, "r")
          properties = dict()
          load_prop(propFile,properties)
          testingPointsPath = os.environ['OPENCMISS_ROOT']+"/cm/examples/42TestingPoints/"+properties['42TestingPointsPATH']
          testpoints = properties['TestingPoint']
          for testpoint in testpoints :
            os.chdir(testingPointsPath + testpoint[0])
            if len(testpoint)<=2 :            
              yield check_run, 'run', system, arch, "gnu", root, testpoint[1], testpoint[0], False
            else :
              yield check_run, 'run', system, arch, "gnu", root, testpoint[1], testpoint[0]
              yield check_output,'check',testpoint[2], testpoint[3]
  
def check_build(status,root,compiler):
  f1 = open("temp0","w")
  command = "make COMPILER=" + compiler
  execCommand = subprocess.Popen(args=command.split(), stdout=f1,stderr=f1)
  f1.close()
  err = os.waitpid(execCommand.pid, 0)[1]
  assert err==0

def check_run(status,system,arch,compiler,masterPath,testArgs,testPath,noCheck=None):
  f1 = open("temp1","w")
  exampleName = masterPath.rpartition("/")[2]
  command = masterPath+"/bin/"+arch+"-"+system+"/mpich2/"+compiler+"/"+exampleName+"Example-debug "+testArgs
  execCommand = subprocess.Popen(args=command.split(), stdout=f1,stderr=f1)
  f1.close()
  err = os.waitpid(execCommand.pid, 0)[1]
  assert err==0

def check_output(status,ndiffDir, outputDir):
  ndiff = os.environ['OPENCMISS_ROOT']+"/cm/utils/ndiff"
  for outputFile in os.listdir(ndiffDir) :
    if outputFile!='.svn' :
      command = ndiff+" --tolerance=1e-10 "+ndiffDir+"/"+outputFile+" "+outputDir+"/"+outputFile
      err = os.system(command)
      assert err==0 
 
if __name__ == '__main__':
  test_example()   
