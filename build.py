import os, mmap, re, sys
from jinja2 import Template, Environment, FileSystemLoader
from datetime import date	
from time import strftime
import json
import fileinput

env = Environment(loader=FileSystemLoader('.'))
template = env.get_template('build.template')
nesiTemplate = env.get_template('nesi1.template')
globalExamplesDir = os.environ['OPENCMISSEXAMPLES_ROOT']
size = 'small' if (not 'SIZE' in os.environ) else os.environ['SIZE']
testSets = ["nightlytest.json","weeklytest.json"] if (size == 'large') else ["nightlytest.json"]
examplesDir = globalExamplesDir if (not 'DIR' in os.environ) else "%s/%s" %(globalExamplesDir,os.environ['DIR'])
rootLogDir = "%s/%s" %(os.environ['OPENCMISS_ROOT'],"build/logs") 
machine = os.environ['MACHINE'] if ('MACHINE' in os.environ) else None
mpi = os.environ['OPENCMISS_MPI_PATH']
masterLogDir = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_%s" %(machine if machine != None else os.environ['archname'])
compiler = 'gnu' if (not 'COMPILER' in os.environ) else os.environ['COMPILER']
compilerVersion = 'gnu_4.6' if (not 'OPENCMISS_COMPILER_PATH' in os.environ) else os.environ['OPENCMISS_COMPILER_PATH']
system = os.uname()[0].lower()
arch = os.uname()[4]

class Class:

  def __init__(self,name,path=None,parent=None):
    self.name = name
    if path==None :
      self.path = "%s/%s" %(parent.path,name)
      parent.addChild(self)
    else :
      self.path = path
      self.parent = None
    self.children = []
    self.fail = 0

  def addChild(self, child):
    self.children.append(child)
    child.parent = self

  def findChild(self, name):
    for child in self.children :
      if child.name == name :
        return child
    return None
  
  def accumulateParentFail(self):
    parent = self.parent
    if (parent!=None) :
      parent.fail = parent.fail+1
      parent.accumulateParentFail() 

  def __repr__(self):
    return self.path

class Example(Class):
  
  def __init__(self,name,dct,parent):
    Class.__init__(self,name=name,parent=parent)
    self.logDir = self.path.replace(globalExamplesDir,rootLogDir)
    self.masterLogDir = self.path.replace(globalExamplesDir,masterLogDir)
    self.ensureDir(self.logDir) 
    self.language = None
    self.tests = []
    self.buildFail = 0
    if dct!=None :
      self.globalTestDir = "" if ("globalTestDir" not in dct) else dct["globalTestDir"]
      self.language = None if  ("language" not in dct) else dct["language"]
      self.script = None if ("script" not in dct) else dct["script"]
      test_dct = dct["test"]
      for test_entry in test_dct :
        if not ((machine=="nesi") ^ ("machine" in test_entry and test_entry["machine"] == "nesi")) :
          self.addTest(Test(test_entry,self))

  def ensureDir(self,path) :
    if not os.path.exists(path):
      self.ensureDir(path[:path.rindex("/")])
      os.makedirs(path)

  def addTest(self, test):
    self.tests.append(test)

  def start(self) :
    if self.language == None :
      self.build()
    if self.buildFail==0 :
      for test in self.tests :
        test.run()
        if test.runFail == 0 and hasattr(test, 'expectedPath'):
          test.check()

  def build(self) :
    cwd = os.getcwd()
    os.chdir(self.path)
    logPath = "%s/nightly_build_%s_%s.log" %(self.logDir,compilerVersion,str(date.today()))
    self.wrapWithPre(logPath,1)
    command = "make >> %s 2>&1" %(logPath)
    self.buildFail = os.system(command)
    self.wrapWithPre(logPath,0)
    self.buildLog = "%s/nightly_build_%s_%s.log" %(self.masterLogDir,compilerVersion,str(date.today()))
    self.buildHistoryLog = "%s/nightly_build_history_%s.log" %(self.masterLogDir, compilerVersion)
    self.buildHistory = self.add_history("%s/nightly_build_history_%s.log" %(self.logDir, compilerVersion),self.buildFail)   
    if self.buildFail != 0 :
      self.fail = 1
      self.accumulateParentFail()
    os.chdir(cwd)


  def invalidConfig(self) :
    self.buildFail = 1
    logPath = "%s/nightly_build_%s_%s.log" %(self.logDir,compilerVersion,str(date.today()))
    f1 = open(logPath,"w")
    f1.write("Invalid JSON configuration.")
    f1.close()
    self.buildLog = "%s/nightly_build_%s_%s.log" %(self.masterLogDir,compilerVersion,str(date.today()))
    self.buildHistoryLog = "%s/nightly_build_history_%s.log" %(self.masterLogDir,compilerVersion)
    self.buildHistory = self.add_history("%s/nightly_build_history_%s.log" %(self.logDir,compilerVersion),self.buildFail)   
    if self.buildFail != 0 :
      self.fail = 1
      self.accumulateParentFail()

  def add_history(self,path,fail) :
    if os.path.exists(path) :
      history = open(path,"a")
    else :
      history = open(path,"w")
      history.write("Completed Time&ensp;Status<br>\n")
    if fail==0 :
      history.write(strftime("%Y-%m-%d %H:%M:%S")+'&ensp;<a class="success">success</a><br>\n')
    else :
      history.write(strftime("%Y-%m-%d %H:%M:%S")+'&ensp;<a class="fail">fail</a><br>\n')
    history.close()
    history = open(path,"r")
    return self.tail(history)


  def tail(self,f,window=5):
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if (bytes - BUFSIZ > 0):
            # Seek back one whole BUFSIZ
            f.seek(block*BUFSIZ, 2)
            # read BUFFER
            data.append(f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.append(f.read(bytes))
        linesFound = data[-1].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return '\n'.join(''.join(data).splitlines()[-window:])

  def __repr__(self):
    return self.path

  def wrapWithPre(self,path,openTag=1) :
    if openTag== 1 :
      f1 = open(path,"w")
      f1.write("<pre>")
    else :
      f1 = open(path,"a")
      f1.write("</pre>")
    f1.close()

class Test(Example):
  def __init__(self, dct, example):
    self.id = dct["id"]
    self.args = "" 
    self.processors = 1 
    self.machine = None
    if ("expectedPath" in dct):
      self.outputPath = "." 
      self.tolerance = 1e-7 
    self.__dict__.update(dct)
    self.path = example.path if ("path" not in dct) else "%s/%s/%s" %(globalExamplesDir,example.globalTestDir,dct["path"])
    self.logDir = example.logDir
    self.exampleName = example.name
    self.masterLogDir = example.masterLogDir
    self.parent = example

  def run(self) :
    cwd = os.getcwd()
    os.chdir(self.path)
    logPath = "%s/nightly_run_%d_%s_%s.log" %(self.logDir,self.id,compilerVersion,str(date.today()))
    self.wrapWithPre(logPath,1)
    if self.parent.language == "python" :
      command = "python %s %s > %s 2>&1" %(self.parent.script, self.args,logPath)
    elif self.machine == "nesi" :
      self.command = "%s/bin/%s-%s/%s/%s/%sExample-debug %s" %(self.parent.path,arch,system,mpi,compilerVersion,self.exampleName,self.args)
      f = open("nesi_%d.ll" %(self.id),"w")
      f.write(nesiTemplate.render(test=self))
      f.close()
      command = "llsubmit -s nesi_%d.ll > %s 2>&1" %(self.id,logPath)
    else :
      command = "%s/bin/%s-%s/%s/%s/%sExample-debug %s > %s 2>&1" %(self.parent.path,arch,system,mpi,compilerVersion,self.exampleName,self.args,logPath)
    self.runFail = os.system(command)
    self.wrapWithPre(logPath,0)
    self.runLog = "%s/nightly_run_%d_%s_%s.log" %(self.masterLogDir,self.id,compilerVersion,str(date.today()))
    self.runHistoryLog = "%s/nightly_run_history_%d_%s.log" %(self.masterLogDir,self.id,compilerVersion)  
    if self.runFail != 0 :
      self.fail = 1
      self.accumulateParentFail() 
    elif self.machine == "nesi" :
      # Find the output log and replace with the submission log
      size = os.stat(logPath).st_size
      f = open(logPath, "r")
      data = mmap.mmap(f.fileno(), size, access=mmap.ACCESS_READ)
      m = re.search(r'\.[0-9]+', data)
      f.close()   
      f1 = open("nesi%s.out" %(m.group(0)), "r")
      output = f1.read()
      f1.close()
      self.wrapWithPre(logPath,1)
      f = open(logPath,"a")
      if output.find("ERROR")>0 :
        self.runFail=1
      f.write(output)
      f.close()
      self.wrapWithPre(logPath,0)
    self.runHistory = self.add_history("%s/nightly_run_history_%d_%s.log" %(self.logDir,self.id,compilerVersion),self.runFail)  
    os.chdir(cwd)


  def check(self):
    cwd = os.getcwd()
    os.chdir(self.path)
    logPath = "%s/nightly_check_%d_%s_%s" %(self.logDir,self.id,compilerVersion,str(date.today()))
    self.wrapWithPre(logPath,1)
    self.checkFail = 0
    ndiff = os.environ['OPENCMISS_ROOT']+"/cm/utils/ndiff"
    for outputFile in os.listdir(self.expectedPath) :
      if outputFile!='.svn' :
        command = "%s --tolerance=%e %s/%s %s/%s >> %s 2>&1" %(ndiff,self.tolerance,self.expectedPath,outputFile,self.outputPath,outputFile,logPath)
        checkFail = os.system(command)
        if checkFail!=0 :
          self.checkFail = 1
    self.wrapWithPre(logPath,0)
    self.checkLog = "%s/nightly_check_%d_%s_%s.log" %(self.masterLogDir,self.id,compilerVersion,str(date.today()))
    self.checkHistoryLog = "%s/nightly_check_history_%d_%s.log" %(self.masterLogDir,self.id,compilerVersion)
    self.checkHistory = self.add_history("%s/nightly_check_history_%d_%s.log" %(self.logDir,self.id,compilerVersion),self.checkFail)  
    if self.checkFail!=0 :
      self.fail = 1
      self.accumulateParentFail()

  def __repr__(self):
    return "processors: %d, args: %s" %(self.processors,self.args)



def object_encode(name,parent,dct=None) :
  global compilerVersion,globalExamplesDir,examplesDir
  if dct==None :
    example = Example(name=name,parent=parent)
  else :
    example = Example(name,dct["example"],parent)
  return example

def fileInTestSets(f,path) :
  global testSets
  for t in testSets :
    if f == t :
      if machine == "nesi" :
        linestring = open("%s/%s" %(path,f), 'r').read()
        if linestring.find("nesi")!=-1:
          return True
      else :
        return True
  return False

root = Class(name="examples", path=examplesDir)
for path, subFolders, files in os.walk(top=root.path,topdown=True) :
  if path.find(".svn")==-1 :	
    for f in files :
      if fileInTestSets(f,path) :
        pathFromRoot = path[len(root.path)+1:path.rfind('/')]
        parent = root
        if len(pathFromRoot.strip()) != 0 :
          for dirToPath in pathFromRoot.split("/") :
            t = parent.findChild(dirToPath)
            if t==None :
              t = Class(name=dirToPath,parent=parent)
            parent = t
        # Example
        os.chdir(path)
        try:
          json_data=open(f).read()
          example = object_encode(name=path[path.rfind('/')+1:],parent=parent,dct=json.loads(json_data))
          example.start()
        except ValueError:
          example = Example(name=path[path.rfind('/')+1:],parent=parent,dct=None)
          example.invalidConfig()
        

os.chdir(globalExamplesDir)

print template.render(examples=root)
if root.fail != 0 :
  exit("ERROR: At least one examples failed")
