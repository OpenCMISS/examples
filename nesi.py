import os
from jinja2 import Template, Environment, FileSystemLoader
from datetime import date	
from time import strftime

env = Environment(loader=FileSystemLoader('.'))
template = env.get_template('nesi.template')
rootLogDir = "%s/../Logs/examples" %(os.getcwd())
masterLogDir = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_nesi/examples"

class Example:

  def __init__(self,name,path=None,parent=None):
    self.name = name
    if path==None :
      self.path = "%s/%s" %(parent.path,name)
      parent.addChild(self)
    else :
      self.path = path
      self.parent = None
    self.logDir = "%s/%s" %(rootLogDir,self.path)
    self.masterLogDir = "%s/%s" %(masterLogDir, self.path)
    self.ensureDir(self.logDir) 
    self.children = []
    self.fail = 0

  def ensureDir(self,path) :
    if not os.path.exists(path):
      self.ensureDir(path[:path.rindex("/")])
      os.makedirs(path)

  
  def addChild(self, child):
    self.children.append(child)
    child.parent = self

  def build(self) :
    cwd = os.getcwd()
    os.chdir(self.path)
    logPath = "%s/nesi_build_%s.log" %(self.logDir,str(date.today()))
    command = "make MPI_DIR=/usr/mpi/gcc/openmpi-1.6 >> %s 2>&1" %(logPath)
    self.buildFail = os.system(command)
    self.buildLog = "%s/nesi_build_%s.log" %(self.masterLogDir,str(date.today()))
    self.buildHistory = self.add_history("%s/nesi_build_history.log" %(self.logDir),self.buildFail)   
    if self.buildFail != 0 :
      self.fail = 1
      parent = self.parent
      while (parent!=None) :
        parent.fail = parent.fail+1
        parent = parent.parent   
    os.chdir(cwd)

  def run(self) :
    cwd = os.getcwd()
    os.chdir(self.path)
    logPath = "%s/nesi_run_%s.log" %(self.logDir,str(date.today()))
    command = "llsubmit nesi.ll %s 2>&1" %(logPath)
    self.runFail = os.system(command)
    self.runLog = "%s/nesi_run_%s.log" %(self.masterLogDir,str(date.today()))
    self.runHistory = self.add_history("%s/nesi_run_history.log" %(self.logDir),self.runFail)
    if self.runFail != 0 :
      self.fail = 1
      parent = self.parent
      while (parent!=None) :
        parent.fail = parent.fail+1
        parent = parent.parent
    os.chdir(cwd)

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

root = Example(name="examples", path=".")
for path, subFolders, files in os.walk(top=root.path,topdown=True) :
  if path.find(".svn")==-1 :	
    for f in files :
      if (f=='nesi.ll') :
        for dirToPath in path.split("/") :
          if dirToPath == root.path :
            parent = root 
          else :
            example = Example(name=dirToPath,parent=parent)
            parent = example
        example.isLeaf = True
        example.build()
        if example.buildFail==0 :
          example.run()

htmlFile = open("aa.html","w")
htmlFile.write(template.render(examples=root))
htmlFile.close()
