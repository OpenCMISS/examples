import os, sys

cwd = os.getcwd();
success = 1;
logDir = cwd + "/../../build/logs";
rootUrl = "https://autotest.bioeng.auckland.ac.nz/opencmiss-admin/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
os.system('mpd.py &')

def testExample(id, path, nodes) :
   global success,compiler,logDir;
   os.putenv('PATH', os.environ['PATH']+':'+cwd+'/../../../opencmissextras/cm/external/x86_64-linux-debug-'+compiler+'/bin')
   newDir = logDir
   for folder in path.split('/') :
     newDir = newDir + '/' + folder
     if not os.path.isdir(newDir):
       os.mkdir(newDir)
   os.chdir(path)
   if os.path.exists(newDir + "/test"+id+"-" + compiler) :
     os.remove(newDir + "/test"+id+"-" + compiler)
   execPath='bin/x86_64-linux/'+path.rpartition('/')[2]+'Example-debug-'+compiler
   err=os.system('mpiexec.py -n ' + nodes + " " + execPath + ' 10 10 0 '+nodes+" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
   if not os.path.exists(execPath) :
     err=-1
   if err==0 :
     print "Testing " + path + id + ': success'
   else :
     success=0
     print "Testing" + path + id + ': fail'
   os.chdir(cwd)
   return;


testExample('1', "ClassicalField/Laplace", '1') 
#testExample('2', "ClassicalField/Laplace", '2')

print "See %slogs_x86_64-linux for detail" % (rootUrl) 
if success==0 :
  raise RuntimeError


