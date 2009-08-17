import os, sys

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
os.putenv('HOME', '/home/autotest')
os.putenv('PATH', os.environ['PATH']+':'+cwd+'/../../../opencmissextras/cm/external/x86_64-linux-debug-'+compiler+'/bin')
os.system('mpd &')

def testExample(id, path, nodes) :
   global compiler,logDir;
   newDir = logDir
   for folder in path.split('/') :
     newDir = newDir + '/' + folder
     if not os.path.isdir(newDir):
       os.mkdir(newDir)
   os.chdir(path)
   if os.path.exists(newDir + "/test"+id+"-" + compiler) :
     os.remove(newDir + "/test"+id+"-" + compiler)
   execPath='bin/x86_64-linux/'+path.rpartition('/')[2]+'Example-debug-'+compiler
   if nodes == '1' :
     err=os.system(execPath + ' 2 2 0 '+nodes+" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
   else :
     err=os.system('mpiexec -n ' + nodes + " " + execPath + ' 2 2 0 '+nodes+" > " + newDir + "/test" + id + "-" + compiler + " 2>&1")
   if not os.path.exists(execPath) :
     err=-1
   if err==0 :
     print "Testing %s%s: <a class='success' href='%slogs_x86_64-linux/%s/test%s-%s'>successed</a><br>" %(path,id,rootUrl,path,id,compiler)
   else :
     print "Testing %s%s: <a class='fail' href='%slogs_x86_64-linux/%s/test%s-%s'>failed</a><br>" %(path,id,rootUrl,path,id,compiler)
   os.chdir(cwd)
   return;


testExample('1', "ClassicalField/Laplace", '1') 
testExample('2', "ClassicalField/Laplace", '2')

os.system('mpdallexit')

