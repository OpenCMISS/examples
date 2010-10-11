import os,subprocess

def test_example():
  rootdir = os.getcwd()
  for root, subFolders, files in os.walk(rootdir) :
    if root.find(".svn")==-1 :
      for file in files :
        if file=="nightlytest.py" :
          yield check_example, root

def check_example(path) :
  temp = os.getcwd()
  os.chdir(path)
  execCommand = subprocess.Popen(args=["nosetests","nightlytest.py"])
  err = os.waitpid(execCommand.pid, 0)[1]
  assert err==0
 
if __name__ == '__main__':
  test_example()   
