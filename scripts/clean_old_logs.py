import os, datetime

for examplePath, subFolders, files in os.walk(".") :
   for f in files :
     now = datetime.datetime.now()
     dur = datetime.timedelta(days=31)
     history = now - dur
     print os.path.getmtime(examplePath+"/"+f)
     if history > datetime.datetime.fromtimestamp(os.path.getmtime(examplePath+"/"+f)) :
             os.remove(examplePath+"/"+f)
