import os, datetime

for examplePath, subFolders, files in os.walk(".") :
   for f in files :
     now = datetime.datetime.now()
     dur = datetime.timedelta(days=31)
     history = now - dur
     if history > datetime.datetime.fromtimestamp(os.path.getmtime(examplePath+"/"+f)) \
          and f.startswith("nightly_") and (not "history" in f) :
               os.remove(examplePath+"/"+f)
