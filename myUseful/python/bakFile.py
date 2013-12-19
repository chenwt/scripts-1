#!/usr/bin/python
#J.HE
#Desp: back up all RUN.sh file under giving root directory
#input: <director name> 
        #<backup director folder> 
        #<search file name> 
        #<backup tile name>
#output: ??? 
#Usage: ----

import os,getopt
import sys
import shutil

def createBackupPath(dirList):
    backupPath = bakDir
    for path in dirList:
        backupPath = backupPath + "/" + path
    return backupPath
    
def createScanPath(dirList):
    scanPath = rootDir
    for path in dirList:
        scanPath = scanPath + "/" + path
    return scanPath

def findFile(dirList):
    cwd = os.getcwd()
    print cwd
    crtFileList  = os.listdir(cwd)
    if not crtFileList:
        return
    for filename in crtFileList:
        if (filename == oldFile and not os.path.isdir(filename)):
            backupDirPath = createBackupPath(dirList)
            if not os.path.exists(backupDirPath):
                try:    
                    os.makedirs(backupDirPath)
                except IOError as e:
                    print "Cannot create folder! " + e
                except OSError as e:
                    print "Cannot copy file!" + e
                except:
                    print "Unknown error"
                finally:
                    pass
            backupFilePath = backupDirPath + "/" + bakFile
            try:
                shutil.copyfile(os.getcwd() + "/" + filename, backupFilePath)
            except IOError as e:
                print "Cannot copy file!" + str(e)
            except OSError as e:
                print "Cannot copy file!" + str(e)
            except:
                print "Unknown error"
            finally:
                pass
        elif (os.path.isdir(filename) and not os.path.islink(filename)):
            dirList.append(filename)
            try:
                os.chdir(dirList[-1])
                findFile(dirList)
                dirList.pop()
                os.chdir(createScanPath(dirList))         
            except OSError as e:
                print "Cannot access to the path!" + str(e)
                dirList.pop()
            finally:
                pass 
            #del dirList[-1]
            
    return

argv = sys.argv[1:]
usage = "Usage : \n python bakFile.py -r <root directory> -b <backup dir> -s <search file> -t <back up file>"
try:
  opts,args = getopt.getopt(argv,"hr:b:s:t:")
except getopt.GetoptError:
  print usage
  sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print usage
    sys.exit()
  elif opt in ("-r"):
    rootDir      = arg
  elif opt in ("-b"):
    bakDir   = arg
  elif opt in ("-s"):
    oldFile  = arg	
  elif opt in ("-t"):
    bakFile  =  arg  
          							
if rootDir == '' or bakDir == '' or oldFile == '' or bakFile == '' :
  print "input parameter ERROR!"
  print usage
  sys.exit(2)

print "backup dir:\t" + rootDir +"\t to \t" + bakDir
print "backup file:\t" + oldFile + "\t as \t" + bakFile

os.chdir(rootDir)
dirList = []
findFile(dirList)

print "#-----END------"



    
    

