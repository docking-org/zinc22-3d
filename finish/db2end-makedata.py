#!/usr/bin/env python

#ryan g. coleman
#this script reads a bunch of db2.gz files specified by the commandline
#it puts them into a directory of directory of directories, etc.
#and edits a plaintext file that tells you where things are located

import string
import sys
import os
import shutil
import glob
import itertools
import gzip
import bz2

db2file = "db.db2.bz2"
outTempDir = 'chembldb2'
outTempDatabase = 'chembldb2.txt'
inDirs = sys.argv[1:]
preDir = os.getcwd()

def writeData(lines, thisid, outDir=outTempDir, datafile=None):
  '''writes these lines to the database in outDir'''
  location = [thisid[count:count + 2] for count in xrange(0, len(thisid), 2)]
  locationDir = os.path.join(outDir, *location[:-1])
  #print thisid, location, locationDir  # uncomment for super slow debugging
  try:
    os.makedirs(locationDir)
  except OSError:
    pass  # this means the location directory exists, which is fine
  outGzName = os.path.join(locationDir, thisid + '.db2.gz')
  if os.path.exists(outGzName):  # means we need to append
    outGz = gzip.GzipFile(outGzName, 'a')
  else:  # otherwise just open new file
    outGz = gzip.GzipFile(outGzName, 'w')
    #also write to this file if first time
    if datafile is not None:  # write the chembl ID & location/filename
      datafile.write(thisid + '\t' + os.path.join(preDir, outGzName) + '\n')
  for line in lines:
    outGz.write(line)
  outGz.close()

outDatabase = open(outTempDatabase, 'a')
for inDir in inDirs:
  for oneDb2file in itertools.chain(
        glob.iglob(os.path.join(inDir, "TEMP*.mol2", db2file)),
        glob.iglob(os.path.join(inDir, "*", "TEMP*.mol2", db2file))): 
    read2file = bz2.BZ2File(oneDb2file, 'r')
    recordLines = []
    currentId = None
    for line in read2file:
      if (line.find('M ') == 0) and (len(recordLines) == 0):  # first M line ID
        currentId = string.split(line, sep=None, maxsplit=2)[1]
        recordLines.append(line)
      elif line.find('M ') > 0:  # this is bad, just dump it. see comment at end
        position = line.find('M ')
        recordLines = [line[position:]] #dump everything before this, start over
        currentId = string.split(recordLines[0], sep=None, maxsplit=2)[1]
      elif line[0] == 'E': #end of one record
        recordLines.append(line)
        writeData(recordLines, currentId, datafile=outDatabase)
        recordLines = []  # reset
        currentId = None  # reset
      else:
        recordLines.append(line)  # default is just extend
    read2file.close()
outDatabase.close()

#due to stupid cluster issues, old code produces these
#X      1003  30    351   +6.8084   -3.1632   -7.8126
# X      1004  31    351   +5.7896   -2.640M     ZINC52473197        98
#find them and remove them

