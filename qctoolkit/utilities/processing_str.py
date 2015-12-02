import re
import fileinput

###################################
# Simple text formating functions #
###################################
def delete_next(target, pattern, line_number):
  itr = 0
  matched = False
  for line in fileinput.input(target, inplace=1):
    if pattern in line:
      matched = True
    if matched and itr < line_number and itr > 0:
      itr += 1
    else:
      print line,

def delete(target, pattern, line_number):
  itr = 0
  matched = False
  for line in fileinput.input(target, inplace=1):
    if pattern in line:
      matched = True
    if matched and itr < line_number:
      itr += 1
    else:
      print line,

def insert(target, pattern, text):
  for line in fileinput.input(target, inplace=1):
    print line,
    if pattern in line:
      print text

def replace(target, pattern, text):
  for line in fileinput.input(target, inplace=1):
    print re.sub(re.compile(pattern), text, line),

def containing(target, pattern):
  result = False
  with open(target,"r") as ftarget:
    for line in ftarget:
      if pattern in line:
        result = True
  return result

def matching(target, pattern):
  result = False
  with open(target,'r') as ftarget:
    for line in ftarget:
      if re.match(re.compile(pattern),line):
        result = True
  return result

def fileStrip(path):
  new_path = re.sub('.*/', '', path)
  return new_path

def pathStrip(path):
  new_path = re.sub('//*','/',path)
  new_path = re.sub(r'([^\.])\.\/',r"\1",new_path)
  return new_path

def getPath(path):
  try:
    out = re.match(r'(.*)/', path).group(1)
  except AttributeError:
    out = '.'
  return out
