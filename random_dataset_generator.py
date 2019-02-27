

import numpy as np
import random


#from samplingMethods.py import RNAsubopt
import sys,os

def RNAsubopt(rna,n):
  L   = [] #list of sec str sampled
  D   = {} # dictionary to allow easy removal of duplicates
  cmd = 'echo %s | RNAsubopt -p %d' % (rna,n)
  file = os.popen(cmd)
  line = file.readline()
  while line:
    secstr    = line.strip()
    D[secstr] = 1
    line = file.readline()
  file.close()
  L = D.keys()
  return L


#from clote_misc import basePairList
import sys,os,tempfile,string,copy


def basePairList(secStr):
  L = []; stack = []
  for i in range(len(secStr)):
    if secStr[i]=="(":
      stack.append(i)
    elif secStr[i]==")":
      x = stack.pop() #return and pop end of stack 
      L.append( (x,i) ) 
  return L




def randomSeqGenerator(length, seed):
    random.seed(seed)
    base = ["A","C","G","U"]
    seq = ""
    for i in range(length):
        seq += base[random.randint(0,3)]
    return seq

def isoLenGenerator(length, dist_lower, dist_upper, n):
    result = [set([]) for i in range(dist_lower, dist_upper+1)]
    seed_count = 0
    succeed = 0
    while succeed < n * (dist_upper - dist_lower + 1):
        seq = randomSeqGenerator(length, seed_count)
        seed_count += 1
        structures = set(RNAsubopt(seq, 2))
        structures.discard(seq)
        if len(structures) != 2: continue
        structures = list(structures)
        bp_set0 = set(basePairList(structures[0]))
        bp_set1 = set(basePairList(structures[1]))
        dist = len(bp_set0 ^ bp_set1)
        if dist < dist_lower or dist_upper < dist: continue
        if len(result[dist - dist_lower]) < n:
            data = (seq, structures[0], structures[1])
            if data in result[dist - dist_lower]: continue
            result[dist - dist_lower].add(data)
            succeed += 1
            #print(str(succeed)+"/"+str(n * (dist_upper - dist_lower + 1)))
    return [list(i) for i in result]

def isoDistGenerator(dist, len_lower, len_upper, step, n):
    result = [set([]) for i in range(len_lower, len_upper+1, step)]
    seed_count = 0
    succeed = 0
    for length in range(len_lower, len_upper+1, step):
        while len(result[(length - len_lower) // step]) < n:
            seq = randomSeqGenerator(length, seed_count)
            seed_count += 1
            structures = set(RNAsubopt(seq, 2))
            structures.discard(seq)
            if len(structures) != 2: continue
            structures = list(structures)
            bp_set0 = set(basePairList(structures[0]))
            bp_set1 = set(basePairList(structures[1]))
            if dist != len(bp_set0 ^ bp_set1): continue
            data = (seq, structures[0], structures[1])
            if data in result[(length - len_lower) // step]: continue
            result[(length - len_lower) // step].add(data)
            succeed += 1
            #print(str(succeed)+"/"+str(n * ((len_upper - len_lower + 1) // step)))

    return [list(i) for i in result]

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("error")
        sys.exit()
    if sys.argv[1] == "isolen":
        result = isoLenGenerator(100,5,20,100)
        for x in result:
            for y in x:
                print(y[0]+" "+y[1]+" "+y[2])
    else:
        result = isoDistGenerator(15,50,150,10,100)
        for x in result:
            for y in x:
                print(y[0]+" "+y[1]+" "+y[2])