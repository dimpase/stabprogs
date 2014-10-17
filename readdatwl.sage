from wlcy import *
def fromdat(fname): # handle .dat format (2 matrices) from http://anjos.mgi.polymtl.ca/qaplib
   with open(fname) as f:
   #   co = int(f.readline())
      n  = int(f.readline())
      a = matrix([[int(x) for x in f.readline().split()] for i in range(n)])
      f.readline()
      b = matrix([[int(x) for x in f.readline().split()] for i in range(n)])
      return wlcy(a), wlcy(b)

def fromwl(fname): # handle stabil/stabcol format
   with open(fname) as f:
      co = int(f.readline())
      n  = int(f.readline())
      a = matrix([[int(x) for x in f.readline().split()] for i in range(n)])
      return wlcy(a)

# instructions: 
# 1) build libwl.so using makefile (just do make)
# 2) build wlcy.so from sage -sh prompt by running
#    python setup.py build_ext --inplace 
# 3) start sage from sage -sh prompt with 
#   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. sage
# 4) load this file into sage, and run fromdat("blah.dat"), etc

