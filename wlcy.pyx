cdef extern from "wl.h":
   ctypedef  int* vlamatrix
   int wl(int, int *, int *)

from libc.stdlib cimport *

def wlcy(a):
   cdef int n = a.nrows()
   cdef int *graph
   cdef int i,j, ncells
   graph = <int *>malloc(n*n*(sizeof(int)))
   for i in range(n):
       for j in range(n):
           graph[n*i+j] = a[i,j]
  # rank = wl(n, &ncells, graph)
   rank = wl(n, &ncells, <vlamatrix>graph)
   if rank < 0:
      print "please check your input!\n"
      free(graph)
      return -1, 0
   free(graph)
   nc = ncells
   return rank, nc
