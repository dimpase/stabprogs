/* An implementation of the Weisfeiler-Leman graph stabilization procedure 
    see http://arxiv.org/abs/1002.1921 and references therein for more info.
    Complete sources can be found at http://www.ntu.edu.sg/home/dima/software.htm

    Copyright (C) 1989-2010 Luitpold Babel, Dmitrii V. Pasechnik, and others 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
    Authors can be contacted at
    e-mail: Luitpold Babel  <luitpold.babel@unibw.de>
    e-mail: Dmitrii Pasechnik <dimpase@gmail.com>                          */

NB!!! c2.wl is an example of 0-1 basis for a commutative unital algebra
of dimension 10 which is not c.c. - it's not closed under transposition.
None of the 3 versions on the program compute stabilization correctly.

TO BE FIXED!!! (note added 2013-10-01)

(this file is outdated!)
PROGRAM IMPLEMENTATIONS OF THE WEISFEILER-LEMAN ALGORITHM

   You have access to two different implementations:
       1. the program   stabil
       2. the program   stcol

   For detailed descriptions of both programs we refer
   to the technical report by
       L. Babel, I.V. Chuvaeva, M. Klin, D.V. Pasechnik:
         'Algebraic Combinatorics in Mathematical Chemistry.
                   Methods and Algorithms.
       II. Program Implementation of the Weisfeiler-Leman Algorithm'

   The programs compute the cellular algebra W which is generated
   by the adjacency matrix of any colored directed or undirected graph Gamma.
 
   Both programs are written in programming language C
   and will work on Unix computers. 

   In the current versions, the programs stabil and stcol will handle 
   graphs with up to 200 respectively 150 vertices.

INPUT:
   
   Both programs require as an input a file containing the
   following information about the colored graph Gamma:
       1. Number of colors
       2. Number of vertices
       3. Adjacency matrix of Gamma

   Note that the colors have to be denoted consecutively by 0,1,2,...
   Furthermore, the colors in the diagonal of the adjacency matrix 
   have to be different from the other colors.
   
   Feasibility of the input is checked within the programs in 
   a preprocessing routine.

OUTPUT:

   Both programs provide as an output the following information
   about the computed cellular algebra W:
       1. Number of colors
       2. Number of cells
       3. Adjacency matrix of W


HOW TO RUN THE PROGRAM:
 
       Start the program with the command
               stabil input_file                     resp.
                stcol input_file 


EXAMPLE:
 
     1. Create a file called input_file 
        (or choose any other name for this file)
        which contains data like the following:


              4
              8
              3  1  2  1  1  2  2  2
              1  0  1  2  2  1  2  2
              2  1  3  1  2  2  1  2
              1  2  1  0  2  2  2  1 
              1  2  2  2  0  1  2  1
              2  1  2  2  1  3  1  2
              2  2  1  2  2  1  0  1
              2  2  2  1  1  2  1  3

     
     2. Run the program using the command:


              stabil input_file 


     3. The output on the screen is:


              number of colors: 8

              number of cells:  2

              adjacency matrix of the cellular algebra:

              0   2   3   2   2   3   5   3 
              4   1   4   6   6   4   6   7 
              3   2   0   2   5   3   2   3 
              4   6   4   1   6   7   6   4 
              4   6   7   6   1   4   6   4 
              3   2   3   5   2   0   2   3 
              7   6   4   6   6   4   1   4 
              3   5   3   2   2   3   2   0 



    _____________________________________________________

