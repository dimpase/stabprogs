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

#include <stdlib.h>
#include  <stdio.h>
#include "wl.h"

int real_main(int vert, FILE *f)
{
int graph[vert][vert];
int rank, v;
int i,j, ncells;
long time; float t;
long start_time,end_time;
for (i=0;i<vert;i++) for(j=0;j<vert;j++) {
   fscanf (f, "%d", &v);
   graph[i][j] = v; }
start_time=clock()/1000;
rank = wl(vert, &ncells, graph);
end_time=clock()/1000-start_time;
if (rank<0) {
   printf("please check your input!\n"); 
return;
}
printf ("\n\n number of colors: %6d",rank);
printf ("\n\n number of cells: %6d", ncells);
/* printf("\n\n%ld msec \n\n",end_time); */
prmat(vert, graph);
}

int main(int narg, char *arg[10])
{
FILE *f;
int rank,vert;
f=fopen(arg[1],"r");        /* char */
if (f == NULL) { printf("\n file does not exist\n"); exit(0);}
fscanf (f,"%d%d",&rank,&vert);
return real_main(vert, f);
}
