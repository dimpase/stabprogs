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
#define memlength (3*(vert*vert))

FILE *f;
struct triple
{int col; int val; struct triple *ptr;};
struct edge
{ int row; int col; struct edge *ptr;};

static struct edge *space;
static struct triple **lines;
static int *memory;
static struct triple *cnst;

int find_value( int *values, int len, int v) {
    int i;
    i=0;
    while(i<len) {
        if(values[i]==v) break;
        i++;
    }
    return i;
}

int standardize(int *graph, int ord) {
    int nv,nd;
    int i,j,k;
    int *values, *dvalues;
 
    values=malloc(ord*ord*sizeof(int));
    dvalues=malloc(ord*sizeof(int));
    nv=0;
    nd=0;

    for(i=0;i<ord;i++) {
        k=find_value(dvalues,nd,graph[i*ord+i]);
        if(k==nd)dvalues[nd++]=graph[i*ord+i];
    }

    for(i=0;i<ord;i++) for(j=0;j<ord;j++) {
        if(i==j) {
            k=find_value(dvalues,nd,graph[i*ord+j]);
            graph[i*ord+j]=k;
        } else {
            k=find_value(values,nv,graph[i*ord+j]);
            if(k==nv)values[nv++]=graph[i*ord+j];
            graph[i*ord+j]=k+nd;
        }
    }
    free(dvalues);
    free(values);
    return nv+nd;
}

int antisymmetrize(int *graph, int ord) {
    int i,j;
    for(i=0;i<ord;i++) for(j=i+1;j<ord;j++) {
        int t, tp;
        t = graph[i*ord+j];
        tp= graph[j*ord+i];
        graph[i*ord+j]+=65536*tp;
        graph[j*ord+i]+=65536*t;
    }
    return 0;
}

int main(int narg, char *arg[10])
{
struct edge **color;
int i,j,rank,vert;
int *graph;
long time; float t;
long start_time,end_time;
f=fopen(arg[1],"r");        /* char */
if (f == NULL) { printf("\n file does not exist\n"); exit(0);}
fscanf (f,"%d%d",&rank,&vert);
graph=malloc(vert*vert*sizeof(int));
lines=malloc(vert*vert*sizeof(struct triple *));
color=malloc(vert*vert*sizeof(struct edge *));
space=malloc(vert*vert*sizeof(struct edge));
memory=malloc(memlength*sizeof(int));
cnst=malloc(vert*sizeof(struct triple));
for (i=0;i<vert;i++) for(j=0;j<vert;j++) fscanf (f, "%d", &graph[i*vert+j]);
antisymmetrize(graph, vert);
rank=standardize(graph,vert);
i=edgepack (graph,rank,vert,color); 
if(i==0) {printf("please check your input!\n"); return;}
printf ("\n\n number of colors: ");

start_time=clock()/1000;
stabil (&rank,vert,graph,color);
end_time=clock()/1000-start_time;

printf("\b\b\b\b\b\b%6d",rank);
for(i=0; i<rank; i++) color[i]->row=-1; j=0;
for(i=0; i<vert; ++i) {
 if(color[graph[i*vert+i]]->row<0) color[graph[i*vert+i]]->row=j++;};
printf ("\n\n number of cells: %6d", j);
for(i=0; i<rank; ++i) if(color[i]->row<0) color[i]->row=j++;
printf("\n\n adjacency matrix of the cellular algebra:\n\n");
 for(i=0; i<vert; ++i) { 
     for(j=0; j<vert; ++j) { 
       printf("%4d ",graph[i*vert+j]);
     };
     printf("\n");
 };
/* printf("\n\n%ld msec \n\n",end_time); */
}

int edgepack (graph,rank,vert,color) 
int *graph, vert, rank;
struct edge **color;
{
int k,i,j;
struct edge *free;

 for(i=0; i<vert; ++i) for(j=0; j<vert; ++j) {
  if(graph[i*vert+j]<0 || graph[i*vert+j]>rank-1) return(0);};
 for(free=space; free<space+rank; free++) free->row=0; 
 for(i=0; i<vert; ++i) space[graph[i*vert+i]].row+=(space[graph[i*vert+i]].row)?0:1;
 for(i=0; i<vert; ++i) for(j=i+1; j<vert; ++j) {
  if(space[graph[i*vert+j]].row==1) return(0); else space[graph[i*vert+j]].row=2;
  if(space[graph[j*vert+i]].row==1) return(0); else space[graph[j*vert+i]].row=2;};
 for(free=space; free<space+rank; free++) if(free->row==0) return(0);

free=&space[0];
for (k=0; k<vert*vert; k++) color[k]=NULL;
for (i=0; i<vert; i++) for (j=0; j<vert; j++) {
  free->row=i; free->col=j;
  if(color[graph[i*vert+j]]==NULL) free->ptr=NULL;
  else free->ptr=color[graph[i*vert+j]];
  color[graph[i*vert+j]]=free++;}; return(1);
}

stabil (arank,vert,graph,color)
int *graph;
int vert,*arank;
struct edge **color;
{
int k,p,i,j,rank,klass,c,s,t,truth,overfl,q,oldq;
int newrank,oldnrank,oldp;
int *gamma;
struct edge *free,*w,*o,*oo;
gamma=&memory[0];
rank=*arank;
printf("%6d",rank); fflush(stdout);
do {   /*  until new colors would not appear */
  truth=0;                 /* new colors were not appear */
  newrank=rank;
  for (k=0; k<rank; k++) {  /* cycle on colors */
    overfl=0;
    klass=0;    /* number of new colors */
    p=0; /*  the begin of newgamma */
    *gamma=k;
    w=color[k];
    o=w;   /* the previous edge of color k  */
    if (w->ptr==NULL) continue;  /* new k  */

    do
    {
     triangl(graph,w->row,w->col,gamma+p,rank,vert);
       oldnrank=newrank;
       oldp=p;
       search(k,gamma,&p,&c,&klass,&s,&newrank,&truth,&q,&oldq);
       if (p>=(memlength-(vert*3+5)) || overfl==1)
	  { p=oldp; newrank=oldnrank; overfl=1;
           if(q==-1) *(gamma+oldq+4)=-1; else{ if(oldq!=-1){
           if(*(gamma+p+3)!=-1) *(gamma+*(gamma+p+3)+2)=*(gamma+p+2);
           if(*(gamma+p+2)!=-1) *(gamma+*(gamma+p+2)+3)=*(gamma+p+3);};};}
       if (oldnrank!=newrank) {
          printf("\b\b\b\b\b\b%6d",newrank); fflush(stdout);}
       if (c!=k)
        {
         o->ptr=w->ptr;
	   w->ptr=color[c];
	   color[c]=w;
	}
       else o=w;
       w=o->ptr;
    }
    while (w!=NULL); /* the last edge of color #k */
    if (overfl==1) newrank++;
    for (i=rank; i<newrank; i++)
     {
     w=color[i];
     while(w!=NULL)
       { graph[w->row*vert+w->col]=i; w=w->ptr; }
     }
     rank=newrank;
   }                 /* next color */
   if (truth==0) break;
   *arank=rank;
 }
while (1);
}

triangl(graph,i,j,newgamma,rank,vert)
int *graph;
int i,j;
int *newgamma;
{ 
struct triple *w;
int s,t,p,numval,q;
struct triple *freemem;
int *nnn; 
numval=0;                /* the number of nonzero const */
freemem=&cnst[0];
for (p=0;p<rank;p++)
    lines[p]=NULL;
for (p=0; p<vert; p++)
    {
     s=graph[i*vert+p];
     t=graph[p*vert+j];
     pack (&lines[s],&freemem,t);
    }
nnn=newgamma+5;
for (s=0; s<rank;  s++)
      {
          w=lines[s];
          while (w!=NULL)
          {
	    numval++;
	    *(nnn++)=s;
	    *(nnn++)=w->col;
	    *(nnn++)=w->val;
            w=w->ptr;
	  }
       }
*(newgamma+1)=numval;
}

pack (struct triple **line, struct triple **free, int t)
{
struct triple *w, *o;
if (*line==NULL)    /* t is the first column in the line s */
   {*line=*free;
    (*line)->col=t;
    (*line)->val=1;
    (*line)->ptr=NULL;
    (*free)++;
   }
else
   {
    if ((*line)->col>t) /* the first column in the line > t */
      {
       (*free)->col=t;
       (*free)->val=1;
       (*free)->ptr=*line;
       *line=*free;
       (*free)++;
      }
     else  /* the first column in the line <=t  */
      {
       w=*line;
       while (w!=NULL)
            {
             if(w->col==t)
               {
                w->val++;
                break;
               }
             if(w->col>t)
               {
		 (*free)->col=t;
		 (*free)->val=1;
		 (*free)->ptr=w;
		 o->ptr=*free;
		 (*free)++;
                 break;
               }
              o=w;
              w=w->ptr;
            }
             if (w==NULL)
               {
		o->ptr=*free;
		(*free)->col=t;
		(*free)->val=1;
		(*free)->ptr=NULL;
		(*free)++;
               }
      }
   }
}

search(k,gamma,ap,ac,aklass,as,anewrank,atruth,aq,aoldq)
int *gamma;
int *ac,*ap,*aklass,*as,*anewrank,*atruth,k,*aq,*aoldq;
{
int q,oldp,oldq,nexte,dl,t,i;
int c,p,klass,newrank,truth;
p=*ap;
klass=*aklass;
truth=*atruth;
newrank=*anewrank;
oldp=p;  /* the begin of newgamma */
q=0;  /* the begin of searching in gamma */
if (klass)
   {
     while (*(gamma+q+1)!=*(gamma+p+1))
     {
      if (*(gamma+q+1)>*(gamma+p+1))
	 { if(*(gamma+q+2)==-1)         /*   prev==-1  */
	     {*(gamma+q+2)=p;
	      *(gamma+p+2)=-1;
	      *(gamma+p)=newrank;
	      c=newrank;
	      newrank++;
	      klass++;
	      *(gamma+p+3)=q;
	      *(gamma+p+4)=-1;
	      truth=1;
	      p=oldp+*(gamma+oldp+1)*3+5;
	      goto IR1;
	     }
	   else
	     { if(*(gamma+p+1)>*(gamma+*(gamma+q+2)+1)) /* between */
		 {*(gamma+p+2)=*(gamma+q+2);
		  *(gamma+p+3)=q;
		  *(gamma+q+2)=p;
      /*neu*/     *(gamma+*(gamma+p+2)+3)=p;
		  *(gamma+p+4)=-1;
		  *(gamma+p)=newrank;
		  c=newrank;
		  newrank++; klass++;
		  truth=1;
		  p=p+*(gamma+p+1)*3+5;
		  goto IR1;
		 }
	       else q=*(gamma+q+2);
	     }
	 }
      else
	 { if(*(gamma+q+3)==-1)          /* next==-1  */
	     {*(gamma+q+3)=p;
	      *(gamma+p+3)=-1;
	      *(gamma+p)=newrank;
	      c=newrank;
	      newrank++;
	      klass++;
	      *(gamma+p+2)=q;
	      *(gamma+p+4)=-1;
              truth=1;
	      p=oldp+*(gamma+oldp+1)*3+5;
	      goto IR1;
	     }
	   else
	     { if (*(gamma+p+1)<*(gamma+*(gamma+q+3)+1)) /* between */
		  {*(gamma+p+3)=*(gamma+q+3);
		   *(gamma+p+2)=q;
		   *(gamma+q+3)=p;
     /*neu*/       *(gamma+*(gamma+p+3)+2)=p;
		   *(gamma+p+4)=-1;
		   *(gamma+p)=newrank;
		   c=newrank;
		   newrank++; klass++;
		   truth=1;
		   p=p+*(gamma+p+1)*3+5;
		   goto IR1;
		  }
	       else q=*(gamma+q+3);
	     }}}
     do
     { oldq=q;
       dl=*(gamma+q+1);
       nexte=*(gamma+q+4);
       q+=5;p+=5;
       for (t=1; t<dl*3; t++,q++,p++)
	   {if (*(gamma+q)!=*(gamma+p)) break;}
       if (t==dl*3)  /* old class  */
	  {
	   c=*(gamma+oldq);  /* colour */
	   p=oldp; oldq=-1;
	   break;
	  }
       if (nexte==-1)     /* create a new class */
	{  klass++;
	   *(gamma+oldp)=newrank;
	   c=newrank;
	   *(gamma+oldp+4)=-1;
	   *(gamma+oldq+4)=oldp;
	   *(gamma+oldp+2)=*(gamma+oldq+2);
	   *(gamma+oldp+3)=*(gamma+oldq+3);
	   newrank++;
	   truth=1;
	   p=oldp+*(gamma+oldp+1)*3+5;
  /*neu*/  q=nexte;
	   break;
	}
	  q=nexte;
	  p=oldp;
     }  while (q!=-1);
   }
    else { klass++;         /* the first edge of colour k  */
	   *(gamma+p+4)=-1;
	   *(gamma+p+3)=-1;
	   *(gamma+p+2)=-1;
	   *(gamma+p)=k;
	   p+=*(gamma+p+1)*3+5;
	   c=k;
	 }
IR1:*ac=c;
    *ap=p;
    *anewrank=newrank;
    *atruth=truth;
    *aklass=klass;
    *aq=q;
    *aoldq=oldq;
}


