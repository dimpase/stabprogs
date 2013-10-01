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

#include  <stdio.h>
#include  <stdlib.h>
#define maxvert 200
#define maxrank maxvert*maxvert
#define memlength maxrank*5
#define mlong int  /* to fix an MSDOS artefact */
FILE *f;
struct triple
{int col; int val; struct triple *ptr;};
struct edge
{ int row; int col; struct edge *ptr;};

void main(int narg, char *arg[10])
{
struct edge *color[maxrank];
int i,j,rank,vert;
int graph[maxvert][maxvert];
long time; float t;
long start_time,end_time;
f=fopen(arg[1],"r");        /* char */
if (f == NULL) { printf("\n file does not exist\n"); exit(0);}
fscanf (f,"%d%d",&rank,&vert);
for (i=0;i<vert;i++) for(j=0;j<vert;j++) fscanf (f, "%d", &graph[i][j]);
i=edgepack (graph,rank,vert,color); if(i==0) return;

start_time=clock()/1000;
printf ("\n\n number of colors: ");
stabil (&rank,vert,graph,color);
end_time=clock()/1000-start_time;
for(i=0; i<rank; i++) color[i]->row=-1; j=0;
for(i=0; i<vert; ++i) {
 if(color[graph[i][i]]->row<0) color[graph[i][i]]->row=j++;};
printf ("\n\n number of cells:  %6d", j);
for(i=0; i<rank; ++i) if(color[i]->row<0) color[i]->row=j++;
printf("\n\n adjacency matrix of the cellular algebra:\n\n");
 for(i=0; i<vert; ++i) { for(j=0; j<vert; ++j) { 
  graph[i][j]=color[graph[i][j]]->row;
  if(graph[i][j] <10) printf("  %d ",graph[i][j]);
  else {if(graph[i][j] <100) printf(" %d ",graph[i][j]);
  else printf("%d ",graph[i][j]);};} printf("\n");};
/* printf("\n\n%ld msec \n\n",end_time);*/
}

int edgepack (graph,rank,vert,color)
int graph[maxvert][maxvert];
int vert,rank;
struct edge *color[maxrank];
{
int k,i,j;
static struct edge space[maxrank];
struct edge *free;

 for(i=0; i<vert; ++i) for(j=0; j<vert; ++j) {
  if(graph[i][j]<0 || graph[i][j]>rank-1) 
   { printf("please check your input!\n"); vert=0; return(0);};};
 for(free=space; free<space+rank; free++) free->row=0; 
 for(i=0; i<vert; ++i) space[graph[i][i]].row+=(space[graph[i][i]].row)?0:1;
 for(i=0; i<vert; ++i) for(j=i+1; j<vert; ++j) {
  if(space[graph[i][j]].row==1) {printf("please check your input!\n");
   vert=0; return(0);} else space[graph[i][j]].row=2;
  if(space[graph[j][i]].row==1) {printf("please check your input!\n");
   vert=0; return(0);} else space[graph[j][i]].row=2;};
 for(free=space; free<space+rank; free++) {
  if(free->row==0) {printf("please check your input!\n");vert=0;return(0);};};

free=&space[0];
for (k=0; k<maxrank; k++)
    color[k]=NULL;
for (i=0; i<vert; i++)
   for (j=0; j<vert; j++)
       {
        free->row=i;
        free->col=j;
        if (color[graph[i][j]]==NULL)
            free->ptr=NULL;
	else  free->ptr=color[graph[i][j]];
        color[graph[i][j]]=free;
        free++;
       }
return(1);
}

stabil (arank,vert,graph,color)
int graph[maxvert][maxvert];
int vert,*arank;
struct edge *color[maxrank];
{
int k,p,i,j,rank,klass,c,s,t,truth,overfl;
/* int actcol[maxrank];  */
int newrank,oldnrank,oldp;
int *gamma;
int memory[memlength];
struct edge *free,*w,*o,*oo;
gamma=&memory[0];
rank=*arank;
printf("%6d",rank); fflush(stdout);
/* for (k=0; k<maxrank; k++)
      actcol[k]=1;  */
do    /*  until new colors would not appear */
 {
  truth=0;                 /* new colors were not appear */
  newrank=rank;
  for (k=0; k<rank; k++)   /* cycle on colors */
   {
    overfl=0;
    klass=0;    /* number of new colors */
/*    edgeprint (rank,color);  */
    p=0; /*  the begin of newgamma */
    *(gamma)=k;
    w=color[k];
    o=w;   /* the previous edge of color k  */
    if (w->ptr==NULL) continue;  /* new k  */
    do
    {
     triangl(graph,w->row,w->col,gamma+p,rank,vert /*,actcol */);
/*     if (p==0)   */    /* the first edge of colour k */
/*      {
      for (i=5; i<*(gamma+1)*3+5; i+=3)
	 {
	  if (actcol[*(gamma+i)]%2==1) break;
	  if (actcol[*(gamma+i+1)]%2==1) break;
	 }
      if (i==*(gamma+1)*3+5) continue; */   /*  new k  */
/*      } */
       oldnrank=newrank;
       oldp=p;
       search(k,gamma,&p,&c,&klass,&s,&newrank,&truth);
       if (p>=(memlength-(vert*3+5)) || overfl==1)
	  { p=oldp; newrank=oldnrank; overfl=1;}
       if (oldnrank!=newrank) {
        printf("\b\b\b\b\b\b%6d",newrank); fflush(stdout);}
       if (c!=k)
        {
         o->ptr=w->ptr;
	/* if (color[c]==NULL)
	  {actcol[c]+=2;
	   actcol[k]+=2;}   */
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
       { graph[w->row][w->col]=i; w=w->ptr; }
     }
     rank=newrank;
   }                 /* next color */
   if (truth==0) break;
   *arank=rank;
/*   for (i=0; i<rank; i++)
      { if (actcol[i]%2==0)
	actcol[i]=1;
	else actcol[i]=0;} */
 }
while (1);
}

triangl(graph,i,j,newgamma,rank,vert /*,actcol */)
/* int actcol[maxrank];  */
int graph[maxvert][maxvert];
int i,j;
mlong *newgamma;
{ static struct triple *lines[maxrank];
struct triple *w;
int s,t,p,numval,q;
 struct triple cnst[maxvert],*freemem;
mlong *nnn; 
numval=0;                /* the number of nonzero const */
freemem=&cnst[0];
for (p=0;p<rank;p++)
    lines[p]=NULL;
for (p=0; p<vert; p++)
    {
     s=graph[i][p];
     t=graph[p][j];
  /*   if (actcol[s]%2==0 && actcol[t]%2==0)  continue;  */
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

search(k,gamma,ap,ac,aklass,as,anewrank,atruth)
mlong *gamma;
int *ac,*ap,*aklass,*as,*anewrank,*atruth,k;
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
		   *(gamma+p+4)=-1;
		   *(gamma+p)=newrank;
		   c=newrank;
		   newrank++; klass++;
		   truth=1;
		   p=p+*(gamma+p+1)*3+5;
		   goto IR1;
		  }
	       else q=*(gamma+q+3);
	     }
	 }
     }
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
	   p=oldp;
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
}


