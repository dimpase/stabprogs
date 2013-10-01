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

#include <stdio.h>
#include <stdlib.h>

#define vm 150
#define vm2 vm*vm
#define vm3 vm*vm2

#define SCAN(r, s) while( r ) { xp=bucket + s ;\
		    if (xp->last==0) xp->first=p;\
		    else xp->last->next=p; xp->last=p;\
		    p=p->next;}

#define SORT(r, s) xp=bucket; while(xp->last==0) xp++; p=xp->first; yp=xp;\
		   while(xp<bucket+ r ) {\
		    if(xp->last) { yp->last->next=xp->first;\
		     q=xp->first;\
		     while(q!=xp->last) { j=q->c; q->c=col;\
		      col+=(j==q->next->c)?0:1; q=q->next;}\
		     q->c=col; s ; yp=xp;}; xp++;}

#define PICK if(pm2->next) {\
	      if(pm2->prev) { pm2->prev->next=pm2->next;\
	       pm2->next->prev=pm2->prev;\
	       if((ewc+pm2->a)->point==pm2) (ewc+pm2->a)->point=pm2->next;}\
	      else { pm1=pm1->next; pm1->prev=0;\
	       (ewc+pm2->a)->point=pm2->next;};}\
	     else { if(pm2->prev) pm2->prev->next=0; else pm1=0;}

#define COLOR(r, s) q=xp->first; while(q!=xp->last) {\
		     j=q-> r ; q-> r =col; col+=(j==q->next-> r )?0:1;\
		     q=q->next;}; q-> r = s ;

#define FARBE(r) while(q!= r ) { q->a=col; ewc[col].m++;\
		  if(q->c!=q->next->c) { col++; ewc[col].point=q->next;\
		  ewc[col].m=0;}; q=q->next;};
FILE *f;
struct triangle{ struct multiset *ij, *jk, *ik;
 struct triangle *next; int l, c;};
struct multiset{ int c, a, l;
 struct triangle *joey; struct multiset *next, *prev;};
struct bucket_triangle{ struct triangle *first, *last;};
struct bucket_multiset{ struct multiset *first, *last;};
struct color_class{ int m;
	struct multiset *point;}ewc[vm2];
int n, maxcolor, newcolor;
int max[vm];



struct triangle *
sort_triangles(struct triangle *p) {
 struct triangle *q;
 static struct bucket_triangle bucket[vm],*xp, *yp;
 int i, j, col, *pc;
 for(i=0; i<4; i++) { for(xp=bucket; xp<bucket+n; xp++) xp->last=0;
  switch(i) {
   case 0: SCAN(p, (p->jk->a)/n);
   case 1: SCAN(p, (p->jk->a)%n);
   case 2: SCAN(p, (p->ij->a)/n);
   case 3: SCAN(p, (p->ij->a)%n);};
  col=0; SORT(n, col++); yp->last->next=0;};
 q=p;
 while(q) { q->l=q->ik->l++; q=q->next;};
 for(xp=bucket; xp<bucket+n; xp++) xp->last=0;
 SCAN(p, p->l); pc=max; col=0;
 SORT(n, *(pc++)=col;col=0); yp->last->next=0; return(p);
}

struct multiset *
sort_multisets(struct multiset *p, struct multiset *px) {
 static struct bucket_multiset bucket[vm2],*xp, *yp;
 int *jp, col, j, i, may, x, hi;
 struct multiset *zp, *q, *mp;
 struct triangle *pt;

 zp=p; while(zp) { zp->c=0; zp=zp->next;};
 for(jp=max+n-1; jp>=max; jp--) {
  for(xp=bucket; xp<=bucket+*jp; xp++) xp->last=0;
  SCAN(p && p->l>(jp-max), p->joey->c; p->joey->c=0);
  zp=p; col=1; yp=bucket; p=yp->first;
  for(xp=bucket; xp<bucket+*jp; xp++) { yp++; xp->last->next=yp->first;
   COLOR(c, col++);}; xp=yp; COLOR(c, col);
  yp->last->next=zp; zp=p;
  while(zp && zp->l>(jp-max)) { pt=zp->joey;
   zp->joey=zp->joey->next; free(pt); zp=zp->next;};};

 newcolor=maxcolor; zp=p;
 while(zp) {zp->joey=0; zp->l=0; zp=zp->next;};
 for(i=0; i<2; i++) {
  for(xp=bucket; xp<bucket+n; xp++) xp->last=0;
  switch(i) {
   case 0: SCAN(p, (p->a)%n);
   case 1: SCAN(p, (p->a)/n);};
  col=0; xp=bucket; while(xp->last==0) xp++; p=xp->first; yp=xp;
  while(xp<bucket+n) {
   if(xp->last) { yp->last->next=xp->first;
    COLOR(l, col++); yp=xp;}; xp++;}; yp->last->next=0;};
 hi=col;
 for(xp=bucket; xp<bucket+hi; xp++) xp->last=0;
 SCAN(p, p->l; p->l=0); p=bucket->first; yp=bucket;
 for(xp=bucket; xp<bucket+hi-1; xp++) { yp++; xp->last->next=yp->first;};
 yp->last->next=px; zp=bucket->first; zp->prev=0;
 while(zp->next!=px) { zp->next->prev=zp; zp=zp->next;};
 if(px) px->prev=zp;

 for(xp=bucket; xp<bucket+hi; xp++) {
  q=xp->first; zp=q; mp=q; i=1; j=1; may=1;
  while(q!=xp->last) {
   if(q->c<q->next->c) { if(j>may) { may=j; mp=zp;};
    zp=q->next; j=0;}; i++; j++; q=q->next;};
  if(j>may) { may=j; mp=zp;}; x=ewc[xp->first->a].m-i;
  if(x<may) { ewc[xp->first->a].m = may;
   if(x) { newcolor++; zp=ewc[xp->first->a].point;
    ewc[newcolor].point=zp; ewc[newcolor].m=x;
    for(j=0; j<x; j++) { zp->a=newcolor; zp=zp->next;};};
   ewc[xp->first->a].point=mp; col=newcolor+1;
   q=xp->first; ewc[col].point=q; ewc[col].m=0; FARBE(mp);
   for(j=1; j<may; j++) q=q->next;
   if(q!=xp->last) {q=q->next; ewc[col].point=q; FARBE(xp->last);
    ewc[col].m++; newcolor=col; q->a=col++;};}
  else { col=newcolor+1; q=xp->first; ewc[col].point=q; ewc[col].m=0;
   FARBE(xp->last); ewc[col].m++; newcolor=col; q->a=col++;};
  newcolor=col-1;}; return(p);
}

struct multiset *
multisets(struct triangle *pt1, struct multiset *pm1) {
 struct multiset *p, *pm2;
 struct triangle *pt2;
 static struct bucket_multiset bucket[vm], *xp, *yp;
 int i, j;
 p=pt1->ik; pt2=pt1->next; pm2=pt1->ik; PICK;
 pm2->next=0; p=pm2; pm2->joey=pt1; pt1=pt2;
 while(pt1) { pt2=pt1->next; pm2=pt1->ik;
  if(pm2->joey) pt1->next=pm2->joey;
  else { PICK; pm2->next=p; p=pm2;};
  pm2->joey=pt1; pt1=pt2;};
 for(xp=bucket; xp<bucket+n; xp++) xp->last=0;
 SCAN(p, (p->l)-1); xp=bucket+n-1; while(xp->last==0) xp--;
 p=xp->first; yp=xp;
 while(xp>=bucket) {
  if(xp->last) { yp->last->next=xp->first; yp=xp;}; xp--;};
 yp->last->next=0; return(sort_multisets(p, pm1));
}

struct multiset *
presort(struct multiset *p) {
 static struct bucket_multiset bucket[vm], *xp, *yp;
 struct color_class *pc;
 for(xp=bucket; xp<bucket+newcolor+1; xp++) xp->last=0;
 SCAN(p, p->a; ewc[p->a].m++;);
 p=bucket->first; pc=ewc; yp=bucket;
 for(xp=bucket; xp<bucket+newcolor; xp++) {
  yp++; (xp->last)->next=yp->first; (pc++)->point=xp->first;};
 (yp->last)->next=0; pc->point=yp->first;
 return(p);
}

struct triangle *
triangles(struct multiset S[vm][vm]) {
 struct triangle help, *pt; struct color_class *pc; struct multiset *pm;
 int i, j, x, y, v;
 pt=&help;
 for(pc=&ewc[maxcolor+1]; pc<=&ewc[newcolor]; pc++) {
  v=(pc-ewc); pm=pc->point;
   for(i=0; i<pc->m; i++) { x=(pm-&S[0][0])/vm; y=(pm-&S[0][0])%vm;
    for(j=0; j<n; j++) { if((pt->next=(struct triangle *)
     malloc(sizeof(struct triangle)))==0) return(0); pt=pt->next;
     pt->ik=pm; pt->ij=&S[x][j]; pt->jk=&S[j][y]; pt->c=1;};
    for(j=0; j<n; j++) { 
     if(S[j][y].a<=maxcolor) {
      if(S[j][x].a<v || (S[j][x].a==v && (&S[j][x]-pm)>0)) {
       if((pt->next=(struct triangle *)
        malloc(sizeof(struct triangle)))==0) return(0); pt=pt->next;
	pt->ik=&S[j][y]; pt->ij=&S[j][x]; pt->jk=pm; pt->c=1;};};};
    for(j=0; j<n; j++) {
     if(S[x][j].a<=maxcolor) {
      if(S[y][j].a<v || (S[y][j].a==v && (&S[y][j]-pm)>0)) {
       if((pt->next=(struct triangle *)
        malloc(sizeof(struct triangle)))==0) return(0); pt=pt->next;
	pt->ik=&S[x][j]; pt->ij=pm; pt->jk=&S[y][j]; pt->c=1;};};};
    pm=pm->next;};};
 pt->next=0; return(help.next);
}

void main(int narg, char *arg[10]){
 long cpu;
 int i, j/*, l, u*/;
 struct color_class *pc;
 struct multiset *pm1, *pm2, S[vm][vm];
 struct triangle *pt;

 f=fopen(arg[1],"r");
 if (f==NULL) { printf("\n file does not exist\n"); exit(0);}
 fscanf (f,"%d%d",&newcolor,&n); newcolor--;
 for(i=0; i<n; ++i) for(j=0; j<n; ++j) { fscanf(f,"%d",&S[i][j].a);
  if(S[i][j].a<0 || S[i][j].a>newcolor) 
   { printf("please check your input!\n"); return;};
  S[i][j].joey=0; S[i][j].l=0; S[i][j].c=0;};
 for(pc=ewc; pc<ewc+n*n; pc++) pc->m=0; 
 for(i=0; i<n; ++i) ewc[S[i][i].a].m+=(ewc[S[i][i].a].m)?0:1;
 for(i=0; i<n; ++i) for(j=i+1; j<n; ++j) {
   if(ewc[S[i][j].a].m==1) { printf("please check your input!\n"); return;}
   else ewc[S[i][j].a].m=2;
   if(ewc[S[j][i].a].m==1) { printf("please check your input!\n"); return;}
   else ewc[S[j][i].a].m=2;};
 for(pc=ewc; pc<=ewc+newcolor; pc++) {
  if(pc->m==0) { printf("please check your input!\n"); return;}
  else pc->m=0;};

 pm1=0;
 for(i=n-1; i>=0; i--) for(j=n-1; j>=0; j--){S[i][j].next=pm1; pm1=&S[i][j];};
 pm1=presort(&S[0][0]);
 pm2=pm1; pm1->prev=0;
 while(pm1->next) {(pm1->next)->prev=pm1; pm1=pm1->next;};
 maxcolor=0; printf("\nnumber of colors: %6d",newcolor+1);
 cpu=clock()/1000; i=0;
 while(newcolor>maxcolor) { i++;
  if((pt=triangles(S))==0) { printf("\nspaceproblems!\n"); newcolor=maxcolor;}
  else{ maxcolor=newcolor;
   pt=sort_triangles(pt);
   pm1=multisets(pt, pm2); pm2=pm1; printf("\b\b\b\b\b\b%6d",newcolor+1);
   fflush(stdout);};};
 cpu=clock()/1000-cpu;
 for(pc=ewc; pc<=ewc+maxcolor; pc++) pc->m=-1; j=0;
 for(i=0; i<n; ++i) { if(ewc[S[i][i].a].m<0) ewc[S[i][i].a].m=j++;};
 printf("\n\nnumber of cells:  %6d",j);
 for(i=0; i<=maxcolor; ++i) if(ewc[i].m<0) ewc[i].m=j++;
/* printf("\n%ld msec",cpu); */
 printf("\n\nadjacency matrix of the cellular algebra:\n\n");	
 for(i=0; i<n; ++i) { for(j=0; j<n; ++j) { S[i][j].a=ewc[S[i][j].a].m;
  if(S[i][j].a <10) printf("  %d ",S[i][j].a);
  else {if(S[i][j].a <100) printf(" %d ",S[i][j].a);
  else printf("%d ",S[i][j].a);};} printf("\n");};
}

