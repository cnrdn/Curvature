#ifndef _ELASTIC_H
#define _ELASTIC_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "libmesh5.h"

#define LS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define LS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )


#define EPSD     1.e-30
#define LONMAX   8192
#define REF 2 // reference of dirichlet sphere

typedef struct {
    double   c[3],count;
    int      ref,s,mark,neword,t;
} Point;
typedef Point * pPoint;

typedef struct {
    int     v[3],ref,mark;

} Tria;
typedef Tria * pTria;

typedef struct {
    double    g[3];
    int     v[4],s,ref,mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
    int      np,na,nt,ne,hmax,hcur,ver,dim,mark,ref,npbn,refb;
    int      *adja,*tab;
    char     *name,cltyp;
    double   min[3],max[3],delta;

    pPoint    point;
    pTria    tria;
    pTetra   tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  int      dim,ver,np,nit,iter,size[2],nmat,nbcl,nps;
  double  *u,*usurf,*p,*d,err,*valp1,*grad,*gradsurf;
  char    *name,cltyp;
  
} Sol;
typedef Sol * pSol;


/* prototypes */
int     loadMesh(pMesh );
int     saveSol(pSol , pMesh);
        int hashTria(pMesh );
int     boulept(pMesh, int, int, int *);
int     mean_curvature(pMesh, double *);




#endif
