#include "curvature.h"

static unsigned char inxt[3]  = {1,2,0};
static unsigned char iprv[3] = {2,0,1};

#define KTA     7
#define KTB     11

unsigned char idir[5]     = {0,1,2,0,1};


/* build hash table for boundary triangles */
int hashTria(pMesh mesh)
{
  pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,iadr;
  int      *hcode,*link,inival,hsize;
  unsigned char  *hvoy,i,ii,i1,i2;
  unsigned int    key;
  
  /* memory alloc */
  hcode = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(hcode);
  link  = mesh->adja;
  hsize = mesh->nt;
  hvoy  = (unsigned char*)hcode;
  
  /* init */
  inival = 2147483647;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;
  
  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if(pt->ref!=mesh->refb){
    if ( !pt->v[0] )  continue;
    
    for (i=0; i<3; i++) {
      i1 = idir[i+1];
      i2 = idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }
      
      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;
      
      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }
  }
  
  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = idir[i+1];
    i2 = idir[i+2];
    pt = &mesh->tria[k];
    
    mins = LS_MIN(pt->v[i1],pt->v[i2]);
    maxs = LS_MAX(pt->v[i1],pt->v[i2]);
    
    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = idir[ii+1];
      i2 = idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }
      
      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  
  free(hcode);
  return(1);
}




/* Return triangular ball of point ip in tria start.
 Results are stored under the form 3*kel + jel , kel = number of the tria, jel = local
 index of p within kel */
int boulept(pMesh mesh, int start, int ip, int * list){
  pTria    pt;
  pPoint   ppt;
  int           *adja,k,ilist;
  char          i,i1,i2;
  
  pt = &mesh->tria[start];
  ppt = &mesh->point[pt->v[ip]];
  ilist = 0;
  
  /* store neighbors */
  k = start;
  i = ip;
  do {
    if ( ilist > 2000)  return(-ilist);
    list[ilist] = 3*k + i;
    ++ilist;
    
    adja = &mesh->adja[3*(k-1)+1];
    i1 = inxt[i];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = inxt[i];
  }
  while ( k && k != start );
  if ( k > 0 ) return(ilist);
  
  /* check if boundary hit */
  k = start;
  i = ip;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = iprv[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;
    i  = adja[i2] % 3;
    i  = iprv[i];
    
    if ( ilist > 2000 )  return(-ilist);
    list[ilist] = 3*k + i;
    ilist++;
  }
  while ( k );
  
  return(ilist);
}











