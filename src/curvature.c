#include "curvature.h"

static const unsigned char inxt2[3] = {1,2,0}; /*!< next vertex of triangle: {1,2,0} */
static const unsigned char iprv2[3] = {2,0,1}; /*!< previous vertex of triangle: {2,0,1} */

static int parsar(int argc,char *argv[],pMesh mesh) {
  int     i;
  
  i = 1;
  while ( i < argc ) {
    
    if ( mesh->name == NULL ) {
      mesh->name = (char*) calloc(strlen(argv[i])+1,sizeof(char));
      strcpy(mesh->name,argv[i]);
    }
    
    i++;
  }
  
  /* check file names */
  if ( mesh->name == NULL ) {
    mesh->name = (char *)calloc(128,sizeof(char));
    assert(mesh->name);
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",mesh->name);
  }
  return(1);
}


/* check if the triangle p0p1p2 is obtuse */
/* return -1 if the triangle is not obtuse
 0,1,2 if the triangle is obtuse at 0,1,2 */
static int check_obtuse_tria(pPoint p0, pPoint p1, pPoint p2) {
  double dot0, dot1, dot2;
  int i;
  
  dot0= 0.0;
  dot1 = 0.0;
  dot2 = 0.0;
  
  for(i=0; i<3; i++){
    dot0 += (p1->c[i] - p0->c[i])*(p2->c[i] - p0->c[i]);
    dot1 += (p0->c[i] - p1->c[i])*(p2->c[i] - p1->c[i]);
    dot2 += (p0->c[i] - p2->c[i])*(p1->c[i] - p2->c[i]);
  }
  
  if(dot0*dot1*dot2>=0) return (-1);
  
  if(dot0<0) return (0);
  if(dot1<0) return (1);
  if(dot2<0) return (2);
  
  return(-2);
}

/* return cotangent of the angle p0p1p2 */
static double cotan(pPoint p0, pPoint p1, pPoint p2) {
  double dot, cross, ux, uy, uz, vx, vy, vz, p[3],cot;
  
  
  ux = p0->c[0] - p1->c[0];
  uy = p0->c[1] - p1->c[1];
  uz = p0->c[2] - p1->c[2];
  
  vx = p2->c[0] - p1->c[0];
  vy = p2->c[1] - p1->c[1];
  vz = p2->c[2] - p1->c[2];
  
  p[0] = uy*vz - uz*vy;
  p[1] = uz*vx - ux*vz;
  p[2] = ux*vy - uy*vx;
  
  dot = ux*vx + uy*vy+uz*vz;
  cross = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  cot = dot/cross;
  
  return cot;
}

/* compute triangle area and unit normal in 3d */
static double area_3d(double *a, double *b, double *c, double *n) {
  double    ux, uy, uz, vx, vy, vz, dd, dd1;
  
  ux = b[0] - a[0];
  uy = b[1] - a[1];
  uz = b[2] - a[2];
  
  vx = c[0] - a[0];
  vy = c[1] - a[1];
  vz = c[2] - a[2];
  
  n[0] = uy*vz - uz*vy;
  n[1] = uz*vx - ux*vz;
  n[2] = ux*vy - uy*vx;
  dd   = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if ( dd > EPSD ) {
    dd1 = 1.0 / dd;
    n[0] *= dd1;
    n[1] *= dd1;
    n[2] *= dd1;
  }
  
  return 0.5*dd;
}


/**************************************************************************/
/*compute the mean curvature at mesh vertices using the cotangent formula */
/**************************************************************************/

int mean_curvature(pMesh mesh, double *curv){
  pPoint ppt,pptj,pptjj,pindex;
  pTria pt;
  int k,i,ilist=0,ip, list[2000],l1,l2,m,flag;
  double cotalpha, cotbeta, curvpt[3],area,n[3];
  
  memset(curv,0,mesh->np*sizeof(double));
  
  for (k = 1; k <= mesh->np; k++) {
    
    ppt = &mesh->point[k];
    pt=&mesh->tria[ppt->t]; // find one triangle sharing the vertex
    /* check the position of ppt in pt */
    for (m=0; m<3; m++){
      pindex = &mesh->point[pt->v[m]];
      if(pindex->c[0]==ppt->c[0]) ip = m;
    }
    assert(ip<3);
    
    /* get the 1 ring of vertex ppt */
    area = 0.0;
    curvpt[0] = 0.0;
    curvpt[1] = 0.0;
    curvpt[2] = 0.0;
    
    ilist = boulept(mesh,ppt->t,ip,list);
    assert(ilist>0);
    
    for (m=0; m<ilist; m++){
      
      pt=&mesh->tria[list[m]/3];
      ip = list[m]%3;
      if(k==pt->v[ip]) continue;
      l1 = inxt2[ip];
      pptj = &mesh->point[pt->v[l1]];
      l2 = iprv2[ip];
      pptjj = &mesh->point[pt->v[l2]];
      assert(ip!=l1);
      assert(l1!=l2);
      assert(l2!=ip);
      
      flag = check_obtuse_tria(ppt,pptj,pptj);
      
      cotalpha = cotan(ppt,pptj,pptjj);
      cotbeta = cotan(ppt,pptjj,pptj);
      
      for(i=0;i<3;i++) {
        curvpt[i] += cotalpha*(ppt->c[i] - pptjj->c[i]);
        curvpt[i] += cotbeta*(ppt->c[i] - pptj->c[i]);
      }
      
      /*compute Voronoi regions */
      if(flag==-1) {
        for(i=0;i<3;i++){
          area += cotalpha*(ppt->c[i] - pptjj->c[i])*(ppt->c[i] - pptjj->c[i])/8;
          area += cotbeta*(ppt->c[i] - pptj->c[i])*(ppt->c[i] - pptj->c[i])/8;
        }
      }
      else if(flag==ip)
        area+= 0.5*area_3d(ppt->c, pptjj->c, pptj->c, n);
      else area += 0.25*area_3d(ppt->c, pptjj->c, pptj->c, n);
      
      
      for(i=0;i<3;i++) curvpt[i]*= 1./(2*area);
      curv[k-1] = curvpt[0]*curvpt[0] + curvpt[1]*curvpt[1] + curvpt[2]*curvpt[2];
      curv[k-1] = 0.5*sqrt(curv[k]);
      
      
    }
  }
  
  return(1);
}


int main(int argc, char**argv){
  
  Mesh     mesh;
  Sol      sol;
  
  /* default values */
  memset(&mesh,0,sizeof(Mesh));
  
  /* command line */
  if ( !parsar(argc, argv, &mesh) )  return 1;
  
  /* load data */
  if ( !loadMesh(&mesh) ) return 1;
  
  sol.u = (double*)calloc(mesh.np, sizeof(double));
  assert(sol.u);
  mesh.adja = (int*)calloc(3*mesh.nt+5,sizeof(int));
  assert(mesh.adja);
  
  if ( !hashTria(&mesh)) return 0;
  if ( !mean_curvature(&mesh,sol.u)) return 0;
  if(! saveSol(&sol,&mesh) ) return 1;
  
  /*free memory */
  free(sol.u);
  
  return 0;
}
