#include "curvature.h"


/* read mesh */
int loadMesh(pMesh mesh) {
    pPoint       ppt;
    pTria        pt1;
    
    float        fp1,fp2,fp3;
    int          k,inm,i;
    char        *ptr,data[256];
    
    strcpy(data,mesh->name);
    ptr = strstr(data,".mesh");
    if ( !ptr ) {
        strcat(data,".meshb");
        if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
            ptr = strstr(data,".mesh");
            *ptr = '\0';
            strcat(data,".mesh");
            if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
                fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
                return(0);
            }
        }
    }
    else if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
    }
    fprintf(stdout,"  %%%% %s OPENED\n",data);
    
    
        fprintf(stdout,"  -- READING DATA FILE %s\n",data);
    
    mesh->np = GmfStatKwd(inm,GmfVertices);
    mesh->nt = GmfStatKwd(inm,GmfTriangles);
    
    
    if ( !mesh->np ) {
        fprintf(stdout,"  ** MISSING DATA\n");
        return(0);
    }
    
    /* memory alloc */
    mesh->point = (pPoint)calloc(mesh->np+1,sizeof(Point));
    assert(mesh->point);
    if ( mesh->ne ) {
        mesh->tetra  = (pTetra)calloc(mesh->ne+1,sizeof(Tetra));
        assert(mesh->tetra);
    }
    if ( mesh->nt ) {
        mesh->tria  = (pTria)calloc(mesh->nt+1,sizeof(Tria));
        assert(mesh->tria);
    }
    
   GmfGotoKwd(inm,GmfVertices);
        for (k=1; k<=mesh->np; k++) {
            ppt = &mesh->point[k];
            if ( mesh->ver == GmfFloat ) {
                GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
                ppt->c[0] = fp1;
                ppt->c[1] = fp2;
                ppt->c[2] = fp3;
            }
            else
                GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
        }
        /* read triangles */
        GmfGotoKwd(inm,GmfTriangles);
        for (k=1; k<=mesh->nt; k++) {
            pt1 = &mesh->tria[k];
            GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
          for (i=0; i<3; i++) {
            ppt = &mesh->point[pt1->v[i]];
            if(!ppt->t) ppt->t = k;
          }
}
            fprintf(stdout,"  %%%% NUMBER OF VERTICES   %8d\n",mesh->np);
            if ( mesh->nt )  fprintf(stdout,"  %%%% NUMBER OF TRIANGLES  %8d\n",mesh->nt);
        
    
    GmfCloseMesh(inm);
    return(1);
}


int saveSol(pSol sol, pMesh mesh) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,ia,i,inm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128];
  
  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  sprintf(data,"%s.sol",data);
  
  sol->ver = 2;
  sol->dim = mesh->dim;
  sol->np = mesh->np;
  
  if ( !(inm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  
  type = 1;
  typtab[0] = GmfSca;
  
  /* write sol */
  GmfSetKwd(inm,GmfSolAtVertices,sol->np,type,typtab);
  if ( sol->ver == GmfFloat ) {
    for (k=0; k<sol->np; k++) {
      ia = sol->dim*k;
      for (i=0; i<sol->dim; i++)
        fbuf[i] = sol->u[ia+i];
      GmfSetLin(inm,GmfSolAtVertices,fbuf);
    }
  }
  else {
    for (k=0; k<sol->np; k++) {
      ia = sol->dim*k;
      for (i=0; i<sol->dim; i++)
        dbuf[i] = sol->u[ia+i];
      GmfSetLin(inm,GmfSolAtVertices,dbuf);
    }
  }
  
  GmfCloseMesh(inm);
  return(1);
}


