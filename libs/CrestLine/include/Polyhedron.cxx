/*
Fast and Robust Detection of Crest Lineson Meshes C++ code
Copyright:(c) Shin Yoshizawa, 2004
E-mail: shin.yoshizawa@mpi-sb.mpg.de
URL: http://www.mpi-sb.mpg.de/~shin
Affiliation: Max-Planck-Institut fuer Informatik: Computer Graphics Group 
 Stuhlsatzenhausweg 85, 66123 Saarbruecken, Germany
 Phone +49 681 9325-408 Fax +49 681 9325-499 

 All right is reserved by Shin Yoshizawa.
This C++ sources are allowed for only primary user of 
research and educational purposes. Don't use secondary: copy, distribution, 
diversion, business purpose, and etc.. 
 */
#include "time.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
//#include<unistd.h>
#include"Point3d.h"
#include"PointTool.h"
#include"IDList.h"
#include"PolarList.h"
#include"IDSet.h"
#include"Polyhedron.h"
#include"SvdSolve.h"
#include"Eigens.h"
#include "Edge.h"
#define PI 3.1415926535897932385
#define getk1t1(a00,a10) 6.0*(a00*(a00*a00*ab[4]+a10*a10*ab[6])+a10*(a00*a00*ab[5]+a10*a10*ab[7]));

void Polyhedron::memoryallocate(int dV,int dF){
  
  if(numberV!=0){
      memorydelete();
    }
  
  
  numberV=dV;
  numberF=dF;
  
  point = new Point3d* [numberV];
  
  Face = new int* [numberF];
  label = new int[numberV];
  IDtool = new IDSet();
  IHead = new IDList* [numberV];
  ITail = new IDList* [numberV];
  INHead = new IDList* [numberV];
  INTail = new IDList* [numberV];
  
  
  VHead = new IDList* [numberV];
  VTail = new IDList* [numberV];
  FHead= new IDList* [numberV];
  FTail= new IDList* [numberV];
  
  neighborF= new int[numberV];
  neighborIN= new int[numberV];
  
  
  boundary = new int[numberV];
  neighborI = new int[numberV]; 
  
  k1 = new double[numberV];
  k2 = new double[numberV];
  ks1 = new double[numberV];
  ks2 = new double[numberV];
  
  
  t1 = new Point3d* [numberV];
  t2 = new Point3d* [numberV];
  normal  = new Point3d* [numberV];

  centroid  = new Point3d* [numberV];
  numboundary=0;
  bc = new Point3d* [10];
  
  dbc = new Point3d* [10];
  int i;
  
  
for(i=0;i<10;i++){
    bc[i] = new Point3d(0.0,0.0,0.0);
    dbc[i] = new Point3d(0.0,0.0,0.0);
    
  }

  
  for(i=0;i<numberV;i++){
      boundary[i]=0;
      neighborIN[i]=0;
      neighborI[i]=0;
      neighborF[i]=0;
      
      ks1[i]=0.0;
      ks2[i]=0.0;
      k1[i]=0.0;
      k2[i]=0.0;
      
    
    }
    

}
void Polyhedron::setPoint(int i,double dx,double dy,double dz){
  point[i] = new Point3d(dx,dy,dz);
  
  normal[i] = new Point3d(0.0,0.0,0.0);
  centroid[i] = new Point3d(0.0,0.0,0.0);

  t1[i] = new Point3d(0.0,0.0,0.0);
  t2[i] = new Point3d(0.0,0.0,0.0);
  k1[i] = 0.0;
  k2[i] = 0.0;

  IHead[i] = new IDList();
  ITail[i] = new IDList();
  IHead[i]->next = ITail[i];
  ITail[i]->back = IHead[i];
  INHead[i] = new IDList();
  INTail[i] = new IDList();
  INHead[i]->next = INTail[i];
  INTail[i]->back = INHead[i];
  
  
  FHead[i] = new IDList();
  FTail[i] = new IDList();
  FHead[i]->next = FTail[i];
  FTail[i]->back = FHead[i];
  
  VHead[i] = new IDList();
  VTail[i] = new IDList();
  VHead[i]->next = VTail[i];
  VTail[i]->back = VHead[i];
  
}
void Polyhedron::SetBoundaryLines(){
    int i=0;numboundary=0;
    if(neighborI!=NULL && neighborF!=NULL && boundary!=NULL)
      for(i=0;i<numberV;i++){
	if(((neighborI[i]) == neighborF[i]) && 
	   (neighborF[i] !=0) && 
	   (neighborI[i] !=0)){
	  boundary[i] = 0;
	}else{
	  boundary[i] = 1;
	  numboundary++;
	}
      }
    printf("number of boundary points = %d\n",numboundary);
    //System.out.println("number of boundary points = "+numboundary);
}



void Polyhedron::setFace(int i,int di,int dj,int dk){
  Face[i] = new int[3];
  Face[i][0] = di;
  Face[i][1] = dj;
  Face[i][2] = dk;
  

  /* One */
  neighborF[di]++;
  
  IDtool->AppendISort(dj,IHead[di],ITail[di],di,neighborI);
  IDtool->AppendISort(dk,IHead[di],ITail[di],di,neighborI);
  
  IDtool->AppendVF(dj,VTail[di]);
  IDtool->AppendVF(dk,VTail[di]);
  /* Two */
  neighborF[dj]++;
  
  IDtool->AppendISort(di,IHead[dj],ITail[dj],dj,neighborI);
  IDtool->AppendISort(dk,IHead[dj],ITail[dj],dj,neighborI);
  
  IDtool->AppendVF(dk,VTail[dj]);
  IDtool->AppendVF(di,VTail[dj]);
  /* Three */
  neighborF[dk]++;
  
  IDtool->AppendISort(di,IHead[dk],ITail[dk],dk,neighborI);
  IDtool->AppendISort(dj,IHead[dk],ITail[dk],dk,neighborI);
  
  IDtool->AppendVF(di,VTail[dk]);
  IDtool->AppendVF(dj,VTail[dk]);
}



/*
void Polyhedron::readmesh(char *filename,char *fileout){
  FILE *in=NULL;
  in = fopen(filename,"r");
  
  int i,j;
  int di=0;
  int dj=0;
  int dk=0;
  double dx=0.0;
  double dy=0.0;
  double dz=0.0;
  int dV=0;
  int dF =0;
  fscanf(in,"%d",&dV);
  fscanf(in,"%d",&dF);
  fscanf(in,"%d",&MAXNEIGHBORLAVEL);
  fscanf(in,"%d",&ridge);
  MAXNEIGHBORLAVEL++;
  
  memoryallocate(dV,dF);
  //printf("End memoryallocate\n");
  
  for(i=0;i<numberV;i++){
    fscanf(in,"%lf %lf %lf",&dx,&dy,&dz);
    //printf("%lf %lf %lf\n",dx,dy,dz);
    setPoint(i,dx,dy,dz);
  }
  for(i=0;i<numberF;i++){
    fscanf(in,"%d %d %d",&di,&dj,&dk);
    //printf("%d %d %d\n",di,dj,dk);
    //printf("%d %d %d\n",di,dj,dk);
    setFace(i,di,dj,dk); 
    //printf("setFace\n");
    IDtool->AppendVF(i,FTail[Face[i][0]]);
    IDtool->AppendVF(i,FTail[Face[i][1]]);
    IDtool->AppendVF(i,FTail[Face[i][2]]);
    
  }
  
  //printf("Read\n");
  fclose(in);
  
  / * feature analysis * /
  //printf("Read\n");
  SetBoundaryLines();
  //printf("SetBoundaryLines\n");

  setProperty1();
  //printf("setProperty1\n");

    SaveRmesh(fileout);
    //printf("SaveRmesh\n");
   
   
}
*/

void Polyhedron::readmesh(int** vadj,int* vadjNum,double* vertices, unsigned vNum, unsigned* faces, unsigned fNum, unsigned neighboreSize, unsigned ridgeTension,
	unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval)
{
	int i,j;
	int di=0;
	int dj=0;
	int dk=0;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	int dV=vNum;
	int dF =fNum;
	MAXNEIGHBORLAVEL=neighboreSize;
	MAXNEIGHBORLAVEL++;
	ridge=ridgeTension;

	memoryallocate(dV,dF);
	//printf("End memoryallocate\n");

	for(i=0;i<numberV;i++){
		dx=vertices[i*3+0];
		dy=vertices[i*3+1];
		dz=vertices[i*3+2];
		setPoint(i,dx,dy,dz);
	}
	for(i=0;i<numberF;i++){
		di=faces[i*3+0];
		dj=faces[i*3+1];
		dk=faces[i*3+2];

		setFace(i,di,dj,dk); 

		IDtool->AppendVF(i,FTail[Face[i][0]]);
		IDtool->AppendVF(i,FTail[Face[i][1]]);
		IDtool->AppendVF(i,FTail[Face[i][2]]);

	}


	/* feature analysis */
	//printf("Read\n");
	SetBoundaryLines();
	//printf("SetBoundaryLines\n");

	setProperty1(vadj,vadjNum,pNum,cNum,eNum,fPoints,fPointType,fFeature,fEdges,interval);
	//printf("setProperty1\n");

	//SaveRmesh(fileout);
	//printf("SaveRmesh\n");


}

void Polyhedron::getCurvature(int*& vids,unsigned vidNum,int**& vadj,int*& vadjNum,double*& vertices, unsigned vNum, unsigned*& faces, unsigned fNum, unsigned neighboreSize,
	double*& mg1,double*& mg2,double*& dr1,double*& dr2)
{
	int i,j;
	int di=0;
	int dj=0;
	int dk=0;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	int dV=vNum;
	int dF =fNum;
	MAXNEIGHBORLAVEL=neighboreSize;
	MAXNEIGHBORLAVEL++;
	//ridge=ridgeTension;

	memoryallocate(dV,dF);
	//printf("End memoryallocate\n");

	for(i=0;i<numberV;i++){
		dx=vertices[i*3+0];
		dy=vertices[i*3+1];
		dz=vertices[i*3+2];
		setPoint(i,dx,dy,dz);
	}
	for(i=0;i<numberF;i++){
		di=faces[i*3+0];
		dj=faces[i*3+1];
		dk=faces[i*3+2];

		setFace(i,di,dj,dk); 

		IDtool->AppendVF(i,FTail[Face[i][0]]);
		IDtool->AppendVF(i,FTail[Face[i][1]]);
		IDtool->AppendVF(i,FTail[Face[i][2]]);

	}

	SetBoundaryLines();

	SmoothingMax();

	MakeNormals(normal);

	mg1 = new double[vidNum];
	mg2 = new double[vidNum];
	dr1 = new double[vidNum*3];
	dr2 = new double[vidNum*3];

	for(int i=0;i<vidNum;i++){
		SVDFit3Fast(vids[i],vadj[i],vadjNum[i]);
		mg1[i] = k1[vids[i]];
		mg2[i] = k2[vids[i]];

		dr1[i*3+0] = t1[vids[i]]->x;dr1[i*3+1] = t1[vids[i]]->y;dr1[i*3+2] = t1[vids[i]]->z;
		dr2[i*3+0] = t2[vids[i]]->x;dr2[i*3+1] = t2[vids[i]]->y;dr2[i*3+2] = t2[vids[i]]->z;
	}
	//SVDFit3Fast(vadj,vadjNum);
}


void Polyhedron::SaveRmesh(char *fileout){
  int i,j;
  FILE *out = fopen(fileout,"w");
  fprintf(out,"%d\n",numberV);
  fprintf(out,"%d\n",numberF);
  for(i=0;i<numberV;i++){
    fprintf(out,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",k1[i],k2[i],t1[i]->x,t1[i]->y,t1[i]->z,t2[i]->x,t2[i]->y,t2[i]->z,normal[i]->x,normal[i]->y,normal[i]->z);
  }
  fclose(out);
}

void Polyhedron::setRidgeRavine(unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval)
{
  int i,j;
  int numberR=0;
  
  double range=0.0;
  double alpha,beta;
  double ksf0,ksf1,ksf2;
  
  Graph *myg = new Graph(numberV);
  Features *myf = new Features();

  
  Graph *myg2 = new Graph(numberV);
  Features *myf2 = new Features();
  
  int checkID=-1;
  int checkID2=-1;
  IDList *now;
  double dx,dy,dz;
  int edge01=0;
  int edge12=0;
  int edge20=0;

  //Modified Ohake's tracing (Maxima/Minima Test) SIGGRAPH'04
  for(i=0;i<numberF;i++){
    
    
    // Ridge: convex crest 
    if((k1[Face[i][0]])>fabs(k2[Face[i][0]])&&
       (k1[Face[i][1]])>fabs(k2[Face[i][1]])&&
       (k1[Face[i][2]])>fabs(k2[Face[i][2]])){
      
      /* ridge-ravine if(k1[Face[i][0]]>0.0&&k1[Face[i][1]]>0.0&&k1[Face[i][2]]>0.0){ */
      
      /* all extrema if(0==0){ */
      edge01=0;
      edge12=0;
      edge20=0;
      
      ksf0 = ks1[Face[i][0]];
      if(PT->InnerProduct(t1[Face[i][0]],t1[Face[i][1]])<0.0){
	ksf1 = -ks1[Face[i][1]];
	bc[5]->x = -t1[Face[i][1]]->x;
	bc[5]->y = -t1[Face[i][1]]->y;
	bc[5]->z = -t1[Face[i][1]]->z;
	
      }else{
	ksf1 = ks1[Face[i][1]];
	bc[5]->x = t1[Face[i][1]]->x;
	bc[5]->y = t1[Face[i][1]]->y;
	bc[5]->z = t1[Face[i][1]]->z;
      }
      
      if(ksf0*ksf1<=0.0){
	
	PT->makeVector(bc[3],point[Face[i][0]],point[Face[i][1]]);
	PT->makeVector(bc[4],point[Face[i][1]],point[Face[i][0]]);
	if(ksf0*PT->InnerProduct(bc[3],t1[Face[i][0]])>0.0||
	   ksf1*PT->InnerProduct(bc[4],bc[5])>0.0
	   
	   ){
	  edge01=1;
	}
      }
      
      ksf1 = ks1[Face[i][1]];
      if(PT->InnerProduct(t1[Face[i][1]],t1[Face[i][2]])<0.0){
	ksf2 = -ks1[Face[i][2]];
	bc[5]->x = -t1[Face[i][2]]->x;
	bc[5]->y = -t1[Face[i][2]]->y;
	bc[5]->z = -t1[Face[i][2]]->z;
      }else{
	ksf2 = ks1[Face[i][2]];
	bc[5]->x = t1[Face[i][2]]->x;
	bc[5]->y = t1[Face[i][2]]->y;
	bc[5]->z = t1[Face[i][2]]->z;
      }
      
      
      if(ksf1*ksf2<=0.0){
	
	PT->makeVector(bc[3],point[Face[i][1]],point[Face[i][2]]);
	PT->makeVector(bc[4],point[Face[i][2]],point[Face[i][1]]);
	
	
	if(ksf1*PT->InnerProduct(bc[3],t1[Face[i][1]])>0.0||
	   ksf2*PT->InnerProduct(bc[4],bc[5])>0.0){
	  edge12=1;
	  }
      }
      
      ksf2 = ks1[Face[i][2]];
      
      if(PT->InnerProduct(t1[Face[i][2]],t1[Face[i][0]])<0.0){
	ksf0 = -ks1[Face[i][0]];
	bc[5]->x = -t1[Face[i][0]]->x;
	bc[5]->y = -t1[Face[i][0]]->y;
	bc[5]->z = -t1[Face[i][0]]->z;
      }else{
	ksf0 = ks1[Face[i][0]];
	bc[5]->x = t1[Face[i][0]]->x;
	bc[5]->y = t1[Face[i][0]]->y;
	bc[5]->z = t1[Face[i][0]]->z;
      }
      if(ksf2*ksf0<=0.0){
	
	PT->makeVector(bc[3],point[Face[i][2]],point[Face[i][0]]);
	PT->makeVector(bc[4],point[Face[i][0]],point[Face[i][2]]);
	if(ksf2*PT->InnerProduct(bc[3],t1[Face[i][2]])>0.0||
	   ksf0*PT->InnerProduct(bc[4],bc[5])>0.0){
	  edge20=1;
	}
      }
      
      if(edge01+edge12+edge20>=2){
	if(edge01==1){
	  checkID = myg->CheckIsThereEdge(Face[i][0],Face[i][1]);
	  if(checkID==-1){
	    alpha = fabs(ksf0);
	    beta = fabs(ksf1);
	    
	    PT->setInter(bc[0],alpha,beta,point[Face[i][0]],point[Face[i][1]]);
	    
	    checkID = myf->AppendP(myg->AppendE(Face[i][0],Face[i][1]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][0],Face[i][1],alpha,beta);
	  }
	  
	  
	  if(edge12==1){
	    checkID2 = myg->CheckIsThereEdge(Face[i][1],Face[i][2]);
	    if(checkID2==-1){
	      alpha = fabs(ksf1);
	      beta = fabs(ksf2);
	      PT->setInter(bc[0],alpha,beta,point[Face[i][1]],point[Face[i][2]]);
	      checkID2 =myf->AppendP(myg->AppendE(Face[i][1],Face[i][2]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][1],Face[i][2],alpha,beta);
	    }
	    
	  }else{
	    checkID2 = myg->CheckIsThereEdge(Face[i][2],Face[i][0]);
	    if(checkID2==-1){
	      alpha = fabs(ksf2);
	      beta = fabs(ksf0);
	      PT->setInter(bc[0],alpha,beta,point[Face[i][2]],point[Face[i][0]]);
	      checkID2 =myf->AppendP(myg->AppendE(Face[i][2],Face[i][0]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][2],Face[i][0],alpha,beta);
	    }
	    
	    
	  }
	  
	}else{
	  checkID = myg->CheckIsThereEdge(Face[i][1],Face[i][2]);
	  if(checkID==-1){
	    alpha = fabs(ksf1);
	    beta = fabs(ksf2);
	    PT->setInter(bc[0],alpha,beta,point[Face[i][1]],point[Face[i][2]]);
	    checkID =myf->AppendP(myg->AppendE(Face[i][1],Face[i][2]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][1],Face[i][2],alpha,beta);
	  }
	  checkID2 = myg->CheckIsThereEdge(Face[i][2],Face[i][0]);
	  if(checkID2==-1){
	    alpha = fabs(ksf2);
	    beta = fabs(ksf0);
	    PT->setInter(bc[0],alpha,beta,point[Face[i][2]],point[Face[i][0]]);
	    checkID2 =myf->AppendP(myg->AppendE(Face[i][2],Face[i][0]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][2],Face[i][0],alpha,beta);
	    
	    
	  }
	}
	myf->AppendE(checkID,checkID2,i);
      }
    }
    
    

    // Ravine: concave crest
    
    /* all extrema if(0==0){ */
    
    /* ravine if(k2[Face[i][0]]<0.0&&k2[Face[i][1]]<0.0&&k2[Face[i][2]]<0.0){*/
          
    if((-fabs(k1[Face[i][0]]))>(k2[Face[i][0]])&&
       (-fabs(k1[Face[i][1]]))>(k2[Face[i][1]])&&
       (-fabs(k1[Face[i][2]]))>(k2[Face[i][2]])){
      
      
      
      edge01=0;
      edge12=0;
      edge20=0;
      ksf0 = ks2[Face[i][0]];
      if(PT->InnerProduct(t2[Face[i][0]],t2[Face[i][1]])<0.0){
	ksf1 = -ks2[Face[i][1]];
	bc[5]->x = -t2[Face[i][1]]->x;
	bc[5]->y = -t2[Face[i][1]]->y;
	bc[5]->z = -t2[Face[i][1]]->z;
      }else{
	ksf1 = ks2[Face[i][1]];
	bc[5]->x = t2[Face[i][1]]->x;
	bc[5]->y = t2[Face[i][1]]->y;
	bc[5]->z = t2[Face[i][1]]->z;
      }
      
      if(ksf0*ksf1<=0.0){
	
	PT->makeVector(bc[3],point[Face[i][0]],point[Face[i][1]]);
	PT->makeVector(bc[4],point[Face[i][1]],point[Face[i][0]]);
	
	
	
	if(ksf0*PT->InnerProduct(bc[3],t2[Face[i][0]])<0.0||
	   ksf1*PT->InnerProduct(bc[4],bc[5])<0.0	   
	   ){
	  edge01=1;
	  }
      }
      
      ksf1 = ks2[Face[i][1]];
      if(PT->InnerProduct(t2[Face[i][1]],t2[Face[i][2]])<0.0){
	ksf2 = -ks2[Face[i][2]];
	bc[5]->x = -t2[Face[i][2]]->x;
	bc[5]->y = -t2[Face[i][2]]->y;
	bc[5]->z = -t2[Face[i][2]]->z;
      }else{
	ksf2 = ks2[Face[i][2]];
	bc[5]->x = t2[Face[i][2]]->x;
	bc[5]->y = t2[Face[i][2]]->y;
	bc[5]->z = t2[Face[i][2]]->z;
      }
      
      
      if(ksf1*ksf2<=0.0){
	
	PT->makeVector(bc[3],point[Face[i][1]],point[Face[i][2]]);
	PT->makeVector(bc[4],point[Face[i][2]],point[Face[i][1]]);
	
	
	if(ksf1*PT->InnerProduct(bc[3],t2[Face[i][1]])<0.0||
	   ksf2*PT->InnerProduct(bc[4],bc[5])<0.0){
	  edge12=1;
	  }
      }
      
      ksf2 = ks2[Face[i][2]];
      
      if(PT->InnerProduct(t2[Face[i][2]],t2[Face[i][0]])<0.0){
	ksf0 = -ks2[Face[i][0]];
	bc[5]->x = -t2[Face[i][0]]->x;
	bc[5]->y = -t2[Face[i][0]]->y;
	bc[5]->z = -t2[Face[i][0]]->z;
      }else{
	ksf0 = ks2[Face[i][0]];
	bc[5]->x = t2[Face[i][0]]->x;
	bc[5]->y = t2[Face[i][0]]->y;
	bc[5]->z = t2[Face[i][0]]->z;
      }
      
      
      if(ksf2*ksf0<=0.0){
	
	PT->makeVector(bc[3],point[Face[i][2]],point[Face[i][0]]);
	PT->makeVector(bc[4],point[Face[i][0]],point[Face[i][2]]);
	
	
	
	if(ksf2*PT->InnerProduct(bc[3],t2[Face[i][2]])<0.0||
	   ksf0*PT->InnerProduct(bc[4],bc[5])<0.0
	   
	   ){
	  edge20=1;
	  }
      }
      
      if(edge01+edge12+edge20>=2){
	if(edge01==1){
	  checkID = myg2->CheckIsThereEdge(Face[i][0],Face[i][1]);
	  if(checkID==-1){
	    alpha = fabs(ksf0);
	    beta = fabs(ksf1);
	    PT->setInter(bc[0],alpha,beta,point[Face[i][0]],point[Face[i][1]]);
	    
	    
	    
	    checkID = myf2->AppendP(myg2->AppendE(Face[i][0],Face[i][1]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][0],Face[i][1],alpha,beta);
	    
	  }
	  
	  
	  if(edge12==1){
	    checkID2 = myg2->CheckIsThereEdge(Face[i][1],Face[i][2]);
	    if(checkID2==-1){
	      alpha = fabs(ksf1);
	      beta = fabs(ksf2);
	      PT->setInter(bc[0],alpha,beta,point[Face[i][1]],point[Face[i][2]]);
	      
	      checkID2 =myf2->AppendP(myg2->AppendE(Face[i][1],Face[i][2]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][1],Face[i][2],alpha,beta);
	    }
	    
	  }else{
	    checkID2 = myg2->CheckIsThereEdge(Face[i][2],Face[i][0]);
	    if(checkID2==-1){
	      alpha = fabs(ksf2);
	      beta = fabs(ksf0);
	      PT->setInter(bc[0],alpha,beta,point[Face[i][2]],point[Face[i][0]]);
	      
	      checkID2 =myf2->AppendP(myg2->AppendE(Face[i][2],Face[i][0]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][2],Face[i][0],alpha,beta);
	    }
	    
	    
	  }
	  
	}else{
	  checkID = myg2->CheckIsThereEdge(Face[i][1],Face[i][2]);
	  if(checkID==-1){
	    alpha = fabs(ksf1);
	    beta = fabs(ksf2);
	    PT->setInter(bc[0],alpha,beta,point[Face[i][1]],point[Face[i][2]]);
	    
	    checkID =myf2->AppendP(myg2->AppendE(Face[i][1],Face[i][2]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][1],Face[i][2],alpha,beta);
	  }
	  checkID2 = myg2->CheckIsThereEdge(Face[i][2],Face[i][0]);
	  if(checkID2==-1){
	    alpha = fabs(ksf2);
	    beta = fabs(ksf0);
	    PT->setInter(bc[0],alpha,beta,point[Face[i][2]],point[Face[i][0]]);
	    
	    checkID2 =myf2->AppendP(myg2->AppendE(Face[i][2],Face[i][0]),bc[0]->x,bc[0]->y,bc[0]->z,Face[i][2],Face[i][0],alpha,beta);
	    
	    
	  }
	}
	myf2->AppendE(checkID,checkID2,i);
      }
    }
           
  }


  // Start additional connection and cleanning
  int ed1,ed2,ed3;
  
  myf->makeEdgesN();
  myf2->makeEdgesN();
 
  
  // make T-junction
  int checktj=0;
  int checktj2=0;
  
  
  for(i=0;i<numberF;i++){
    ed1 = myg->CheckIsThereEdge(Face[i][0],Face[i][1]);
    ed2 = myg->CheckIsThereEdge(Face[i][1],Face[i][2]);
    ed3 = myg->CheckIsThereEdge(Face[i][2],Face[i][0]);
    
    if(ed1!=-1&&ed2!=-1&&ed3!=-1){
      PT->setCenter(bc[0],point[Face[i][0]],point[Face[i][1]],point[Face[i][2]]);
      checktj=1;
      myf->PushT(ed1,ed2,ed3,bc[0]->x,bc[0]->y,bc[0]->z,i);
    }
    
    ed1 = myg2->CheckIsThereEdge(Face[i][0],Face[i][1]);
    ed2 = myg2->CheckIsThereEdge(Face[i][1],Face[i][2]);
    ed3 = myg2->CheckIsThereEdge(Face[i][2],Face[i][0]);
    
    if(ed1!=-1&&ed2!=-1&&ed3!=-1){
      
      PT->setCenter(bc[0],point[Face[i][0]],point[Face[i][1]],point[Face[i][2]]);
      checktj2=1;
      myf2->PushT(ed1,ed2,ed3,bc[0]->x,bc[0]->y,bc[0]->z,i);
    } 
    
  }
  if(checktj!=0)myf->makeEdgesNIDchange();
  if(checktj2!=0)myf2->makeEdgesNIDchange();
  
  
  
  
  // clean up 
  
  
  
  
  
  // Start modeified Version
  
  /* Clean Up really near points: if neighbor ridge points are closer than 0.0000001 */
  
  myf->RemoveDiplicatesBoundary(0.0000001,bc);
  myf2->RemoveDiplicatesBoundary(0.0000001,bc);
  
  int *check1v = new int[numberV];
  int *check2v = new int[numberV];
  
  
  myf->Connect1RingPoint(numberV,point,bc,check1v,check2v,PT);
  myf2->Connect1RingPoint(numberV,point,bc,check1v,check2v,PT);
  
  myf->Connect1RingEdge((0.5*PI),(PI/3.0),numberV,point,bc,check1v,PT,IHead,ITail);
  
  myf2->Connect1RingEdge((0.5*PI),(PI/3.0),numberV,point,bc,check1v,PT,IHead,ITail);
  
  
  myf->RemoveDiplicatesInner(0.000001);
  myf2->RemoveDiplicatesInner(0.000001);
  
  
  // Filtering
  myf->Labeling();
  myf2->Labeling();
  double shd=0.0;
  // Nor only Ridgeness Ohtake et al. SIGGRAPH'04 but also
  //Window functioned Shapeindex Koenderink, 1990
  double *shapeindex1 = new double[numberV];
  //MVS functional Moreton and Sequin, SIGGRAPH'92
  double *shapeindex2 = new double[numberV];
    
  
  
  for(i=0;i<numberV;i++){
    if(k1[i]-k2[i]<=0.0000000000001){
      shd = 1.0;
    }else{
      
      shd = fabs((2.0/PI)*atan(((k1[i]+k2[i])/(k1[i]-k2[i]))));
    }
    if(shd<=0.5){
      shapeindex1[i] = 2.0*shd;
      
    }else{
      shapeindex1[i] = -2.0*shd+2.0;
    }
        
    
    shapeindex2[i] = sqrt(ks2[i]*ks2[i]+ks1[i]*ks1[i]);
    
        
  }
  
  myf->setRidgeness(k1,shapeindex1,shapeindex2,Face);
  myf2->setRidgeness(k2,shapeindex1,shapeindex2,Face);
    
  
  delete [] shapeindex1;
  delete [] shapeindex2;
  
  delete [] check1v;
  delete [] check2v;

/*
  myf->PrintConnect("ridges.txt",ridge);
  delete myg;
  delete myf;
  myf2->PrintConnect("ravines.txt",ridge);
  delete myg2;
  delete myf2;
*/

  unsigned pNum1,cNum1,eNum1;
  double* fPoints1;
  unsigned* fPointType1;
  double* fFeature1;
  unsigned* fEdges1;
  myf->PrintConnect(pNum1,cNum1,eNum1,fPoints1,fPointType1,fFeature1,fEdges1,ridge);
  delete myg;
  delete myf;

  unsigned pNum2,cNum2,eNum2;
  double* fPoints2;
  unsigned* fPointType2;
  double* fFeature2;
  unsigned* fEdges2;
  myf2->PrintConnect(pNum2,cNum2,eNum2,fPoints2,fPointType2,fFeature2,fEdges2,ridge);
  delete myg2;
  delete myf2;
  
  pNum=pNum1+pNum2;
  cNum=cNum1+cNum2;
  eNum=eNum1+eNum2;
  fPoints=new double[pNum*3];
  fPointType=new unsigned[pNum];
  for(unsigned i=0;i<pNum1;i++){
	  fPoints[i*3+0]=fPoints1[i*3+0];fPoints[i*3+1]=fPoints1[i*3+1];fPoints[i*3+2]=fPoints1[i*3+2];
	  fPointType[i]=fPointType1[i];
  }
  unsigned strp=pNum1*3;
  for(unsigned i=0;i<pNum2;i++){
	  fPoints[strp+i*3+0]=fPoints2[i*3+0];fPoints[strp+i*3+1]=fPoints2[i*3+1];fPoints[strp+i*3+2]=fPoints2[i*3+2];
	  fPointType[pNum1+i]=fPointType2[i]+cNum1;
  }

  fFeature=new double[cNum*3];
  for(unsigned i=0;i<cNum1;i++){
	  fFeature[i*3+0]=fFeature1[i*3+0];fFeature[i*3+1]=fFeature1[i*3+1];fFeature[i*3+2]=fFeature1[i*3+2];
  }
  strp=cNum1*3;
  for(unsigned i=0;i<cNum2;i++){
	  fFeature[strp+i*3+0]=fFeature2[i*3+0];fFeature[strp+i*3+1]=fFeature2[i*3+1];fFeature[strp+i*3+2]=fFeature2[i*3+2];
  }

  fEdges=new unsigned[eNum*3];
  for(unsigned i=0;i<eNum1;i++){
	  fEdges[i*3+0]=fEdges1[i*3+0];fEdges[i*3+1]=fEdges1[i*3+1];fEdges[i*3+2]=fEdges1[i*3+2];
  }
  strp=eNum1*3;
  for(unsigned i=0;i<eNum2;i++){
	  fEdges[strp+i*3+0]=fEdges2[i*3+0]+pNum1;fEdges[strp+i*3+1]=fEdges2[i*3+1]+pNum1;fEdges[strp+i*3+2]=fEdges2[i*3+2];
  }

  interval = eNum1;
}


void Polyhedron::CleanDiplicateEdges(Features *myf,double ep){
  int i;
  if(myf->V==0)return;
  
  myf->setEdgePointerNext();
  int dOV1=-1;
  int dOV2=-1;
  int edgeFID=-1;
  int orivID1=-1;
  int orivID2=-1;
  int neiFID=-1;
  myf->setcheckFace(numberF);
  int checkEremove=0;
  int checkremove=0;
  while(myf->setEdgeIDsNext(&dOV1,&dOV2,&edgeFID,&orivID1,&orivID2)==1){
    if(myf->pBox[dOV1]->tconner!=1&&
       myf->pBox[dOV2]->tconner!=1){
      
      int df1 = get1FaceIDs(edgeFID,Face[edgeFID][0],Face[edgeFID][1]);
      int df2 = get1FaceIDs(edgeFID,Face[edgeFID][1],Face[edgeFID][2]);
      int df3 = get1FaceIDs(edgeFID,Face[edgeFID][2],Face[edgeFID][0]);
      checkremove=0;
      
      if(df1!=-1){
	if(myf->pointerCheck[df1]!=NULL){
	  if(myf->isDiplicateEdge(edgeFID,df1,ep)==1)
	    checkremove=1;
	}
	
      }
      if(df2!=-1){
	if(myf->pointerCheck[df2]!=NULL){
	    if(myf->isDiplicateEdge(edgeFID,df2,ep)==1)checkremove=1;
	}
	
      }
      if(df3!=-1){
	if(myf->pointerCheck[df3]!=NULL){
	  
	    if(myf->isDiplicateEdge(edgeFID,df3,ep)==1)checkremove=1;
	}
	
      }
      if(checkremove==1){
	
	myf->RemoveEdge(edgeFID);
	checkEremove=1;
      }
    }
  }
  
  if(checkEremove==1){
    myf->RemoveV();
    myf->makeEdgesNIDchange(); 
  }
  myf->DeletecheckFace();
  
}


   int Polyhedron::get1FaceIDs(int tFID1,int ov1,int ov2){
     IDList *now =  NULL;
     if(boundary[ov1]==1&&boundary[ov2]==1)return -1;
     now =  FHead[ov1];
     while(next(now)!=FTail[ov1]){
       now = next(now);
       if(now->ID!=((tFID1))&&((Face[now->ID][0]==ov1&&Face[now->ID][1]==ov2)||
			       (Face[now->ID][0]==ov2&&Face[now->ID][1]==ov1)||
			       (Face[now->ID][1]==ov1&&Face[now->ID][2]==ov2)||
			       (Face[now->ID][1]==ov2&&Face[now->ID][2]==ov1)||
			       (Face[now->ID][2]==ov1&&Face[now->ID][0]==ov2)||
			       (Face[now->ID][2]==ov2&&Face[now->ID][0]==ov1))){
	 return now->ID;
       }
     }
     return -1;

   }
  void Polyhedron::set2FaceIDs(int *tFID1,int *tFID2,int ov1,int ov2){
    IDList *now =  NULL;
    if(boundary[ov1]==1&&boundary[ov2]==1){
      (*tFID2)=-1;
      now =  FHead[ov1];
      while(next(now)!=FTail[ov1]){
	now = next(now);
	if((Face[now->ID][0]==ov1&&Face[now->ID][1]==ov2)||
	   (Face[now->ID][0]==ov2&&Face[now->ID][1]==ov1)||
	   (Face[now->ID][1]==ov1&&Face[now->ID][2]==ov2)||
	   (Face[now->ID][1]==ov2&&Face[now->ID][2]==ov1)||
	   (Face[now->ID][2]==ov1&&Face[now->ID][0]==ov2)||
	  (Face[now->ID][2]==ov2&&Face[now->ID][0]==ov1)){
	  (*tFID1) = now->ID;
	  return;
	}
      }
    }else{
      
      now =  FHead[ov1];
      while(next(now)!=FTail[ov1]){
	now = next(now);
	if((Face[now->ID][0]==ov1&&Face[now->ID][1]==ov2)||
	   (Face[now->ID][0]==ov2&&Face[now->ID][1]==ov1)||
	   (Face[now->ID][1]==ov1&&Face[now->ID][2]==ov2)||
	   (Face[now->ID][1]==ov2&&Face[now->ID][2]==ov1)||
	   (Face[now->ID][2]==ov1&&Face[now->ID][0]==ov2)||
	  (Face[now->ID][2]==ov2&&Face[now->ID][0]==ov1)){
	  (*tFID1) = now->ID;
	  break;
	}
      }
      now =  FHead[ov1];
      while(next(now)!=FTail[ov1]){
	now = next(now);
	if(now->ID!=((*tFID1))&&((Face[now->ID][0]==ov1&&Face[now->ID][1]==ov2)||
	   (Face[now->ID][0]==ov2&&Face[now->ID][1]==ov1)||
	   (Face[now->ID][1]==ov1&&Face[now->ID][2]==ov2)||
	   (Face[now->ID][1]==ov2&&Face[now->ID][2]==ov1)||
	   (Face[now->ID][2]==ov1&&Face[now->ID][0]==ov2)||
	  (Face[now->ID][2]==ov2&&Face[now->ID][0]==ov1))){
	  (*tFID2) = now->ID;
	  return;
	}
      }
          
    }

  }


/* smoothing via neighbor centroids */
void Polyhedron::SmoothingMax(){
  int i,j;
  IDList *now=NULL;
 
  for(i=0;i<numberV;i++){
    centroid[i]->x = point[i]->x;
    centroid[i]->y = point[i]->y;
    centroid[i]->z = point[i]->z;
  }
  
  for(i=0;i<numberV;i++){
    point[i]->x=0.0;
    point[i]->y=0.0;
    point[i]->z=0.0;
    
    
    now = VHead[i];
    while(next(now)!=VTail[i]){
      now = next(now);
      
      PT->setCenter(bc[0],centroid[i],centroid[now->ID],centroid[next(now)->ID]);
      
      
      point[i]->x += (bc[0]->x);
      point[i]->y += (bc[0]->y);
      point[i]->z += (bc[0]->z);
      
      now = next(now);
    }
    if(neighborF[i]!=0){
      point[i]->x /= ((double)(neighborF[i]));
      point[i]->y /= ((double)(neighborF[i]));
      point[i]->z /= ((double)(neighborF[i]));
    } 
    
  }
  
}




/* N. Max, Graphical Tools 4(2), normal approximation */
void Polyhedron::MakeNormals(Point3d **dNorma){
  int i,j;
  IDList *now=NULL;
  double dummytemp=0.0;
  double angle=0.0;
  double dsize1=0.0;
  double dsize2=0.0;
  double weight=0.0;
  for(i=0;i<numberV;i++){
    now = VHead[i];
    bc[3]->x=0.0;
    bc[3]->y=0.0;
    bc[3]->z=0.0;
    
    while(next(now)!=VTail[i]){
      now = next(now);
      PT->makeVector(bc[0],point[i],point[next(now)->ID]);dsize1=PT->Point3dSize(bc[0]);if(dsize1==0.0)dsize1=1.0;
      PT->makeVector(bc[1],point[i],point[now->ID]);dsize2=PT->Point3dSize(bc[1]);if(dsize2==0.0)dsize2=1.0;
      
      PT->CrossVector(bc[2],bc[1],bc[0]);
      weight = 1.0/(dsize1*dsize1*dsize2*dsize2);
      bc[3]->x += weight*bc[2]->x;
      bc[3]->y += weight*bc[2]->y;
      bc[3]->z += weight*bc[2]->z;
      

      now = next(now);
    }
    dummytemp = PT->Point3dSize(bc[3]);
    if(dummytemp != 0.0){
      dNorma[i]->x = ((bc[3]->x)/dummytemp);
      dNorma[i]->y = ((bc[3]->y)/dummytemp);
      dNorma[i]->z = ((bc[3]->z)/dummytemp);
    }
    
  }
}


/* normal-enhanced cubic fitting: Goldfeather and Interrante TOG 23(1) */
void Polyhedron::SVDFit3Fast(int**& vadj,int*& vadjNum){
/*
	for(int vid=0;vid<numberV;vid++){
		SVDFit3Fast(vid,vadj[vid],vadjNum[vid]);
	}
*/

  int i,j;
  IDList *now;
  
  
  SvdSolve *mysvd = new SvdSolve();
  double dz = 0.0;
  double dx=0.0;
  double dy=0.0;
  
  int sN=8;
  double *w = new double[sN];
  double *ab = new double[sN];
  double **v = new double* [sN];
  for(i=0;i<sN;i++)v[i] = new double[sN];
  double **a=NULL;
  double *dkk = NULL;
  
  double dsize1=0.0;
  double dsize2=0.0;
  double wmax=0.0;
  double wmin=0.0;
  
  
  double unb=0.0;
  double eid[3];
  double **eia = new double* [3];
  double **eiv = new double* [3];
  int nrot=0;
  Eigens *myeigen = new Eigens();
  for(i=0;i<3;i++){
    eia[i] = new double[3];
    eiv[i] = new double[3];
  }
  
  
  Point3d *bt1 = new Point3d(0.0,0.0,0.0);
  Point3d *bt2 = new Point3d(0.0,0.0,0.0);
  double dUU=0.0;
  double dVV=0.0;
  double UU,UV,VV;
  double dist=0.0;
  double acc=0.0;
  double dw=0.0;

  //my test code
/*
  for(i=0;i<numberV;i++){
	  if(vadjNum[i]!=neighborIN[i])
		printf("i:%d maj:%d caj:%d\t",i,vadjNum[i],neighborIN[i]);
  }
*/
  //end of test

  for(i=0;i<numberV;i++){

/*
	  printf("i=%d  ",i);
		for(int idn=0;idn<vadjNum[i];idn++){
			printf("%f ",vadjWt[i][idn]);		
		}
*/


//change code:   
/*
    a = new double* [3*neighborIN[i]+1];
    for(j=0;j<3*neighborIN[i]+1;j++)a[j] = new double[sN];
    dkk = new double[3*neighborIN[i]+1];
*/
	  if(vadj==NULL){
		  a = new double* [3*neighborIN[i]+1];
		  for(j=0;j<3*neighborIN[i]+1;j++)a[j] = new double[sN];
		  dkk = new double[3*neighborIN[i]+1];
	  }
	  else{
		  a = new double* [3*vadjNum[i]+1];
		  for(j=0;j<3*vadjNum[i]+1;j++)a[j] = new double[sN];
		  dkk = new double[3*vadjNum[i]+1];
	  }
//end of change

    j=0;
    PT->setVitrualTangents(normal[i],bt1,bt2);

  
//change code:

	int count = 0;
	if(vadj==NULL){
		 now = INHead[i];
		while(next(now)!=INTail[i]){
		  now = next(now);
      
		  PT->makeVector(bc[0],point[i],point[now->ID]);
		  dz = PT->InnerProduct(bc[0],normal[i]);
      
		  dUU = PT->InnerProduct(bc[0],bt1);
		  dVV = PT->InnerProduct(bc[0],bt2);
		  UU = dUU*dUU;
		  UV = dUU*dVV;
		  VV = dVV*dVV;
      
  
		  a[3*j+1][1] = 0.5*UU;
		  a[3*j+1][2] = UV;
		  a[3*j+1][3] = 0.5*VV;
		  a[3*j+1][4] = dUU*UU;
		  a[3*j+1][5] = UU*dVV;
		  a[3*j+1][6] = dUU*VV;
		  a[3*j+1][7] = VV*dVV;
      
       
      
      
		  dkk[3*j+1] = dz;
     
      
		  a[3*j+2][1] = dUU;
		  a[3*j+2][2] = dVV;
		  a[3*j+2][3] = 0.0;
		  a[3*j+2][4] = 3.0*UU;
		  a[3*j+2][5] = 2.0*UV;
		  a[3*j+2][6] = VV;
		  a[3*j+2][7] = 0.0;
      
      
		  dx = PT->InnerProduct(normal[now->ID],bt1);
		  dz = PT->InnerProduct(normal[now->ID],normal[i]);
		  dy = PT->InnerProduct(normal[now->ID],bt2);
      
		  dkk[3*j+2] = -dx/dz;
      
		  a[3*j+3][1] = 0.0;
		  a[3*j+3][2] = dUU;
		  a[3*j+3][3] = dVV;
		  a[3*j+3][4] = 0.0;
		  a[3*j+3][5] = UU;
		  a[3*j+3][6] = 2.0*UV;
		  a[3*j+3][7] = 3.0*VV;
      
      
      
		  dkk[3*j+3] = -dy/dz;
      

		  j++;
		}
	}
	else{
		while(count<vadjNum[i]){
			int adjNode = vadj[i][count];

			PT->makeVector(bc[0],point[i],point[adjNode]);
			dz = PT->InnerProduct(bc[0],normal[i]);

			dUU = PT->InnerProduct(bc[0],bt1);
			dVV = PT->InnerProduct(bc[0],bt2);
			UU = dUU*dUU;
			UV = dUU*dVV;
			VV = dVV*dVV;


			a[3*j+1][1] = 0.5*UU;
			a[3*j+1][2] = UV;
			a[3*j+1][3] = 0.5*VV;
			a[3*j+1][4] = dUU*UU;
			a[3*j+1][5] = UU*dVV;
			a[3*j+1][6] = dUU*VV;
			a[3*j+1][7] = VV*dVV;




			dkk[3*j+1] = dz;


			a[3*j+2][1] = dUU;
			a[3*j+2][2] = dVV;
			a[3*j+2][3] = 0.0;
			a[3*j+2][4] = 3.0*UU;
			a[3*j+2][5] = 2.0*UV;
			a[3*j+2][6] = VV;
			a[3*j+2][7] = 0.0;


			dx = PT->InnerProduct(normal[adjNode],bt1);
			dz = PT->InnerProduct(normal[adjNode],normal[i]);
			dy = PT->InnerProduct(normal[adjNode],bt2);

			dkk[3*j+2] = -dx/dz;

			a[3*j+3][1] = 0.0;
			a[3*j+3][2] = dUU;
			a[3*j+3][3] = dVV;
			a[3*j+3][4] = 0.0;
			a[3*j+3][5] = UU;
			a[3*j+3][6] = 2.0*UV;
			a[3*j+3][7] = 3.0*VV;



			dkk[3*j+3] = -dy/dz;

			count++;
			j++;
		}
	}
//end of change


//change code:
	//mysvd->svdcmp(a,3*neighborIN[i],(sN-1),w,v);
	if(vadj==NULL)
		mysvd->svdcmp(a,3*neighborIN[i],(sN-1),w,v);
	else
		mysvd->svdcmp(a,3*vadjNum[i],(sN-1),w,v);
//end of change

    wmax=0.0;
    for(j=1;j<=(sN-1);j++){
      if(w[j]>wmax)wmax=w[j];
    }
    wmin = wmax*0.000001;
    for(j=1;j<=(sN-1);j++)if(w[j]<wmin)w[j]=0.0;
	
//change code:
	//mysvd->svbksb(a,w,v,3*neighborIN[i],(sN-1),dkk,ab);
	if(vadj==NULL)
		mysvd->svbksb(a,w,v,3*neighborIN[i],(sN-1),dkk,ab);
	else
		mysvd->svbksb(a,w,v,3*vadjNum[i],(sN-1),dkk,ab);
//end of change
   
    double fxx = ab[1];
    double fxy = ab[2];
    double fyy = ab[3];
    eia[1][1] = ab[1];
    eia[1][2] = ab[2];
    eia[2][1] = ab[2];
    eia[2][2] = ab[3];
    myeigen->jacobi(eia,2,eid,eiv,&nrot);
    ab[1] = fxx;
    ab[2] = fxy;
    ab[3] = fyy;
    
    double t1x,t1y,t2x,t2y;

    if(eid[1]<eid[2]){
      k1[i] = eid[2];
      k2[i] = eid[1];
      t2[i]->x = (bt1->x)*eiv[1][1] + (bt2->x)*eiv[2][1];
      t2[i]->y = (bt1->y)*eiv[1][1] + (bt2->y)*eiv[2][1];
      t2[i]->z = (bt1->z)*eiv[1][1] + (bt2->z)*eiv[2][1];
      t1[i]->x = (bt1->x)*eiv[1][2] + (bt2->x)*eiv[2][2];
      t1[i]->y = (bt1->y)*eiv[1][2] + (bt2->y)*eiv[2][2];
      t1[i]->z = (bt1->z)*eiv[1][2] + (bt2->z)*eiv[2][2];
      
      t1x = eiv[1][2];
      t1y = eiv[2][2];
      t2x = eiv[1][1];
      t2y = eiv[2][1];
    }else{
      k1[i] = eid[1];
      k2[i] = eid[2];
            
      t1[i]->x = (bt1->x)*eiv[1][1] + (bt2->x)*eiv[2][1];
      t1[i]->y = (bt1->y)*eiv[1][1] + (bt2->y)*eiv[2][1];
      t1[i]->z = (bt1->z)*eiv[1][1] + (bt2->z)*eiv[2][1];
      t2[i]->x = (bt1->x)*eiv[1][2] + (bt2->x)*eiv[2][2];
      t2[i]->y = (bt1->y)*eiv[1][2] + (bt2->y)*eiv[2][2];
      t2[i]->z = (bt1->z)*eiv[1][2] + (bt2->z)*eiv[2][2];
      t2x = eiv[1][2];
      t2y = eiv[2][2];
      t1x = eiv[1][1];
      t1y = eiv[2][1]; 
    }
    PT->Normalize3D(t1[i]);
    PT->Normalize3D(t2[i]);
    
    // computing first derivative of principal curvature w.r.t their direction
    ks1[i] = getk1t1(t1x,t1y);
    ks2[i] = getk1t1(t2x,t2y);

    delete [] dkk;dkk=NULL;

//change code:
	//for(j=0;j<3*neighborIN[i]+1;j++)delete a[j];
	if(vadj==NULL)
		for(j=0;j<3*neighborIN[i]+1;j++)delete a[j];
	else
		for(j=0;j<3*vadjNum[i]+1;j++)delete a[j];
//end of change

    delete [] a;
    a = NULL;
  }
  delete mysvd;
  delete [] w;
  delete [] ab;
  for(i=0;i<sN;i++)delete v[i];
  for(i=0;i<3;i++){
    delete eia[i];
    
    delete eiv[i];
  }

  delete [] eiv;
  delete [] eia;
  delete [] v;
  delete myeigen;
  delete bt1;
  delete bt2;
}
void Polyhedron::SVDFit3Fast(int& vid,int*& vadj,int& vadjNum){
  int i,j;
  
  SvdSolve *mysvd = new SvdSolve();
  double dz = 0.0;
  double dx=0.0;
  double dy=0.0;
  
  int sN=8;
  double *w = new double[sN];
  double *ab = new double[sN];
  double **v = new double* [sN];
  for(i=0;i<sN;i++)v[i] = new double[sN];
  double **a=NULL;
  double *dkk = NULL;
  
  double dsize1=0.0;
  double dsize2=0.0;
  double wmax=0.0;
  double wmin=0.0;
  
  
  double unb=0.0;
  double eid[3];
  double **eia = new double* [3];
  double **eiv = new double* [3];
  int nrot=0;
  Eigens *myeigen = new Eigens();
  for(i=0;i<3;i++){
    eia[i] = new double[3];
    eiv[i] = new double[3];
  }  
  
  Point3d *bt1 = new Point3d(0.0,0.0,0.0);
  Point3d *bt2 = new Point3d(0.0,0.0,0.0);
  double dUU=0.0;
  double dVV=0.0;
  double UU,UV,VV;
  double dist=0.0;
  double acc=0.0;
  double dw=0.0;

  ///////////////////////////
	 a = new double* [3*vadjNum+1];
	for(j=0;j<3*vadjNum+1;j++)a[j] = new double[sN];
		  dkk = new double[3*vadjNum+1];

     i=vid; j=0;
    PT->setVitrualTangents(normal[i],bt1,bt2);

	while(j<vadjNum){
		int adjNode = vadj[j];

		PT->makeVector(bc[0],point[i],point[adjNode]);
		dz = PT->InnerProduct(bc[0],normal[i]);

		dUU = PT->InnerProduct(bc[0],bt1);
		dVV = PT->InnerProduct(bc[0],bt2);
		UU = dUU*dUU;
		UV = dUU*dVV;
		VV = dVV*dVV;


		a[3*j+1][1] = 0.5*UU;
		a[3*j+1][2] = UV;
		a[3*j+1][3] = 0.5*VV;
		a[3*j+1][4] = dUU*UU;
		a[3*j+1][5] = UU*dVV;
		a[3*j+1][6] = dUU*VV;
		a[3*j+1][7] = VV*dVV;




		dkk[3*j+1] = dz;


		a[3*j+2][1] = dUU;
		a[3*j+2][2] = dVV;
		a[3*j+2][3] = 0.0;
		a[3*j+2][4] = 3.0*UU;
		a[3*j+2][5] = 2.0*UV;
		a[3*j+2][6] = VV;
		a[3*j+2][7] = 0.0;


		dx = PT->InnerProduct(normal[adjNode],bt1);
		dz = PT->InnerProduct(normal[adjNode],normal[i]);
		dy = PT->InnerProduct(normal[adjNode],bt2);

		dkk[3*j+2] = -dx/dz;

		a[3*j+3][1] = 0.0;
		a[3*j+3][2] = dUU;
		a[3*j+3][3] = dVV;
		a[3*j+3][4] = 0.0;
		a[3*j+3][5] = UU;
		a[3*j+3][6] = 2.0*UV;
		a[3*j+3][7] = 3.0*VV;



		dkk[3*j+3] = -dy/dz;

		j++;
	}

	mysvd->svdcmp(a,3*vadjNum,(sN-1),w,v);

    wmax=0.0;
    for(j=1;j<=(sN-1);j++){
      if(w[j]>wmax)wmax=w[j];
    }
    wmin = wmax*0.000001;
    for(j=1;j<=(sN-1);j++)if(w[j]<wmin)w[j]=0.0;
	
	mysvd->svbksb(a,w,v,3*vadjNum,(sN-1),dkk,ab);

    
    double fxx = ab[1];
    double fxy = ab[2];
    double fyy = ab[3];
    eia[1][1] = ab[1];
    eia[1][2] = ab[2];
    eia[2][1] = ab[2];
    eia[2][2] = ab[3];
    myeigen->jacobi(eia,2,eid,eiv,&nrot);
    ab[1] = fxx;
    ab[2] = fxy;
    ab[3] = fyy;
    
    double t1x,t1y,t2x,t2y;
    
    if(eid[1]<eid[2]){
      k1[i] = eid[2];
      k2[i] = eid[1];
      t2[i]->x = (bt1->x)*eiv[1][1] + (bt2->x)*eiv[2][1];
      t2[i]->y = (bt1->y)*eiv[1][1] + (bt2->y)*eiv[2][1];
      t2[i]->z = (bt1->z)*eiv[1][1] + (bt2->z)*eiv[2][1];
      t1[i]->x = (bt1->x)*eiv[1][2] + (bt2->x)*eiv[2][2];
      t1[i]->y = (bt1->y)*eiv[1][2] + (bt2->y)*eiv[2][2];
      t1[i]->z = (bt1->z)*eiv[1][2] + (bt2->z)*eiv[2][2];
      
      t1x = eiv[1][2];
      t1y = eiv[2][2];
      t2x = eiv[1][1];
      t2y = eiv[2][1];
    }else{
      k1[i] = eid[1];
      k2[i] = eid[2];
            
      t1[i]->x = (bt1->x)*eiv[1][1] + (bt2->x)*eiv[2][1];
      t1[i]->y = (bt1->y)*eiv[1][1] + (bt2->y)*eiv[2][1];
      t1[i]->z = (bt1->z)*eiv[1][1] + (bt2->z)*eiv[2][1];
      t2[i]->x = (bt1->x)*eiv[1][2] + (bt2->x)*eiv[2][2];
      t2[i]->y = (bt1->y)*eiv[1][2] + (bt2->y)*eiv[2][2];
      t2[i]->z = (bt1->z)*eiv[1][2] + (bt2->z)*eiv[2][2];
      t2x = eiv[1][2];
      t2y = eiv[2][2];
      t1x = eiv[1][1];
      t1y = eiv[2][1]; 
    }
    PT->Normalize3D(t1[i]);
    PT->Normalize3D(t2[i]);
    
    // computing first derivative of principal curvature w.r.t their direction
    ks1[i] = getk1t1(t1x,t1y);
    ks2[i] = getk1t1(t2x,t2y);
        
    
    delete [] dkk;dkk=NULL;

	for(j=0;j<3*vadjNum+1;j++)delete a[j];

    delete [] a;
    a = NULL;

  delete mysvd;
  delete [] w;
  delete [] ab;
  for(i=0;i<sN;i++)delete v[i];
  for(i=0;i<3;i++){
    delete eia[i];
    
    delete eiv[i];
  }

  delete [] eiv;
  delete [] eia;
  delete [] v;
  delete myeigen;
  delete bt1;
  delete bt2;

}

void Polyhedron::OriginalCoordinate(){
  int i;
  for(i=0;i<numberV;i++){
    point[i]->x = centroid[i]->x;
    point[i]->y = centroid[i]->y;
    point[i]->z = centroid[i]->z;
  }
  
}

void Polyhedron::setProperty1(int**& vadj,int*& vadjNum,unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval)
{
 // clock_t timestr = clock();

 // clock_t time1 = clock();
  SmoothingMax();
 // time1 = clock() - time1;

 // clock_t time2 = clock();
  MakeNormals(normal);
 // time2 = clock() - time2;
   
 // clock_t time3 = clock();
  if(vadj==NULL) setNeighbor(); //this part takes most of times when k increases.
  //time3 = clock() - time3;

 // clock_t time4 = clock();
  SVDFit3Fast(vadj,vadjNum);
  //time4 = clock() - time4;

 // clock_t time5 = clock();
  OriginalCoordinate();
 // time5 = clock() - time5;

 // clock_t time6 = clock();
  if(ridge==1)
    setRidgeRavine(pNum,cNum,eNum,fPoints,fPointType,fFeature,fEdges,interval);
 // time6 = clock() - time6;
    
 // timestr = clock() - timestr;

  //printf("compute crestline. time(%d) t1(%d) t2(%d) t3(%d) t4(%d) t5(%d) t6(%d) \n",timestr,time1,time2,time3,time4,time5,time6);
}

  
int Polyhedron::setNring(int clevel,int maxlevel,int cID,int labelID){
  
  if(maxlevel<=clevel)return 0;
  if(cID!=labelID){
    if(label[cID]!=labelID){
      if(PT->InnerProduct(normal[labelID],normal[cID])>0.0){
      	IDtool->AppendVF(cID,INTail[labelID]);
	neighborIN[labelID]++;
      }
    }
    label[cID]=labelID;
    IDList *now = IHead[cID];
    while(next(now)!=ITail[cID]){
      now = next(now);
      setNring((clevel+1),maxlevel,now->ID,labelID);
    }
  }else{
    if(label[cID]!=labelID){
      label[cID]=labelID;
      IDList *now = IHead[cID];
      while(next(now)!=ITail[cID]){
	now = next(now);
	setNring((clevel+1),maxlevel,now->ID,labelID);
      }
    }
  
  }
  
  return 0;
}
void Polyhedron::setNeighbor(){
  int i=0;
  label = new int[numberV];
  for(i=0;i<numberV;i++){
    label[i] = -1;
  }
  
  for(i=0;i<numberV;i++){
    neighborIN[i]=0;
    setNring(0,MAXNEIGHBORLAVEL,i,i);
    
    if(neighborIN[i]==0){
      IDList *now = IHead[i];
      while(next(now)!=ITail[i]){
	now = next(now);
	IDtool->AppendVF(now->ID,INTail[i]);
	neighborIN[i]++;
      }
    }
    
  }
  delete [] label;
}


void Polyhedron::memorydelete(){
  int i;
  
  
  for(i=0;i<10;i++){
    delete bc[i];
    delete dbc[i];
  }
  delete [] bc;
  bc = NULL;
  delete [] dbc;
  dbc = NULL;
  
  
  if(IDtool!=NULL){
    if(FHead!=NULL){
      IDtool->CleanNeighborLL(FHead,FTail,numberV,neighborF);
          
    }
    if(INHead!=NULL){
      IDtool->CleanNeighborLL(INHead,INTail,numberV,neighborIN);
      
    }
  if(IHead!=NULL&&ITail!=NULL){
      IDtool->CleanNeighborL(IHead,ITail,numberV);
      IHead=NULL;
      ITail=NULL;
    }
    if(VHead!=NULL&&VTail!=NULL){
      IDtool->CleanNeighborL(VHead,VTail,numberV);
      VHead=NULL;
      VTail=NULL;
    }
      
    
    delete IDtool;
    IDtool=NULL; 
  }

  if(point!=NULL){
    if(numberV!=0){
      for(i=0;i<numberV;i++){
	delete point[i];
	delete centroid[i];
	delete t1[i];
	delete t2[i];
	delete normal[i];
	
      }
    }
    delete [] ks1;
    delete [] ks2;
    
    
    delete [] k1;
    delete [] k2;
        
    delete [] t1;
    delete [] t2;
    t1 = NULL;
    t2 = NULL;
    delete [] point;
    point = NULL;
    
    delete [] normal;
    delete [] centroid;
    
    normal = NULL;
  }
  if(Face!=NULL){
    if(numberF!=0){
      for(i=0;i<numberF;i++)delete [] Face[i];
    }
    delete [] Face;
    Face=NULL;
  }
  if(neighborI!=NULL){
    delete [] neighborI;
  }
  /* bug
 if(neighborF!=NULL){
    delete [] neighborF;
  }
  */
 if(boundary!=NULL){
    delete [] boundary;
  }

 delete PT;

 
}
  

