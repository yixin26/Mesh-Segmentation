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

class Point3d;
class IDList;
class PolarList;
class IDSet;
class PointTool;
class Features;
class Polyhedron{
 public:
  int numberV;
  int numberF;
  int featurestrip;
  int featureconnect1ring;
  Point3d **point;
  int MAXNEIGHBORLAVEL;
  int **Face;
  double *ks1;
  double *ks2;
    
  double *k1;
  double *k2;
  Point3d **t1;
  Point3d **t2;
  Point3d **normal;
  Point3d **centroid;
  
  
    
  IDList **IHead;
  IDList **ITail;
  IDList **INHead;
  IDList **INTail;
  int *neighborIN;

  IDList **VHead;
  IDList **VTail;
  IDList **FHead;
  IDList **FTail;
  
  IDSet *IDtool;
  int *boundary;
  
  int *neighborI;
  int *neighborF;
  
    
  int numboundary;
  Point3d **bc;
  Point3d **dbc;
  PointTool *PT;
  int ridge;
  int *label;
  
  Polyhedron(){
    featurestrip=1;featureconnect1ring=1;
    
    ridge=1;
    
    MAXNEIGHBORLAVEL=1;
    point=NULL;
    
    bc=NULL;
    Face=NULL;
    k1=NULL;
    k2=NULL;
    normal=NULL;
    t1=NULL;
    t2=NULL;
    numberV=0;
    numberF=0;
    
    IHead=NULL;
    ITail=NULL;
    IDtool=NULL;
    
    PT=new PointTool();
    
    FHead=NULL;
    FTail=NULL;
    
    neighborF=NULL;
    
  }
  virtual ~Polyhedron(){
    if(numberV!=0&&point!=NULL){
      memorydelete();
    }
  }
  
  int setNring(int clevel,int maxlevel,int cID,int labelID);
  void setNeighbor();
  void set2FaceIDs(int *tFID1,int *tFID2,int ov1,int ov2);
  int get1FaceIDs(int tFID1,int ov1,int ov2);
  void CleanDiplicateEdges(Features *myf,double ep);
  void setRidgeRavine(unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval);
  
  void SmoothingMax();
  
  
  void OriginalCoordinate();
  
  void MakeNormals(Point3d **dNorma);
  
  
  
  
  void setPoint(int i,double dx,double dy,double dz);
  void setFace(int i,int di,int dj,int dk);
 
  void q_sort( double *a , const long left , const long right )
{
        long    somewhere , i , j;
        double pivot,work;
	if( left < right ) {
                somewhere = ( left + right ) / 2;
                pivot = a[ somewhere ];
                i = left; j = right;
                do
                {
                        while( a[ i ] < pivot ) i++;
                        while( a[ j ] > pivot ) j--;
                        if( i <= j ) {
                                work = a[ j ];
                                a[ j-- ] = a[ i ];
                                a[ i++ ] = work;
                        }
                }
                 while( i <= j );
                 q_sort( a , left , j );
                 q_sort( a , i , right );
        }
}

  void setProperty1(int**& vadj,int*& vadjNum,unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval);
  
  void SVDFit3Fast(int**& vadj,int*& vadjNum);
  void SVDFit3Fast(int& vid,int*& vadj,int& vadjNum);

  void readmesh(int** vadj,int* vadjNum,double* vertices, unsigned vNum, unsigned* faces, unsigned fNum, unsigned neighboreSize, unsigned ridgeTension,
	  unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval);
  
  //void readmesh(char *filename,char *fileout);

  void SaveRmesh(char *fileout);
  void SetBoundaryLines();

  void getCurvature(int*& vids,unsigned vidNum,int**& vadj,int*& vadjNum,double*& vertices, unsigned vNum, unsigned*& faces, unsigned fNum, unsigned neighboreSize,double*& mg1,double*& mg2,double*& dr1,double*& dr2);
  
 private: 
  
  void memoryallocate(int, int );
  void memorydelete();
  
  Polyhedron(const Polyhedron& rhs);
  const Polyhedron &operator=(const Polyhedron& rhs);
  
};

