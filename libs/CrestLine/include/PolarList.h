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
#define next(p) p->next
#define back(p) p->back

class PolarList{
 public:
  int ID;
  double r;
  double theta;
  double lambda;
  double sigma;
  double cotw;
  //double omega;
  //int backID;
  //int nextID;
  PolarList *next;
  PolarList *back;
  //int nextID;
  //PolarList *nextPP;
  int tR,tL;
  

  //PolarList *Dnext;
  //PolarList *Dback;
  //int right,left;
  PolarList(){next=NULL;back=NULL;}
  PolarList(int dv){
    //Dback=NULL;Dnext=NULL;
    ID=dv;next=NULL;back=NULL;r=0.0;theta=0.0;lambda=0.0;}
  PolarList(int dv,double dx,double dy){
    //Dback=NULL;Dnext=NULL;
    ID=dv;next=NULL;back=NULL;r=dx;theta=dy;lambda=0.0;}
   PolarList(int dv,double dx,double dy,double dz,double cotwij){
     //Dback=NULL;Dnext=NULL;
     ID=dv;next=NULL;back=NULL;r=dx;theta=dy;lambda=dz;cotw=cotw;}
   /*
    PolarList(int dv,double dx,double dy,double dz,int dnextID,int dbackID){
     nextID = dnextID;backID = dbackID;
     ID=dv;next=NULL;back=NULL;r=dx;theta=dy;lambda=dz;}
   */
   PolarList(int dv,double dx,double dy,double dz,int dR,int dL,int dnextID){
     
     //nextID = dnextID;
     ID=dv;next=NULL;back=NULL;r=dx;theta=dy;lambda=dz;tR = dR;tL = dL;//omega=0.0;nextPP=NULL;
   }
   
  virtual ~PolarList(){ID = -1; next=NULL;back=NULL;}
 private:
  PolarList(const PolarList& rhs);
  const PolarList &operator=(const PolarList& rhs);
   
};
