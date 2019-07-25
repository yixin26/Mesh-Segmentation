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
class Edge{
 public:
  int ID;
  Edge *next;
  Edge *back;  
  int v1,v2;
  int fID;
  int ff1,ff2;
  Edge(){
    v1=-1;
    v2=-1;
    ID=-1;
    fID=-1;
    ff1=-1;
    ff2=-1;
  }
  Edge(int dv1,int dv2){
    v1 = dv1;
    v2 = dv2;
    ID=-1;
    fID=-1;
    ff1=-1;
    ff2=-1;
  }
  Edge(int dID,int dv1,int dv2){
    v1 = dv1;
    v2 = dv2;
    ID=dID;
    fID=-1;
    ff1=-1;
    ff2=-1;
  }
  
  
  virtual ~Edge(){}

 private:
  Edge(const Edge& rhs);
  const Edge &operator=(const Edge& rhs);
  
};

class EList{
 public:
  int ID;
  Edge *e;
  EList *next;
  EList *back;
  EList(){
    ID=-1;
    e=NULL;
  }
  EList(int dID,Edge *de){
    e = de;
    ID = dID;
  }
  virtual ~EList(){}
 private:
  EList(const EList& rhs);
  const EList &operator=(const EList& rhs);

};

class Graph{
 public:
  
  int E;
  int numberV;
  Edge *EHead;
  Edge *ETail;
  EList **ELhead;
  EList **ELtail;
  
  Graph(int dV){
    int i;
    numberV = dV;
    E = 0;
    EHead = new Edge();
    ETail = new Edge();
    EHead->next = ETail;
    ETail->back = EHead;
    ELhead = new EList *[numberV];
    ELtail = new EList *[numberV];
    for(i=0;i<numberV;i++){
      ELhead[i] = new EList();
      ELtail[i] = new EList();
      ELhead[i]->next = ELtail[i];
      ELtail[i]->back = ELhead[i];
    }
    
  }
  virtual ~Graph(){
    CleanEList();
    CleanEdge();
  }
 
  int CheckIsThereEdge(int dv1,int dv2){
    EList *now = ELhead[dv1];
    while(next(now)!=ELtail[dv1]){
      now = next(now);
      if(now->ID==dv2)return now->e->ID;
    }
    return -1;
  }



  
  void CleanEList(){
    EList *now=NULL;
    EList *dummy=NULL;
    int i;
    for(i=0;i<numberV;i++){
    
      now = ELhead[i]->next;
      if(now!=ELtail[i])
	while(now->next!=ELtail[i]){
	  dummy = now->next;
	  
	  delete now;
	  
	  now = dummy;
	}
    
      delete ELhead[i];
      delete ELtail[i]; 
      
    }
    delete [] ELhead;
    delete [] ELtail;
    
    
  }
  int AppendE(int dv1,int dv2){
    
    Edge *now = new Edge(E,dv1,dv2);
    Edge *dummy = ETail->back;
    now->next = ETail;
    ETail->back = now;
    now->back = dummy;
    dummy->next =now;
    AppendEL(now,dv1,dv2);
    AppendEL(now,dv2,dv1);
    E++;
    return now->ID;
  }
  void AppendEL(Edge *de,int dv1,int dID){
    EList *now = new EList(dID,de);
    EList *dummy = ELtail[dv1]->back;
    now->next = ELtail[dv1];
    ELtail[dv1]->back = now;
    now->back = dummy;
    dummy->next =now;
  }
  
  
  void CleanEdge(){
    Edge *dummy; 
    Edge *now;
    int i;
    now = EHead->next;
    if(now!=ETail)
    while(now->next!=ETail){
      dummy = now->next;
      delete now;
      now = dummy;
    }
    delete EHead;
    delete ETail;
    
  }
 private:
  Graph(const Graph& rhs);
  const Graph &operator=(const Graph& rhs);
  
};
class Point3d;


class Features{
 public:
  int V;
  int E;
  int NV;
  int connectN;
  int *connectID;
  //int *Pindex;
  //int *RPindex;
  

  double *connectVal;
  
     double *shapeval;
     double *cylinderval;
     
     double *symetryval;
  Edge *head;
  Edge *tail;
  IDList **NEhead;
  IDList **NEtail;
  int *neighborE;
  IDSet *mylist;
  Edge *pointer;
  Edge **pointerCheck;
  class Point{
  public:
    int ID;
    int ov1,ov2;
    double val1,val2;
    Point *next;
    Point *back;
    double x,y,z,kval,si,sy;
    
    int fID;
    int tconner;

    Point(){
      ID=-1;
      si=0.0;
      sy = 0.0;
      next = NULL;
      back = NULL;
      
      x=0.0;
      y=0.0;
      z=0.0;
      kval=0.0;
      val1=0.0;
      val2=0.0;
      ov1=-1;
      fID=-1;
      ov2=-1;tconner=0;
    }
    Point(int dID,double dx,double dy,double dz,int dfID){
      ID=dID;
      fID = dfID;
      x=dx;
      y=dy;
      z=dz;
      tconner=1;
    }
    Point(int dID,double dx,double dy,double dz,int dov1,int dov2,double dval1,double dval2){
      ID=dID;
      x=dx;
      y=dy;
      z=dz;
      ov1=dov1;
      ov2=dov2;
      val1 = dval1;
      val2 = dval2;tconner=0;
    }
    
    virtual ~Point(){
      
      
    }
    void setPoint3D(Point3d *out){
      out->x = x;
      out->y = y;
      out->z = z;
      
    }
    
  private:
    Point(const Point& rhs);
    const Point &operator=(const Point& rhs);
  };
  Point *phead;
  Point *ptail;
  Point **pBox;
  int NT;
  Features(){
    connectN=0;pointer=NULL;
    connectID=NULL;
    connectVal=NULL;
     shapeval = NULL;
    cylinderval = NULL;
    
    symetryval  = NULL;
    //Pindex=NULL;
    //RPindex=NULL;
    pointerCheck=NULL;
    mylist = new IDSet();
    V = 0;pBox=NULL;
    E=0;
    NT=0;
    NV=0;
    head = new Edge();
    tail = new Edge();
    head->next = tail;
    tail->back = head;
    phead = new Point();
    ptail = new Point();
    phead->next = ptail;
    ptail->back = phead;
    NEhead=NULL;
    NEtail=NULL;
    neighborE=NULL;
  }
  
  virtual ~Features(){
    //if(Pindex!=NULL)delete [] Pindex;
    //if(RPindex!=NULL)delete [] RPindex;
    
    if(connectID!=NULL)delete [] connectID;
    if(connectVal!=NULL)delete [] connectVal;
    if(shapeval!=NULL)delete [] shapeval;
    if(cylinderval!=NULL)delete [] cylinderval;
    
    CleanEdge();
    CleanPoint();
    
    mylist->CleanNeighborLL(NEhead,NEtail,(V-NT),neighborE);
    delete [] pBox;
    
    delete mylist;
    
  }
  void setNeighborIDs(int ID,int *dID1,int *dID2){
    if(neighborE[ID]==1){
      (*dID1) = NEhead[ID]->next->ID;
      (*dID2) = -1;
    }else{
    (*dID1) = NEhead[ID]->next->ID;
    (*dID2) = NEhead[ID]->next->next->ID;
        
    }
  
  }
  void setEdgePointer(){
    pointer = head;
  
  }
  void setEdgePointerNext(){
    pointer = next(head);
    
  }
  int setEdgeIDs(int *vID1,int *vID2,int *edgeFID,int *orivID1,int *orivID2){
    pointer = next(pointer);
    if(pointer==tail)return 0;
    (*vID1) = pointer->v1;
    (*vID2) = pointer->v2;
    (*edgeFID) = pointer->fID;
    (*orivID1) = pointer->ff1;
    (*orivID2) = pointer->ff2;
    
    return 1;
  
  }
  int setEdgeIDsNext(int *vID1,int *vID2,int *edgeFID,int *orivID1,int *orivID2){
    (*vID1) = pointer->v1;
    (*vID2) = pointer->v2;
    (*edgeFID) = pointer->fID;
    (*orivID1) = pointer->ff1;
    (*orivID2) = pointer->ff2;
    
    pointer = next(pointer);
    if(pointer==tail)return 0;
        
    return 1;
  
  }
  void startLabel(int dID,int label){
    if(connectID[dID]==-1){
      IDList *msHead = new IDList();
      IDList *msTail = new IDList();
      msHead->next = msTail;
      msTail->back = msHead;
      mylist->AppendVF(dID,msTail);
      
      
      connectID[dID]=label;
      while(msHead->next!=msTail){
	IDList *now = NEhead[msHead->next->ID];
	while(next(now)!=NEtail[msHead->next->ID]){
	  now = next(now);
	  if(connectID[now->ID]==-1){
	    connectID[now->ID]=label;
	    mylist->AppendVF(now->ID,msTail);
	  }
	}
	now = msHead->next; 
	IDList *nextI = now->next; 
	msHead->next = nextI;
	nextI->back = msHead;
	delete now;
      }
      delete msTail;
      delete msHead;
      
    }
  }
  int CleanUp(double ep){
    
    int i,j,dv1,dv2;
    int checkep=0;
    
    Edge *now = head;
    while(next(now)!=tail){
      now = next(now);
      if(Distance(now->v1,now->v2)<=ep){
	dv1 = now->v1;
	dv2 = now->v2;
	checkep=1;
	Edge *dback = now->back;
	Edge *dnext = now->next;
	dback->next = dnext;
	dnext->back = dback;
	delete now;
	now = dback;
	E--;
	Edge *now2 = head;
	while(next(now2)!=tail){
	  now2 = next(now2);
	  if(now2->v1==dv1)now2->v1 = dv2;
	  if(now2->v2==dv1)now2->v2 = dv2;
	}
      }
    }
    if(checkep==0)return 0;
    int *checkv = new int[V];
    for(i=0;i<V;i++)checkv[i]=-1;
   now = head;
    while(next(now)!=tail){
      now = next(now);
      checkv[now->v1]=1;
      checkv[now->v2]=1;
    }
    NV = 0;
    Point *nowp = phead;
    while(next(nowp)!=ptail){
      nowp = next(nowp);
      if(checkv[nowp->ID]!=1){
	Point *dback = nowp->back;
	Point *dnext = nowp->next;
	dback->next = dnext;
	dnext->back = dback;
	delete nowp;
	nowp = dback;
	
	NV++;
	if(next(nowp)==ptail)break;
      }
    }
    delete [] checkv;
    
    
    //printf("V = %d E = %d\n",V,E);
    
    return 1;
    
    
  }
  void Remove1Seg(){
    int i;
    int checkep=0;
    Edge *now = head;
    while(next(now)!=tail){
      now = next(now);
      if(neighborE[now->v1]==1&&neighborE[now->v2]==1){
	checkep=1;
	Edge *dback = now->back;
	Edge *dnext = now->next;
	dback->next = dnext;
	dnext->back = dback;
	delete now;
	now = dback;
	E--;
	
      }
    }
    
    if(checkep==0)return;
    int *checkv = new int[V];
    for(i=0;i<V;i++)checkv[i]=-1;
    now = head;
    while(next(now)!=tail){
      now = next(now);
      checkv[now->v1]=1;
      checkv[now->v2]=1;
    }
    NV = 0;
    Point *nowp = phead;
    while(next(nowp)!=ptail){
      nowp = next(nowp);
      if(checkv[nowp->ID]!=1){
	Point *dback = nowp->back;
	Point *dnext = nowp->next;
	dback->next = dnext;
	dnext->back = dback;
	delete nowp;
	nowp = dback;
	
	NV++;
	if(next(nowp)==ptail)break;
      }
    }
    delete [] checkv;
    makeEdgesNIDchange();
    
  
  }
  int CleanUp1Seg(double ep){
    
    int i,j,dv1,dv2;
    int checkep=0;
    Edge *now = head;
    while(next(now)!=tail){
      now = next(now);
      if(neighborE[now->v1]==1&&neighborE[now->v2]==1){
	checkep=1;
	Edge *dback = now->back;
	Edge *dnext = now->next;
	dback->next = dnext;
	dnext->back = dback;
	delete now;
	now = dback;
	E--;
	
      }else{
	if(Distance(now->v1,now->v2)<=ep){
	  dv1 = now->v1;
	  dv2 = now->v2;
	  checkep=1;
	  Edge *dback = now->back;
	  Edge *dnext = now->next;
	  dback->next = dnext;
	  dnext->back = dback;
	  delete now;
	  now = dback;
	  E--;
	  Edge *now2 = head;
	  while(next(now2)!=tail){
	    now2 = next(now2);
	    if(now2->v1==dv1)now2->v1 = dv2;
	  }
	}
      }
    }
    if(checkep==0)return 0;
    int *checkv = new int[V];
    for(i=0;i<V;i++)checkv[i]=-1;
    now = head;
    while(next(now)!=tail){
      now = next(now);
      checkv[now->v1]=1;
      checkv[now->v2]=1;
    }
    NV = 0;
    Point *nowp = phead;
    while(next(nowp)!=ptail){
      nowp = next(nowp);
      if(checkv[nowp->ID]!=1){
	Point *dback = nowp->back;
	Point *dnext = nowp->next;
	dback->next = dnext;
	dnext->back = dback;
	delete nowp;
	nowp = dback;
	
	NV++;
	if(next(nowp)==ptail)break;
      }
    }
    
    
    delete [] checkv;
    //printf("V = %d E = %d\n",V,E);
    
    return 1;
    
    
  }


  void Labeling(){
    if(V==0)return;
    connectN=0;
    if(connectID!=NULL)delete [] connectID;
    connectID = new int[V];
    int i;
    for(i=0;i<V;i++)connectID[i] = -1;
    for(i=0;i<V;i++){
      if(connectID[i]==-1){
	startLabel(i,connectN);
	connectN++;
      }
    }
  }  
  void setRidgeness(double *dk1,double *si,double *sy,int **Face){
    if(V==0||connectN<=1)return;
    int i;
    connectVal = new double[connectN];
    shapeval = new double[connectN];
    cylinderval = new double[connectN];
   
     
     int *coN = new int[connectN];
     double *coD = new double[connectN];
     for(i=0;i<connectN;i++){
       connectVal[i]=0.0;
       shapeval[i] = 0.0;coD[i]=0.0;
       cylinderval[i] = 0.0;coN[i]=0;
     }
    Point *nowp;
    nowp = phead;
    double sy1val;
    while(next(nowp)!=ptail){
      nowp = next(nowp);
      if(nowp->tconner==0){
	nowp->si = (nowp->val1*si[nowp->ov2]+nowp->val2*si[nowp->ov1])/(nowp->val1+nowp->val2);
	

	nowp->sy = (nowp->val1*sy[nowp->ov2]+nowp->val2*sy[nowp->ov1])/(nowp->val1+nowp->val2);
	

	nowp->kval = fabs(nowp->val1*dk1[nowp->ov2]+nowp->val2*dk1[nowp->ov1])/(nowp->val1+nowp->val2);
	
      }else{
	nowp->si = (si[Face[nowp->fID][0]]+si[Face[nowp->fID][1]]+si[Face[nowp->fID][2]])/3.0;
	
	nowp->sy = (sy[Face[nowp->fID][0]]+sy[Face[nowp->fID][1]]+sy[Face[nowp->fID][2]])/3.0;
	
	
	nowp->kval = fabs(dk1[Face[nowp->fID][0]]+dk1[Face[nowp->fID][1]]+dk1[Face[nowp->fID][2]])/3.0;
      }    
    }
    
   
    Edge *now = head;
    double dlenval=0.0;
    while(next(now)!=tail){
      now = next(now);
      dlenval = Distance(now->v1,now->v2);
      
      connectVal[connectID[now->v1]]+= 0.5*(pBox[now->v1]->kval+pBox[now->v2]->kval)*dlenval;
      
      
      shapeval[connectID[now->v1]] += 0.5*(pBox[now->v1]->si+pBox[now->v2]->si);
      
      cylinderval[connectID[now->v1]] += 0.5*(pBox[now->v1]->sy+pBox[now->v2]->sy)*dlenval;
      coD[connectID[now->v1]] += dlenval;
      coN[connectID[now->v1]]++;
    }
    
 
    for(i=0;i<connectN;i++){
      shapeval[i] /= ((double)(coN[i]));
 
      
      cylinderval[i] *= (coD[i]);
    }
        
    delete [] coN;
    delete [] coD;
  }
  
  void removeRiges(double TH,double TH1,double TH2){
    if(V==0||connectN<=1)return;
    int i,j;
    int *checkremove = new int[connectN];
    
    for(i=0;i<connectN;i++){
      checkremove[i]=0;
      if(!(connectVal[i]>=TH&&shapeval[i]>=TH1&&cylinderval[i]>=TH2)){
	
	checkremove[i]=1;
      }
      
    }
    
    //printf("V = %d E = %d\n",V,E);
    Point *nowp = phead;
    while(next(nowp)!=ptail){
      nowp = next(nowp);
      if(checkremove[connectID[nowp->ID]]==1){
	Point *dback = nowp->back;
	Point *dnext = nowp->next;
	dback->next = dnext;
	dnext->back = dback;
	delete nowp;
	nowp = dback;
	
	NV++;
	if(next(nowp)==ptail)break;
      }
    }
    Edge *now = head;
    while(next(now)!=tail){
      now = next(now);
      if(checkremove[connectID[now->v1]]==1){
	Edge *dback = now->back;
	Edge *dnext = now->next;
	dback->next = dnext;
	dnext->back = dback;
	delete now;
	now = dback;
	E--;
      }
    }//printf("V = %d E = %d\n",V,E);
  }
  double Distance(int v1,int v2){
    double dx = pBox[v1]->x - pBox[v2]->x;
    double dy = pBox[v1]->y - pBox[v2]->y;
    double dz = pBox[v1]->z - pBox[v2]->z;
    return sqrt(dx*dx+dy*dy+dz*dz);
  
  }
  void setcheckFace(int numberF){
    pointerCheck = new Edge* [numberF];
    int i;
    for(i=0;i<numberF;i++)pointerCheck[i] = NULL;
    Edge *now = head;
    while(next(now)!=tail){
      now = next(now);
      if(now->fID!=-1){
	pointerCheck[now->fID] = now;
      }
    }
    
  }
  void RemoveEdge(int fID){
    E--;
    Edge *now = pointerCheck[fID];
    Edge *next = now->next;
    Edge *back = now->back;
    next->back = back;
    back->next = next;
    delete now;
    pointerCheck[fID]=NULL;
  }
  void RemoveV(){
    int i;
    int *checkv = new int[V];
    for(i=0;i<V;i++)checkv[i]=-1;
    Edge *now = head;
    while(next(now)!=tail){
      now = next(now);
      checkv[now->v1]=1;
      checkv[now->v2]=1;
    }
    NV = 0;
    Point *nowp = phead;
    while(next(nowp)!=ptail){
      nowp = next(nowp);
      if(checkv[nowp->ID]!=1){
	Point *dback = nowp->back;
	Point *dnext = nowp->next;
	dback->next = dnext;
	dnext->back = dback;
	delete nowp;
	nowp = dback;
	NV++;
	if(next(nowp)==ptail)break;
      }
    }
        
    delete [] checkv; 
  }
  int isDiplicateEdge(int dOV1,int dOV2,int dv1,int dv2,double ep){
    double dis1 = Distance(dOV1,dv1);
    double dis2 = Distance(dOV1,dv2);
    double dis3 = Distance(dOV2,dv1);
    double dis4 = Distance(dOV2,dv2);
    printf("%lf %lf %lf %lf\n",dis1,dis2,dis3,dis4);
    if((dis1<ep&&dis3<ep)||
       (dis1<ep&&dis4<ep)||
       (dis2<ep&&dis3<ep)||
      (dis2<ep&&dis4<ep)){
      return 1;
    }else{
      return 0;
    }
  
  
  }

  int isDiplicateEdge(int df1,int df2,double ep){
    int dOV1 = pointerCheck[df1]->v1;
    int dOV2 = pointerCheck[df1]->v2;
    int dv1 = pointerCheck[df2]->v1;
    int dv2 = pointerCheck[df2]->v2;
    
    double dis1 = Distance(dOV1,dv1);
    double dis2 = Distance(dOV1,dv2);
    double dis3 = Distance(dOV2,dv1);
    double dis4 = Distance(dOV2,dv2);
    printf("%lf %lf %lf %lf\n",dis1,dis2,dis3,dis4);
    if((dis1<ep&&dis3<ep)||
       (dis1<ep&&dis4<ep)||
       (dis2<ep&&dis3<ep)||
      (dis2<ep&&dis4<ep)){
      return 1;
    }else{
      return 0;
    }
  
  
  }
  int EdgeFeature(int dOV1,int dOV2,double ep,int *orivID1,int *orivID2){
    double alpha1 = pBox[dOV1]->val1/(pBox[dOV1]->val1+pBox[dOV1]->val2);
    double beta1 = 1.0-alpha1;
    double alpha2 = pBox[dOV2]->val1/(pBox[dOV2]->val1+pBox[dOV2]->val2);
    double beta2 = 1.0-alpha2;
    
    if((alpha1<=ep&&alpha2<=ep)||
       (alpha1<=ep&&beta2<=ep)||
       (beta1<=ep&&alpha2<=ep)||
       (beta1<=ep&&beta2<=ep)){
      
      if(alpha1<=ep){
	(*orivID1) = pBox[dOV1]->ov1;
      }else{
	(*orivID1) = pBox[dOV1]->ov2;
      }
      if(alpha2<=ep){
	(*orivID2) = pBox[dOV2]->ov1;
      }else{
	(*orivID2) = pBox[dOV2]->ov2;
      }
      
      
      return 1;
    }else{
      return 0;
    }
  
  }
  void DeletecheckFace(){
    delete [] pointerCheck;
    pointerCheck = NULL;

  }
  void setVec(Point3d *out,int pID){
    //if(NEhead[pID]->next==NEtail[pID])printf("Error in setVec !!\n");
    out->x = pBox[pID]->x-pBox[NEhead[pID]->next->ID]->x;
    out->y = pBox[pID]->y-pBox[NEhead[pID]->next->ID]->y;
    out->z = pBox[pID]->z-pBox[NEhead[pID]->next->ID]->z;
    
  
  }
  
  int Connect(int p1,int p2,int o1,int o2){
    if(NEhead[p1]->next->ID==p2||
       NEhead[p2]->next->ID==p1)return -1;
    if(NEhead[p1]->next->ID==NEhead[p2]->next->ID)return -1;
    
      
      AppendE(p1,p2,o1,o2);
      
    //mylist->AppendVF(p1,NEtail[p2]);neighborE[p2]++;
    //mylist->AppendVF(p2,NEtail[p1]);neighborE[p1]++;
      return 0;
  }
  void makeEdgesN(){
 
    if(V!=0){
      if(NEhead!=NULL){
	delete [] pBox;
	mylist->CleanNeighborLL(NEhead,NEtail,(V-NT),neighborE);
	NT=0;
	V = V-NV;
	NV=0;
      }else{
	
      }
      pBox = new Point* [V];
      NEhead = new IDList *[V];
      NEtail = new IDList *[V];
      neighborE = new int[V];
      int i;
      for(i=0;i<V;i++){
	NEhead[i] = new IDList();
	NEtail[i] = new IDList();
	NEhead[i]->next = NEtail[i];
	NEtail[i]->back = NEhead[i];
	neighborE[i]=0;
      }
      
      Point *snow = phead;i=0;
      while(next(snow)!=ptail){
	snow= next(snow);
	//snow->ID=i;
	pBox[snow->ID]=snow;
	//i++;
      }
      
      
      Edge *now = head;
      while(next(now)!=tail){
	now = next(now);
	//now->v1 = now->v1;
	//now->v2 = now->v2;
	
	mylist->AppendVF(now->v1,NEtail[now->v2]);neighborE[now->v2]++;
	mylist->AppendVF(now->v2,NEtail[now->v1]);neighborE[now->v1]++;
      }
      
    }
    
  }
  void makeEdgesNIDchange(){
    int *index;
    if(V!=0){
      if(NEhead!=NULL){
	delete [] pBox;
	mylist->CleanNeighborLL(NEhead,NEtail,(V-NT),neighborE);
	NT=0;
	index = new int[V];
	V = V-NV;
	NV=0;
      }else{
	index = new int[V];
      }
      
      
      pBox = new Point* [V];
      NEhead = new IDList *[V];
      NEtail = new IDList *[V];
      neighborE = new int[V];
      
      int i;
      for(i=0;i<V;i++){
	NEhead[i] = new IDList();
	NEtail[i] = new IDList();
	NEhead[i]->next = NEtail[i];
	NEtail[i]->back = NEhead[i];
	neighborE[i]=0;
      }
      
      
      Point *snow = phead;i=0;
      while(next(snow)!=ptail){
	snow= next(snow);
	index[snow->ID]=i;
	snow->ID=i;
	pBox[i]=snow;
	i++;
      }
      
 
      Edge *now = head;
      while(next(now)!=tail){
	now = next(now);
	now->v1 = index[now->v1];
	now->v2 = index[now->v2];
	
	mylist->AppendVF(now->v1,NEtail[now->v2]);neighborE[now->v2]++;
	mylist->AppendVF(now->v2,NEtail[now->v1]);neighborE[now->v1]++;
      }
      
      

      delete [] index;
     
      
      
    }
    
  }
  
  
  int IstEdge(int dv1,int dv2){
    IDList *now = NEhead[dv1];
    while(next(now)!=NEtail[dv1]){
      now = next(now);
      if(now->ID==dv2){
	return 1;
      }
    }
    return 0;
  }
  void setCandidateVector(Point3d *out,int dv1,int dv2){
    out->x = pBox[dv2]->x - pBox[dv1]->x;
    out->y = pBox[dv2]->y - pBox[dv1]->y;
    out->z = pBox[dv2]->z - pBox[dv1]->z;
  }
  void PushT(int dv1,int dv2,int dv3,double dx,double dy,double dz,int fID){
   
    int newpid = AppendP(V,dx,dy,dz,fID);
    if(IstEdge(dv1,dv2)==1)removeEdge(dv1,dv2);
    if(IstEdge(dv2,dv3)==1)removeEdge(dv2,dv3);
    if(IstEdge(dv3,dv1)==1)removeEdge(dv3,dv1);
    AppendE(newpid,dv1,fID);mylist->AppendVF(newpid,NEtail[dv1]);neighborE[dv1]++;
    AppendE(newpid,dv2,fID);mylist->AppendVF(newpid,NEtail[dv2]);neighborE[dv2]++;
    AppendE(newpid,dv3,fID);mylist->AppendVF(newpid,NEtail[dv3]);neighborE[dv3]++;
    // now conner neighbor edges are not assigned 
    NT++;
  }
  void removeEdge(int dv1,int dv2){
    Edge *now = head;
    Edge *dback;
    Edge *dnext;
    while(next(now)!=tail){
      now = next(now);
      if((now->v1==dv1&&now->v2==dv2)||
	(now->v2==dv1&&now->v1==dv2)){
	dback = now->back;
	dnext = now->next;
	dback->next = dnext;
	dnext->back = dback;
	delete now;
	break;
      }
    }
    IDList *inow = NEhead[dv1];
    IDList *iback;
    IDList *inext;
    while(next(inow)!=NEtail[dv1]){
      inow = next(inow);
      if(inow->ID==dv2){
	iback = inow->back;
	inext = inow->next;
	iback->next = inext;
	inext->back = iback;
	neighborE[dv1]--;
	delete inow;
	break;
      }
    }
    inow = NEhead[dv2];
    while(next(inow)!=NEtail[dv2]){
      inow = next(inow);
      if(inow->ID==dv1){
	iback = inow->back;
	inext = inow->next;
	iback->next = inext;
	inext->back = iback;
	neighborE[dv2]--;
	delete inow;
	break;
      }
    }
    E--;
  
  }
  void Connect1RingEdge(double Phi,double Theta,int numberV,Point3d **point,Point3d **bc,int *check1v,PointTool *PT,IDList **IHead,IDList **ITail){
    int i,cv1,cv2,checktj;
    if(V==0)return;
    int target=0;
    double anglec=0.0;
    double dan=0.0;
    double anglec2=0.0;
    double anglec3=0.0;
    for(i=0;i<numberV;i++)check1v[i]=-1;
    for(i=0;i<V;i++){
      if(neighborE[i]==1){
	
	cv1 = pBox[i]->ov1;
	cv2 = pBox[i]->ov2;
	bc[0]->x = point[cv1]->x-pBox[i]->x;
	bc[0]->y = point[cv1]->y-pBox[i]->y;
	bc[0]->z = point[cv1]->z-pBox[i]->z;
	bc[1]->x = point[cv2]->x-pBox[i]->x;
	bc[1]->y = point[cv2]->y-pBox[i]->y;
	bc[1]->z = point[cv2]->z-pBox[i]->z;
	 if(PT->Point3dSize(bc[0])<PT->Point3dSize(bc[1])){
	   check1v[cv1]=i;
	   
	 }else{
	   check1v[cv2]=i;
	   
	 }
	 
       }
    }
    IDList *now;
    for(i=0;i<numberV;i++){
      if(check1v[i]!=-1){
	now = IHead[i];
	checktj=0;
	setVec(bc[0],check1v[i]);
	PT->Normalize3D(bc[0]);
	while(next(now)!=ITail[i]){
	  now = next(now);
	  if(check1v[now->ID]!=-1){
	    if(checktj==0){
	      target=now->ID;
	      checktj=1;
	      setVec(bc[1],check1v[now->ID]);
	      PT->Normalize3D(bc[1]);
	      anglec = (PT->InnerProduct(bc[0],bc[1]));
	    }else{
	      setVec(bc[1],check1v[now->ID]);
	      PT->Normalize3D(bc[1]);
	      dan = (PT->InnerProduct(bc[0],bc[1]));
	      if(anglec>dan){
		anglec = dan;
		target=now->ID;
	      }
	    }
	    
	  }
	  
	}
	if(checktj==1){
	  setVec(bc[1],check1v[target]);PT->Normalize3D(bc[1]);
	  setCandidateVector(bc[2],check1v[i],check1v[target]);PT->Normalize3D(bc[2]);
	  setCandidateVector(bc[3],check1v[target],check1v[i]);PT->Normalize3D(bc[3]);
	  
	  anglec2 = acos(PT->InnerProduct(bc[1],bc[3]));
	  anglec3 = acos(PT->InnerProduct(bc[0],bc[2]));
	  
	  
	  if(acos(anglec)>=Phi&&anglec3<=Theta&&anglec2<=Theta){
	    
	    Connect(check1v[i],check1v[target],i,target);
	    check1v[i]=-1;
	  }
	}
      }
            
    }
        
    makeEdgesNIDchange();
    
  }
  void Connect1RingPoint(int numberV,Point3d **point,Point3d **bc,int *check1v,int *check2v,PointTool *PT){
    int i,cv1,cv2;double d1,d2;
    if(V==0)return;
    double *dis1 = new double[numberV];
    double *dis2 = new double[numberV];
    
    for(i=0;i<numberV;i++){
      check1v[i]=-1;
      check2v[i]=-1;
      dis1[i]=0.0;
      dis2[i]=0.0;
    }
    for(i=0;i<V;i++){
      if(neighborE[i]==1){
	
	cv1 = pBox[i]->ov1;
	cv2 = pBox[i]->ov2;
	bc[0]->x = point[cv1]->x-pBox[i]->x;
	bc[0]->y = point[cv1]->y-pBox[i]->y;
	bc[0]->z = point[cv1]->z-pBox[i]->z;
	bc[1]->x = point[cv2]->x-pBox[i]->x;
	bc[1]->y = point[cv2]->y-pBox[i]->y;
	bc[1]->z = point[cv2]->z-pBox[i]->z;
	d1 = PT->Point3dSize(bc[0]);
	d2 = PT->Point3dSize(bc[0]);
	
	if(d1<d2){
	  if(check1v[cv1]==-1){
	    check1v[cv1]=i;
	    dis1[cv1]=d1;
	  }else{
	    if(check2v[cv1]==-1){
	      check2v[cv1]=i;
	      dis2[cv1]=d1;
	    }else{
	      if(d1<check1v[cv1]&&check1v[cv1]<check2v[cv1]){
		check2v[cv1]=i;
		dis2[cv1]=d1;
		
	      }else if(d1<check2v[cv1]&&check2v[cv1]<check1v[cv1]){
		check1v[cv1]=i;
		dis1[cv1]=d1;
		
	      }
	      
	    }
	  }
	}else{
	  if(check1v[cv2]==-1){
	    check1v[cv2]=i;
	    dis1[cv2]=d2;
	  }else{
	    if(check2v[cv2]==-1){
	      check2v[cv2]=i;
	      dis2[cv2]=d2;
	    }else{
	      if(d2<check1v[cv2]&&check1v[cv2]<check2v[cv2]){
		check1v[cv2]=i;
		dis1[cv2]=d2;
		
	      }else if(d2<check2v[cv2]&&check2v[cv2]<check1v[cv2]){
		check2v[cv2]=i;
		dis2[cv2]=d2;
		
	      }
	      
	    }
	  }
	}
	
      }
    }
    int check = 0;
    for(i=0;i<numberV;i++){
      if(check1v[i]!=-1&&check2v[i]!=-1){
	//printf("%d %d\n",check1v[i],check2v[i]);
	check = 1;
	Connect(check1v[i],check2v[i],i,i);
      }
    }
    if(check == 1){
      makeEdgesNIDchange();
    }
    delete [] dis1;
    delete [] dis2;
    
  }


  double checkBoundaydist(Point3d *out,int pID){
    
    out->x = pBox[pID]->x-pBox[NEhead[pID]->next->ID]->x;
    out->y = pBox[pID]->y-pBox[NEhead[pID]->next->ID]->y;
    out->z = pBox[pID]->z-pBox[NEhead[pID]->next->ID]->z;
    return sqrt((out->x*out->x+out->y*out->y+out->z*out->z));
  }
  double checkBoundaydist(Point3d *out,int pID,int backID){
    
    if(NEhead[pID]->next->ID==backID){
      out->x = pBox[pID]->x-pBox[NEhead[pID]->next->next->ID]->x;
      out->y = pBox[pID]->y-pBox[NEhead[pID]->next->next->ID]->y;
      out->z = pBox[pID]->z-pBox[NEhead[pID]->next->next->ID]->z;
    }else{
      out->x = pBox[pID]->x-pBox[NEhead[pID]->next->ID]->x;
      out->y = pBox[pID]->y-pBox[NEhead[pID]->next->ID]->y;
      out->z = pBox[pID]->z-pBox[NEhead[pID]->next->ID]->z;
    }
    return sqrt((out->x*out->x+out->y*out->y+out->z*out->z));
  }
  
  
  void RemoveDiplicatesInner(double ep){
    if(CleanUp(ep)==1){
      makeEdgesNIDchange(); 
    }
  }
  void RemoveDiplicatesInnerWith1Seg(double ep){
    if(CleanUp1Seg(ep)==1){
      makeEdgesNIDchange(); 
    }
  }
  
  void RemoveDiplicatesBoundary(double ep,Point3d **bc){
    int checktj=0;
    int i;
    NV=0;
    double distd=0.0;
    int dnextid=-1;
    int dnextidtmp=-1;
    int dbackID=-1;
    for(i=0;i<V;i++){
      if(neighborE[i]==1){
	distd = checkBoundaydist(bc[0],i);
	if(distd<ep){
	  checktj=1;
	  dnextid = RemoveBoundaryV(i,-1);
	  dbackID = i;
	  while(dnextid!=-1){
	    
	    distd = checkBoundaydist(bc[0],dnextid,dbackID);
	    if(distd<ep){
	      dnextidtmp = dnextid;
	      dnextid = RemoveBoundaryV(dnextidtmp,dbackID);
	      dbackID = dnextidtmp;
	    }else{
	      break;
	      
	    }
	    
	  }
	}
	
      }
    }
    
    
    if(checktj==1){
      RemoveBoundaryE();
      makeEdgesNIDchange();
    } 
    
    
    
  }
  int RemoveBoundaryV(int pID,int backID){
    NV++;
    Point *now = pBox[pID];
    
    Point *next = now->next;
    Point *back = now->back;
    next->back = back;
    back->next = next; 
    delete now;
    int nextID = NEhead[pID]->next->ID;
    if(nextID==backID){
      nextID = NEhead[pID]->next->next->ID;
    }
    
    if(neighborE[nextID]==1){
      NV++;
      neighborE[nextID]=0;
      now = pBox[nextID];
      
      next = now->next;
      back = now->back;
      next->back = back;
      back->next = next; 
      delete now;  
      
      pBox[nextID]=NULL;
      pBox[pID]=NULL;
      return -1;
    }
    pBox[pID]=NULL;
    return nextID;
    
    
  }
  void RemoveBoundaryE(){
    
    Edge *now = head;
    Edge *next = NULL;
    Edge *back = NULL;
    while(next(now)!=tail){
      now = next(now);
      
      if(pBox[now->v1]==NULL||pBox[now->v2]==NULL){
	
	E--;
	next = now->next;
	back = now->back;
	next->back = back;
	back->next = next;
	
	delete now;
	
	now = back;
      }
    }
  }
  void AppendE(int dv1,int dv2){
    
    
    Edge *now = new Edge(dv1,dv2);
    Edge *dummy = tail->back;
    now->next = tail;
    tail->back = now;
    now->back = dummy;
    dummy->next =now;
    E++;
   
    
  }
  void AppendE(int dv1,int dv2,int o1,int o2){
    
    
    Edge *now = new Edge(dv1,dv2);
    now->ff1 = o1;
    now->ff2 = o2;
    Edge *dummy = tail->back;
    now->next = tail;
    tail->back = now;
    now->back = dummy;
    dummy->next =now;
    E++;
   
    
  }
  
  void AppendE(int dv1,int dv2,int dfID){
    
    
    Edge *now = new Edge(dv1,dv2);
    now->fID = dfID;
    Edge *dummy = tail->back;
    now->next = tail;
    tail->back = now;
    now->back = dummy;
    dummy->next =now;
    E++;
   
    
  }
  void AppendEadd(int dv1,int dv2){
    
    
    Edge *now = new Edge(dv1,dv2);
    Edge *dummy = tail->back;
    now->next = tail;
    tail->back = now;
    now->back = dummy;
    dummy->next =now;
    E++;
    mylist->AppendVF(dv2,NEtail[dv1]);neighborE[dv1]++;
    mylist->AppendVF(dv1,NEtail[dv2]);neighborE[dv2]++;
    
  }
  
  int AppendP(int ID,double dx,double dy,double dz,int fID){
    
    Point *now = new Point(ID,dx,dy,dz,fID);
    Point *dummy = ptail->back;
    now->next = ptail;
    ptail->back = now;
    now->back = dummy;
    dummy->next =now;
    V++;
    return ID;
  }
  int AppendP(int ID,double dx,double dy,double dz,int ov1,int ov2,double dval1,double dval2){
    
    Point *now = new Point(ID,dx,dy,dz,ov1,ov2,dval1,dval2);
    Point *dummy = ptail->back;
    now->next = ptail;
    ptail->back = now;
    now->back = dummy;
    dummy->next =now;
    V++;
    return ID;
  }
  
  
  void CleanEdge(){
    Edge *dummy; 
    Edge *now;
    int i;
    now = head->next;
    if(now!=tail)
      while(now->next!=tail){
	dummy = now->next;
	delete now;
	now = dummy;
      }
    delete head;
    delete tail;
  }
  
  void PrintConnect(char* filename,int rcheck){
   
    FILE *out = fopen(filename,"w");
    if(rcheck==0){
      fprintf(out,"%d\n",0);
      fprintf(out,"%d\n",0);
      fprintf(out,"%d\n",0);
      fclose(out);
      return;
    }
    
    fprintf(out,"%d\n",V);
    fprintf(out,"%d\n",E);
    fprintf(out,"%d\n",connectN);
    
    Point *now = phead;
    int i=0;
    
    while(next(now)!=ptail){
      now= next(now);
      fprintf(out,"%lf %lf %lf %d\n",now->x,now->y,now->z,connectID[i]);
      i++;
    }
    
    
    for(i=0;i<connectN;i++)fprintf(out,"%lf %lf %lf\n",connectVal[i],shapeval[i],cylinderval[i]);
    
    Edge *nowe = head;
    while(next(nowe)!=tail){
      nowe = next(nowe);
      fprintf(out,"%d %d %d\n",nowe->v1,nowe->v2,nowe->fID);
    }
    
    fclose(out);
    
  }
  void PrintConnect(unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges,int rcheck){

	  if(rcheck==0){
		  return;
	  }

	  pNum=V;
	  eNum=E;
	  cNum=connectN;
	  fPoints=new double[V*3];
	  fPointType=new unsigned[V];
	  fFeature=new double[connectN*3];
	  fEdges=new unsigned[E*3];

	  Point *now = phead;
	  int i=0;

	  while(next(now)!=ptail){
		  now= next(now);
		 // fprintf(out,"%lf %lf %lf %d\n",now->x,now->y,now->z,connectID[i]);
		  fPoints[i*3+0]=now->x;fPoints[i*3+1]=now->y;fPoints[i*3+2]=now->z;
		  fPointType[i]=connectID[i];
		  i++;
	  }


	  for(i=0;i<connectN;i++){
		 // fprintf(out,"%lf %lf %lf\n",connectVal[i],shapeval[i],cylinderval[i]);
		  //fFeature[i*3+0]=connectVal[i];fFeature[i*3+1]=shapeval[i];fFeature[i*3+2]=cylinderval[i];
		  fFeature[i*3+0]=fFeature[i*3+1]=fFeature[i*3+2]=1;
	  }

	  Edge *nowe = head;
	  i=0;
	  while(next(nowe)!=tail){
		  nowe = next(nowe);
		  //fprintf(out,"%d %d %d\n",nowe->v1,nowe->v2,nowe->fID);
		  fEdges[i*3+0]=nowe->v1;fEdges[i*3+1]=nowe->v2;fEdges[i*3+2]=nowe->fID;
		  i++;
	  }

  }


  void CleanPoint(){
    Point *dummy; 
    Point *now;
    int i;
    now = phead->next;
    if(now!=ptail)
      while(now->next!=ptail){
	dummy = now->next;
	delete now;
	now = dummy;
      }
    delete phead;
    delete ptail;
  }
  
 private:
  Features(const Features& rhs);
  const Features &operator=(const Features& rhs);
};
