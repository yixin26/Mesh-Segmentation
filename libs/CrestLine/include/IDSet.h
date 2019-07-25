

class IDSet{
  
 public:
  IDSet(){}
  virtual ~IDSet(){}
  void AppendVF(int,IDList *);
  //void AppendVF(int,CotList *);
  
  int SearchI(int dID,IDList *dIHead,IDList *dITail);
  void AppendI(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum);
  void AppendI(int dID,IDList *dIHead,IDList *dITail,int nowID);
  
  void AppendIF(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum);
  void AppendISort(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum);
  void AppendISort(int *dnum,int dID,IDList *dIHead,IDList *dITail,int nowID);
  
  void AppendI(int dID,IDList *dIHead,IDList *dITail);
  void AppendISort(int dID,IDList *dIHead,IDList *dITail,int *num);
  void CleanNeighbor(IDList*,IDList*);
  
  void CleanNeighborLL(IDList **,IDList **,int ,int *);
  void CleanNeighborL(IDList **,IDList **,int );
  //void CleanNeighborL(CotList **,CotList **,int );
  void CleanNeighborL(PolarList **,PolarList **,int );
  void CleanNeighborNot(IDList* dHead,IDList* dTail);
  void Clean(IDList **dFHead,IDList **dFTail,int numberSV,int *dneighborN);
  void AppendICl(int dID,IDList *dIHead,IDList *dITail,int *dnum);
  void AppendPolarI(int dID,PolarList *dITail,double dx,double dy);
 private:
  IDSet(const IDSet& rhs);
  const IDSet &operator=(const IDSet& rhs);

};


