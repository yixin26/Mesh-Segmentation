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

class IDList{
 public:
  int ID;
  IDList *next;
  IDList *back;
  IDList(){ID=-1;next=NULL;back=NULL;}
  IDList(int dv){ID=dv;next=NULL;back=NULL;}
  virtual ~IDList(){ID = -1; next=NULL;back=NULL;}
 private:
  IDList(const IDList& rhs);
  const IDList &operator=(const IDList& rhs);
   
};
