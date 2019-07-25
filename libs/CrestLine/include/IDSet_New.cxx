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
#include<stdio.h>
#include"IDList.h"
#include"PolarList.h"
#include"IDSet.h"

void IDSet::AppendPolarI(int dID,PolarList *dITail,double dx,double dy){
  PolarList *now = new PolarList(dID,dx,dy);
  PolarList *dummy = dITail->back;
  now->next = dITail;
  dITail->back = now;
  now->back = dummy;
  dummy->next =now;
}
void IDSet::AppendVF(int myID,IDList *dVTail)
    {
      IDList *now = new IDList(myID);
      
      IDList *dummy = dVTail->back;
      now->next = dVTail;
      dVTail->back = now;
      now->back = dummy;
      dummy->next =now;
    }
/*
void IDSet::AppendVF(int myID,CotList *dVTail)
    {
      CotList *now = new CotList(myID);
      
      CotList *dummy = dVTail->back;
      now->next = dVTail;
      dVTail->back = now;
      now->back = dummy;
      dummy->next =now;
    }
*/
int IDSet::SearchI(int dID,IDList *dIHead,IDList *dITail){
    
    IDList *now = dIHead;
    while(now->next!=dITail){
      now = now->next;
      if((now->ID==dID)){
	return 1;
      }
    }
    return 0;
}
 
  void IDSet::AppendI(int dID,IDList *dIHead,IDList *dITail,int nowID){
      if(dID!=nowID){
	if(SearchI(dID,dIHead,dITail)==0){
	  IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	  
	}
      }else{
	
      }
    }
 void IDSet::AppendI(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum)
    {
      if(dID!=nowID){
	if(SearchI(dID,dIHead,dITail)==0){
	  IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	  dnum[nowID]++;
	}
      }else{
	
      }
    }
void IDSet::AppendI(int dID,IDList *dIHead,IDList *dITail){
      
	if(SearchI(dID,dIHead,dITail)==0){
	  IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	}
 
}
void IDSet::AppendISort(int dID,IDList *dIHead,IDList *dITail,int *num){
  
  
    
    if(SearchI(dID,dIHead,dITail)==0){
      
      IDList *now = new IDList(dID);
      IDList *dummy = dITail->back;
      if(dIHead->next!=dITail)
	while(dID<dummy->ID){
	  dummy = dummy->back;
	  if(dummy==dIHead)break;
	}
      dummy = dummy->next;
      IDList *dummyback=dummy->back;
      
      now->next = dummy;
      dummy->back = now;
      now->back = dummyback;
      dummyback->next =now;
      (*num)++;
    }
    
  
  
  
}
void IDSet::AppendISort(int *dnum,int dID,IDList *dIHead,IDList *dITail,int nowID){
  if(dID!=nowID){
    if(SearchI(dID,dIHead,dITail)==0){
     
     IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	  
      (*dnum)++;
    }
  }else{
    
  }
  
}

void IDSet::AppendISort(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum){
  if(dID!=nowID){
    if(SearchI(dID,dIHead,dITail)==0){
      
      IDList *now = new IDList(dID);
      IDList *dummy = dITail->back;
      while(dID<dummy->ID){
	dummy = dummy->back;
	if(dummy==dIHead)break;
      }dummy = dummy->next;
      IDList *dummyback=dummy->back;
      
      now->next = dummy;
      dummy->back = now;
      now->back = dummyback;
      dummyback->next =now;
      dnum[nowID]++;
    }
  }else{
    
  }
  
}
void IDSet::AppendIF(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum)
    {
      
	if(SearchI(dID,dIHead,dITail)==0){
	  IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	  dnum[nowID]++;
	}
    }
void IDSet::AppendICl(int dID,IDList *dIHead,IDList *dITail,int *dnum)
    {
      
	if(SearchI(dID,dIHead,dITail)==0){
	  IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	  (*dnum)++;
	}
    }


void IDSet::CleanNeighbor(IDList* dHead,IDList* dTail){
  IDList *dummy; IDList *now;
  int i;
  now = dHead->next;
  if(now!=dTail){
    while(now->next!=dTail){
      dummy = now->next;
      delete now;
      now = dummy;
    }
    now = dTail->back;
    delete now;
  }
  delete dHead;dHead=NULL;
  delete dTail;dTail=NULL;
}
void IDSet::CleanNeighborNot(IDList* dHead,IDList* dTail){
  IDList *dummy; IDList *now;
  int i;
  now = dHead->next;
  if(now!=dTail){
    while(now->next!=dTail){
      dummy = now->next;
      delete now;
      now = dummy;
    }
    now = dTail->back;
    delete now;
  }
  
}

void IDSet::Clean(IDList **dFHead,IDList **dFTail,int numberSV,int *dneighborN){
  IDList *now=NULL;
  IDList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++){
    if(dneighborN[i]!=0){
      now = dFHead[i]->next;
      if(now!=dFTail[i]){
	while(now->next!=dFTail[i]){
	  dummy = now->next;
	  
	  delete now;
	  
	  now = dummy;
	}
	now = dFTail[i]->back;
	delete now;
      }
    }
    
    dFHead[i]->next = dFTail[i];
    dFTail[i]->back = dFHead[i];
    dneighborN[i]=0;
  }
  
}
  
void IDSet::CleanNeighborLL(IDList **dFHead,IDList **dFTail,int numberSV,int *dneighborN){
  IDList *now=NULL;
  IDList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++){
    if(dneighborN[i]!=0){
      now = dFHead[i]->next;
      if(now!=dFTail[i]){
	while(now->next!=dFTail[i]){
	  dummy = now->next;
	  
	  delete now;
	  
	  now = dummy;
	}
	now = dFTail[i]->back;
	delete now;
      }
    }
    delete dFHead[i];
    delete dFTail[i]; 
    
  }
  delete [] dneighborN;dneighborN=NULL;
  delete [] dFHead;
  delete [] dFTail;
  }
void IDSet::CleanNeighborL(IDList **dFHead,IDList **dFTail,int numberSV){
  IDList *now=NULL;
  IDList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++){
    
    now = dFHead[i]->next;
    if(now!=dFTail[i]){
      while(now->next!=dFTail[i]){
	dummy = now->next;
	
	delete now;
	
	now = dummy;
      }
      now = dFTail[i]->back;
      delete now;
    }
    
    delete dFHead[i];
    delete dFTail[i]; 
    
  }
  delete [] dFHead;
  delete [] dFTail;
  }
/*
void IDSet::CleanNeighborL(CotList **dFHead,CotList **dFTail,int numberSV){
  CotList *now=NULL;
  CotList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++){
    
    now = dFHead[i]->next;
    
    while(now->next!=dFTail[i]){
      dummy = now->next;
      
      delete now;
      
      now = dummy;
    }
    
    delete dFHead[i];
    delete dFTail[i]; 
    
  }
  delete [] dFHead;
  delete [] dFTail;
  }
*/
void IDSet::CleanNeighborL(PolarList **dFHead,PolarList **dFTail,int numberSV){
  PolarList *now=NULL;
  PolarList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++){
    
    now = dFHead[i]->next;
    if(now!=dFTail[i]){
      while(now->next!=dFTail[i]){
	dummy = now->next;
	
	delete now;
	
	now = dummy;
      }
      now = dFTail[i]->back;
      delete now;
    }
    
    delete dFHead[i];
    delete dFTail[i]; 
    
  }
  delete [] dFHead;
  delete [] dFTail;
  }




