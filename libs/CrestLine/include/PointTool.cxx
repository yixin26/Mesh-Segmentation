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
#include<math.h>
#include"Point3d.h"
#include"PointTool.h"
void PointTool::setInverse(double out[3][3],double in[3][3]){
  double div = (in[0][2]*(-(in[1][1]*in[2][0]) + in[1][0]*in[2][1]) + in[0][1]*(in[1][2]*in[2][0] - in[1][0]*in[2][2]) + in[0][0]*(-(in[1][2]*in[2][1]) + in[1][1]*in[2][2]));if(div==0.0)div=1.0;
  out[0][0] = (-(in[1][2]*in[2][1]) + in[1][1]*in[2][2])/div;
  out[0][1] = (in[0][2]*in[2][1] - in[0][1]*in[2][2])/div;
  out[0][2] = (-(in[0][2]*in[1][1]) + in[0][1]*in[1][2])/div;
  out[1][0] = (in[1][2]*in[2][0] - in[1][0]*in[2][2])/div;
  out[1][1] =  (-(in[0][2]*in[2][0]) + in[0][0]*in[2][2]) /div;
  out[1][2] = (in[0][2]*in[1][0] - in[0][0]*in[1][2])/div;
  out[2][0] =  (-(in[1][1]*in[2][0]) + in[1][0]*in[2][1])/div;
  out[2][1] = (in[0][1]*in[2][0] - in[0][0]*in[2][1])/div;
  out[2][2] = (-(in[0][1]*in[1][0]) + in[0][0]*in[1][1])/div;
  

}
double PointTool::getTraceJTJ(double in[3][3]){
  return (in[0][0]*in[0][0]  + in[1][0]*in[1][0]  + in[2][0]*in[2][0]+
	  in[0][1]*in[0][1]  + in[1][1]*in[1][1]  + in[2][1]*in[2][1]+
	  in[0][2]*in[0][2]  + in[1][2]*in[1][2]  + in[2][2]*in[2][2]);


}
void PointTool::setJTJ(double out[3][3],double in[3][3]){
  out[0][0] = in[0][0]*in[0][0]  + in[1][0]*in[1][0]  + in[2][0]*in[2][0];
  out[0][1] = in[0][0]*in[0][1] + in[1][0]*in[1][1] + in[2][0]*in[2][1];
  out[0][2] = in[0][0]*in[0][2] + in[1][0]*in[1][2] + in[2][0]*in[2][2];
  out[1][0] = out[0][1];
  out[1][1] = in[0][1]*in[0][1]  + in[1][1]*in[1][1]  + in[2][1]*in[2][1];
  out[1][2] = in[0][1]*in[0][2] + in[1][1]*in[1][2] + in[2][1]*in[2][2];
  out[2][0] = out[0][2];
  out[2][1] = out[1][2];
  out[2][2] = in[0][2]*in[0][2]  + in[1][2]*in[1][2]  + in[2][2]*in[2][2];


}
void PointTool::setNormalF(Point3d* out,double dU,double dV,double *Func){
  out->x = Func[0]*dU+Func[1]*dV+3.0*Func[3]*dU*dU+2.0*Func[4]*dU*dV+Func[5]*dV*dV;
  out->y = Func[1]*dU+Func[2]*dV+Func[4]*dU*dU+2.0*Func[5]*dU*dV+3.0*Func[6]*dV*dV;
  
  
  out->z = -1.0;

  


}
double PointTool::DetCC(double **in){
  return (-in[0][2]*in[1][1]*in[2][0] + in[0][1]*in[1][2]*in[2][0] + in[0][2]*in[1][0]*in[2][1] - in[0][0]*in[1][2]*in[2][1] - in[0][1]*in[1][0]*in[2][2] + in[0][0]*in[1][1]*in[2][2] + in[0][2]*in[1][1]*in[3][0] - in[0][1]*in[1][2]*in[3][0] - in[0][2]*in[2][1]*in[3][0] + in[1][2]*in[2][1]*in[3][0] + in[0][1]*in[2][2]*in[3][0] - in[1][1]*in[2][2]*in[3][0] - in[0][2]*in[1][0]*in[3][1] + in[0][0]*in[1][2]*in[3][1] + in[0][2]*in[2][0]*in[3][1] - in[1][2]*in[2][0]*in[3][1] - in[0][0]*in[2][2]*in[3][1] + in[1][0]*in[2][2]*in[3][1] + in[0][1]*in[1][0]*in[3][2] - in[0][0]*in[1][1]*in[3][2] - in[0][1]*in[2][0]*in[3][2] + in[1][1]*in[2][0]*in[3][2] +  in[0][0]*in[2][1]*in[3][2] - in[1][0]*in[2][1]*in[3][2]);



}
void PointTool::setRodriguesM(double out[3][3], Point3d *dr, double costheta){
  int i,j;
  double ds = sqrt(fabs(1.0-costheta*costheta));
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)out[i][j] = 0.0;
  out[0][0] = costheta;
  out[1][1] = costheta;
  out[2][2] = costheta;
  out[0][0] += (1.0-costheta)*dr->x*dr->x;
  out[0][1] += (1.0-costheta)*dr->x*dr->y;
  out[0][2] += (1.0-costheta)*dr->x*dr->z;
  out[1][0] += (1.0-costheta)*dr->y*dr->x;
  out[1][1] += (1.0-costheta)*dr->y*dr->y;
  out[1][2] += (1.0-costheta)*dr->y*dr->z;
  out[2][0] += (1.0-costheta)*dr->z*dr->x;
  out[2][1] += (1.0-costheta)*dr->z*dr->y;
  out[2][2] += (1.0-costheta)*dr->z*dr->z;
  
  out[0][1] += -ds*dr->z;
  out[0][2] += ds*dr->y;
  out[1][0] += ds*dr->z;
  out[1][2] += -ds*dr->x;
  out[2][0] += -ds*dr->y;
  out[2][1] += ds*dr->x;
  
  
  


}
void PointTool::setMatrixVector(Point3d *out,double A[3][3],Point3d *dv){
  out->x = A[0][0]*dv->x+A[0][1]*dv->y+A[0][2]*dv->z;
  out->y = A[1][0]*dv->x+A[1][1]*dv->y+A[1][2]*dv->z;
  out->z = A[2][0]*dv->x+A[2][1]*dv->y+A[2][2]*dv->z;
  
  

}
double PointTool::LineDistance(Point3d *cen,Point3d *in1,Point3d *in2){
  
  double dist=0.0;
  makeVector(ddv1,cen,in1);
  makeVector(ddv2,in1,in2);dist = Point3dSize(ddv2);
  double st =-((InnerProduct(ddv1,ddv2))/dist);
  if(st<=0.0)return dist;
  if(st>=1.0)return Distance(cen,in2);
  CrossVector(ddv3,ddv2,ddv1);
  return (Point3dSize(ddv3)/dist);
}
  
void PointTool::setFaceNormal(Point3d *out,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc){
  makeVector(ddv1,v1,v2);
  makeVector(ddv2,v2,v3);
  CrossVector(out,ddv1,ddv2);
  Normalize3D(out);
}
void PointTool::setFaceNormal(Point3d *out,Point3d *v1,Point3d *v2,Point3d *v3){
  makeVector(ddv1,v1,v2);
  makeVector(ddv2,v2,v3);
  CrossVector(out,ddv1,ddv2);
  Normalize3D(out);
}
void PointTool::setFaceNormalArea(Point3d *out,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc){
  makeVector(bc[0],v1,v2);
  makeVector(bc[1],v2,v3);
  CrossVector(out,bc[0],bc[1]);
}

void PointTool::setInter(Point3d *out,double alpha,double beta,Point3d *dv1,Point3d *dv2){
  out->x = (alpha*dv2->x+beta*dv1->x)/(alpha+beta);
  out->y = (alpha*dv2->y+beta*dv1->y)/(alpha+beta);
  out->z = (alpha*dv2->z+beta*dv1->z)/(alpha+beta);
}
void PointTool::setMeanValue(double w[4],Point3d *ev,Point3d *v1,Point3d *v2,Point3d *v3,Point3d *v4,Point3d **bc){
  double leng[4],angle[4],alpha[4],wsum,ww1[4],ww2[4];
  wsum = 0.0;

  makeVector(bc[0],ev,v1);
  leng[0] = Point3dSize(bc[0]);if(leng[0]==0.0)leng[0]=1.0;
  makeVector(bc[1],ev,v2);
  leng[1] = Point3dSize(bc[1]);if(leng[1]==0.0)leng[1]=1.0;
  angle[0] = InnerProduct(bc[0],bc[1])/(leng[0]*leng[1]);
  
  makeVector(bc[2],ev,v3);
  leng[2] = Point3dSize(bc[2]);if(leng[2]==0.0)leng[2]=1.0;
  angle[1] = InnerProduct(bc[1],bc[2])/(leng[1]*leng[2]);
  
  
  makeVector(bc[3],ev,v4);
  leng[3] = Point3dSize(bc[3]);if(leng[3]==0.0)leng[3]=1.0;
  angle[2] = InnerProduct(bc[2],bc[3])/(leng[2]*leng[3]);
  angle[3] = InnerProduct(bc[3],bc[0])/(leng[3]*leng[0]);
  
  
  
  
  int i;
  int checkedge=1;
  int checkID=0;
  for(i=0;i<4;i++){
    if(1.0+angle[i]==0.0){
      alpha[i] = 0.0;checkedge=0;checkID=i;
      break;
    }else{
      alpha[i] = sqrt(((1.0-angle[i])/(1.0+angle[i])));
      
    }
  }
  if(checkedge==0){
    if(checkID==0){
      w[2]=0.0;w[3]=0.0;
      w[0] = leng[1]/(leng[0]+leng[1]);
      w[1] = leng[0]/(leng[0]+leng[1]);
      
      
    }else if(checkID==1){
      w[0]=0.0;w[3]=0.0;
      w[1] = leng[2]/(leng[1]+leng[2]);
      w[2] = leng[1]/(leng[1]+leng[2]);
      

    }else if(checkID==2){
      w[0]=0.0;w[1]=0.0;
      w[2] = leng[3]/(leng[2]+leng[3]);
      w[3] = leng[2]/(leng[2]+leng[3]);
      
    }else{
      w[1]=0.0;w[2]=0.0;
      w[3] = leng[0]/(leng[0]+leng[3]);
      w[0] = leng[3]/(leng[0]+leng[3]);
      
    }
    
  }else{
  ww1[0] =alpha[0]/leng[0];
  ww2[0] =alpha[0]/leng[1];
  
  wsum += ww1[0] + ww2[0];
  ww1[1] =alpha[1]/leng[1];
  ww2[1] =alpha[1]/leng[2];
  
  wsum += ww1[1] + ww2[1];
  
  ww1[2] =alpha[2]/leng[2];
  ww2[2] =alpha[2]/leng[3];
  
  wsum += ww1[2] + ww2[2];
  ww1[3] =alpha[3]/leng[3];
  ww2[3] =alpha[3]/leng[0];
  
  wsum += ww1[3] + ww2[3];
  if(wsum==0.0)wsum=1.0;
  w[0] = (ww1[0] + ww2[3])/wsum;
  w[1] = (ww1[1] + ww2[0])/wsum;
  w[2] = (ww1[2] + ww2[1])/wsum;
  w[3] = (ww1[3] + ww2[2])/wsum;
  }
}


void PointTool::setTIIntersect(Point3d *out,Point3d *fn,Point3d *fp,Point3d *lt,Point3d *lp){
  double d = -InnerProduct(fn,fp);
  double div = InnerProduct(fn,lt);if(div==0.0)div=1.0;
  double t = -((d+InnerProduct(fn,lp))/div);
  out->x = lp->x + t*lt->x;
  out->y = lp->y + t*lt->y;
  out->z = lp->z + t*lt->z;
  

}
void PointTool::Normalize3D(Point3d* inout){
  double dsize = this->Point3dSize(inout);
  if(dsize==0.0)dsize=1.0;
  inout->x /= dsize;
  inout->y /= dsize;
  inout->z /= dsize;
  

}
double PointTool::Point2DSize(double d1x,double d1y,double d2x,double d2y){
  double ddx = d1x-d2x;
  double ddy = d1y-d2y;
  return sqrt((ddx*ddx+ddy*ddy));

}
double PointTool::Point2DSizeSq(double d1x,double d1y,double d2x,double d2y){
  double ddx = d1x-d2x;
  double ddy = d1y-d2y;
  return (ddx*ddx+ddy*ddy);

}
double PointTool::HessGrad(Point3d* Suu,Point3d *Su,Point3d* Suv,Point3d *Sv){

  return (InnerProduct(Suu,Su)+InnerProduct(Suv,Sv));
}
void PointTool::ScalarVector(Point3d *out,double dv,Point3d *in){
  out->x = dv*in->x;
  out->y = dv*in->y;
  out->z = dv*in->z;
  

}
void PointTool::setParametricDs(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double t1,double t2,double t3,double A){
  double dt1 = t2-t3;
  double dt2 = t3-t1;
  double dt3 = t1-t2;
  out->x = ((q1->x)*dt1 + (q2->x)*dt2 + (q3->x)*dt3)/(2.0*A);
  out->y = ((q1->y)*dt1 + (q2->y)*dt2 + (q3->y)*dt3)/(2.0*A);
  out->z = ((q1->z)*dt1 + (q2->z)*dt2 + (q3->z)*dt3)/(2.0*A);
  
}
void PointTool::setParametricDt(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double s1,double s2,double s3,double A){
  double ds1 = s3-s2;
  double ds2 = s1-s3;
  double ds3 = s2-s1;
  out->x = ((q1->x)*ds1 + (q2->x)*ds2 + (q3->x)*ds3)/(2.0*A);
  out->y = ((q1->y)*ds1 + (q2->y)*ds2 + (q3->y)*ds3)/(2.0*A);
  out->z = ((q1->z)*ds1 + (q2->z)*ds2 + (q3->z)*ds3)/(2.0*A);

}
void PointTool::setParametricDss(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double t1,double t2,double t3){
  double dt1 = t2-t3;
  double dt2 = t3-t1;
  double dt3 = t1-t2;
  out->x = ((q1->x)*dt1 + (q2->x)*dt2 + (q3->x)*dt3);
  out->y = ((q1->y)*dt1 + (q2->y)*dt2 + (q3->y)*dt3);
  out->z = ((q1->z)*dt1 + (q2->z)*dt2 + (q3->z)*dt3);
  
}
void PointTool::setParametricDtt(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double s1,double s2,double s3){
  double ds1 = s3-s2;
  double ds2 = s1-s3;
  double ds3 = s2-s1;
  out->x = ((q1->x)*ds1 + (q2->x)*ds2 + (q3->x)*ds3);
  out->y = ((q1->y)*ds1 + (q2->y)*ds2 + (q3->y)*ds3);
  out->z = ((q1->z)*ds1 + (q2->z)*ds2 + (q3->z)*ds3);

}

double PointTool::getParametricA(double t1,double t2,double t3,double s1,double s2,double s3){
  double da = ((s2-s1)*(t3-t1)-(s3-s1)*(t2-t1))/2.0;
  if(da==0.0)return 1.0;
  return da;
}  
void PointTool::setBCVec(Point3d *out,Point3d *ev,Point3d *v1,Point3d *v2,Point3d *v3){
  double div = (-((v1->z)*(v2->y)*(v3->x))+((v1->y)*(v2->z)*(v3->x))+((v1->z)*(v2->x)*(v3->y))-((v1->x)*(v2->z)*(v3->y))-((v1->y)*(v2->x)*(v3->z))+((v1->x)*(v2->y)*(v3->z)));if(div==0.0)div=1.0;
  out->x = -(((v2->z)*(v3->y)*(ev->x))-((v2->y)*(v3->z)*(ev->x))-((v2->z)*(v3->x)*(ev->y))+((v2->x)*(v3->z)*(ev->y))+((v2->y)*(v3->x)*(ev->z))-((v2->x)*(v3->y)*(ev->z)))/div;
  out->y = -(((v1->z)*(v3->y)*(ev->x))-((v1->y)*(v3->z)*(ev->x))-((v1->z)*(v3->x)*(ev->y))+((v1->x)*(v3->z)*(ev->y))+((v1->y)*(v3->x)*(ev->z))-((v1->x)*(v3->y)*(ev->z)))/div;
  out->z = -(((v1->z)*(v2->y)*(ev->x))-((v1->y)*(v2->z)*(ev->x))-((v1->z)*(v2->x)*(ev->y))+((v1->x)*(v2->z)*(ev->y))+((v1->y)*(v2->x)*(ev->z))-((v1->x)*(v2->y)*(ev->z)))/div;
}
double PointTool::Distance(Point3d *in1,Point3d *in2){
  return sqrt(((((in1->x)-(in2->x))*((in1->x)-(in2->x)))
	       +(((in1->y)-(in2->y))*((in1->y)-(in2->y)))
	       +(((in1->z)-(in2->z))*((in1->z)-(in2->z)))));
		      
  }
void PointTool::makeVector(Point3d *out,Point3d *in1,Point3d *in2){
    out->x = ((in2->x) - (in1->x));
    out->y = ((in2->y) - (in1->y));
    out->z = ((in2->z) - (in1->z));
}
void PointTool::CrossVector(Point3d *out,Point3d *in1,Point3d *in2){
  out->x = ((in1->y)*(in2->z) - (in2->y)*(in1->z));
  out->y = ((in1->z)*(in2->x) - (in2->z)*(in1->x));
  out->z = ((in1->x)*(in2->y) - (in2->x)*(in1->y));
  }
double PointTool::InnerProduct(Point3d *in1,Point3d *in2){
    return ((in1->x)*(in2->x) + (in1->y)*(in2->y) + (in1->z)*(in2->z));
  }
void PointTool::setCenter(Point3d *out,Point3d *in1,Point3d *in2,Point3d *in3){
  out->x = (in1->x+in2->x+in3->x)/3.0;
  out->y = (in1->y+in2->y+in3->y)/3.0;
  out->z = (in1->z+in2->z+in3->z)/3.0;
  
}
double PointTool::Point3dSize(Point3d *in){
    return sqrt(((in->x)*(in->x)+(in->y)*(in->y)+(in->z)*(in->z)));
  }
void PointTool::Projection(Point3d *out,Point3d *in,Point3d *normal){
  double dv = InnerProduct(in,normal);
  out->x = in->x - dv*normal->x;
  out->y = in->y - dv*normal->y;
  out->z = in->z - dv*normal->z;
} 
void PointTool::setVitrualTangents(Point3d *normal,Point3d *outt1,Point3d *outt2){
  
  if(fabs(normal->x)<fabs(normal->y)){
    outt1->x = 0.0;
    outt1->y = -normal->z;
    outt1->z = normal->y;
    
  }else{
    outt1->x = normal->z;
    outt1->y = 0.0;
    outt1->z = -normal->x;
  }
  this->Normalize3D(outt1);
  this->CrossVector(outt2,normal,outt1);
  this->Normalize3D(outt2);
  
  /*
  outt1->x =  normal->y;
  outt1->y = -(normal->x);
  outt1->z = 0.0;
  this->Normalize3D(outt1);
  this->CrossVector(outt2,normal,outt1);
  this->Normalize3D(outt2);
  */
}

int PointTool::getBraycentricCQuad(Point3d *evaluation,Point3d *baryparam,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc){
  double dsize=0.0;
  double InterRadi=0.0;
  
      
  this->makeVector(bc[0],v1,v2);
  this->makeVector(bc[1],v2,v3);
  this->makeVector(bc[2],v3,v1);
  this->CrossVector(bc[3],bc[0],bc[1]);
  this->Normalize3D(bc[3]);
  
  
  makeVector(bc[4],v1,evaluation);
  this->CrossVector(bc[5],bc[4],bc[0]);
  if(this->InnerProduct(bc[5],bc[3])>0.0)return 1;
  makeVector(bc[4],v2,evaluation);
  this->CrossVector(bc[5],bc[4],bc[1]);
  if(this->InnerProduct(bc[5],bc[3])>0.0)return 1;
  makeVector(bc[4],v3,evaluation);
  this->CrossVector(bc[5],bc[4],bc[2]);
  if(this->InnerProduct(bc[5],bc[3])>0.0)return 1;
  

  
  return 0;

}
/*
int PointTool::getBraycentricCQuad(Point3d *evaluation,Point3d *baryparam,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc){
  double dsize=0.0;
  double InterRadi=0.0;
  double darea=0.0;
      
  this->makeVector(bc[0],v1,v2);
  this->makeVector(bc[1],v2,v3);
  this->makeVector(bc[2],v3,v1);
  this->CrossVector(bc[3],bc[0],bc[1]);
  dsize = Point3dSize(bc[3]);
  if(dsize==0.0)dsize=1.0;
  darea = dsize;
    

  bc[8]->x = (bc[3]->x/dsize);
  bc[8]->y = (bc[3]->y/dsize);
  bc[8]->z = (bc[3]->z/dsize);
  makeVector(bc[7],v1,evaluation);
  dsize = this->InnerProduct(bc[7],bc[8]);
  bc[9]->x = (evaluation->x - dsize*bc[8]->x);
  bc[9]->y = (evaluation->y - dsize*bc[8]->y);
  bc[9]->z = (evaluation->z - dsize*bc[8]->z);
  this->makeVector(bc[4],v3,bc[9]);
  this->CrossVector(bc[3],bc[1],bc[4]);
  dsize = this->Point3dSize(bc[3]);
  if(this->InnerProduct(bc[8],bc[3])<0.0)dsize *=-1.0;
  baryparam->x = (dsize/darea);
  this->makeVector(bc[5],v1,bc[9]);
  this->CrossVector(bc[3],bc[2],bc[5]);
  dsize = this->Point3dSize(bc[3]);
  if(this->InnerProduct(bc[8],bc[3])<0.0)dsize *=-1.0;
  baryparam->y = (dsize/darea);
  this->makeVector(bc[6],v2,bc[9]);
  this->CrossVector(bc[3],bc[0],bc[6]);
  dsize = this->Point3dSize(bc[3]);
  if(this->InnerProduct(bc[8],bc[3])<0.0)dsize *=-1.0;
  baryparam->z = (dsize/darea);
  
  if((baryparam->x >= 0.0) && (baryparam->y >=0.0) && (baryparam->z >=0.0)){
    return 1;
      }else{
      return 0;
    }
  return 0;

  }*/


int PointTool::getBraycentricC(Point3d *evaluation,Point3d *baryparam,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc){
  double tA=0.0;
  double A1=0.0;
  double A2=0.0;
  double A3=0.0;
    
  this->makeVector(bc[0],evaluation,v1);
  this->makeVector(bc[1],evaluation,v2);
  this->makeVector(bc[2],evaluation,v3);
  
  this->CrossVector(bc[3],bc[0],bc[1]);
  A3 = this->Point3dSize(bc[3]);
  this->CrossVector(bc[3],bc[0],bc[2]);
  A2 = this->Point3dSize(bc[3]);
  this->CrossVector(bc[3],bc[1],bc[2]);
  A1 = this->Point3dSize(bc[3]);
  tA = A1+A2+A3;if(tA==0.0)tA=1.0;
  baryparam->x = A1/tA;
  baryparam->y = A2/tA;
  baryparam->z = A3/tA;
    
  if((baryparam->x >= 0.0) && (baryparam->y >=0.0) && (baryparam->z >=0.0)){
    return 1;
  
    }else{
      return 0;
    }
  return 0;

}

void PointTool::setFloydC(int sx,int sy,int i,int j,double **tmp,double error){
  if(j+1<sx){
    if(i<sy-1){
      tmp[i][j+1] += ((7.0/16.0)*error);
    }else{
      tmp[i][j+1] += error;
    }
  }else{
    if(i<sy-1)
      tmp[i+1][0] += ((7.0/16.0)*error);
  }
  if(i<sy-1){
    tmp[i+1][j] += ((5.0/16.0)*error);
    if(j+1<sx){
      tmp[i+1][j+1] += ((1.0/16.0)*error);
    }else{
      
      if(i<sy-2){
	tmp[i+2][0] += ((1.0/16.0)*error);
      }else{
	if(j+1<sx)
	  tmp[i][j+1] += ((1.0/16.0)*error);
      }
    }
    if(j>=1){
      tmp[i+1][j-1] += ((3.0/16.0)*error);
    }else{
      tmp[i][sx-1] += ((3.0/16.0)*error);
    }
  }
}
/*
void PointTool::set3Dpoint(Point3d* out,Point3d *in,double *pU,double *pV,Point3d **point,int **Face,IDrendering *myrender, Point3d **bc, Point3d **dbc,int imgN){
  int i = ((int)((-in->y+((double)(imgN))-1.0+0.5)));
  int j = ((int)((in->x-0.5)));
    
  if(i<0)i=0;
  if(i>=imgN)i = imgN-1;
  if(j<0)j=0;
  if(j>=imgN)j = imgN-1;
  
  
  bc[1]->x = pU[Face[myrender->imgF[i][j]][0]];
  bc[1]->y = pV[Face[myrender->imgF[i][j]][0]];
  bc[1]->z = 0.0;
  bc[2]->x = pU[Face[myrender->imgF[i][j]][1]];
  bc[2]->y = pV[Face[myrender->imgF[i][j]][1]];
  bc[2]->z = 0.0;
  bc[3]->x = pU[Face[myrender->imgF[i][j]][2]];
  bc[3]->y = pV[Face[myrender->imgF[i][j]][2]];
  bc[3]->z = 0.0;
  getBraycentricC(in,bc[4],bc[1],bc[2],bc[3],dbc);
  out->x = (bc[4]->x)*point[Face[myrender->imgF[i][j]][0]]->x
    + (bc[4]->y)*point[Face[myrender->imgF[i][j]][1]]->x
    + (bc[4]->z)*point[Face[myrender->imgF[i][j]][2]]->x;
  out->y = (bc[4]->x)*point[Face[myrender->imgF[i][j]][0]]->y
    + (bc[4]->y)*point[Face[myrender->imgF[i][j]][1]]->y
    + (bc[4]->z)*point[Face[myrender->imgF[i][j]][2]]->y;
  out->z = (bc[4]->x)*point[Face[myrender->imgF[i][j]][0]]->z
    + (bc[4]->y)*point[Face[myrender->imgF[i][j]][1]]->z
    + (bc[4]->z)*point[Face[myrender->imgF[i][j]][2]]->z;
}
*/ 
double PointTool::getArea(Point3d* dv1,Point3d* dv2,Point3d* dv3){
  makeVector(ddv1,dv2,dv1);
  makeVector(ddv2,dv3,dv1);
  CrossVector(ddv3,ddv1,ddv2);
  return (Point3dSize(ddv3)/2.0);

}
