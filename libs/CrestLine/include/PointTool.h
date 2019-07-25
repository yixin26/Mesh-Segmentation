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

class PointTool{
 public:
  Point3d *ddv1;
  Point3d *ddv2;
  Point3d *ddv3;
  
  PointTool(){
    ddv1 = new Point3d(0.0,0.0,0.0);
    ddv2 = new Point3d(0.0,0.0,0.0);
    ddv3 = new Point3d(0.0,0.0,0.0);
  }
  virtual ~PointTool(){
    delete ddv1;
    delete ddv2;
    delete ddv3;
  }
 private:
  PointTool(const PointTool& rhs);
  const PointTool &operator=(const PointTool& rhs);
 public:
  void setInverse(double out[3][3],double in[3][3]);
  void setJTJ(double out[3][3],double in[3][3]);
  double getTraceJTJ(double in[3][3]);
  double DetCC(double **in);
  void setRodriguesM(double out[3][3], Point3d *dr, double costheta);
  void setNormalF(Point3d* out,double dU,double dV,double *Func);
  void setMatrixVector(Point3d *out,double A[3][3],Point3d *dv);
  void setFaceNormal(Point3d *out,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc);
  void setFaceNormal(Point3d *out,Point3d *v1,Point3d *v2,Point3d *v3);
  void setFaceNormalArea(Point3d *out,Point3d *v1,Point3d *v2,Point3d *v3,Point3d **bc);
  
  void setInter(Point3d *out,double alpha,double beta,Point3d *dv1,Point3d *dv2);
 void setMeanValue(double w[4],Point3d *ev,Point3d *v1,Point3d *v2,Point3d *v3,Point3d *v4,Point3d **bc);
 //void set3Dpoint(Point3d* out,Point3d *in,double *pU,double *pV,Point3d **point,int **Face,IDrendering *myrender, Point3d **bc, Point3d **dbc,int imgN);
  double getArea(Point3d* dv1,Point3d* dv2,Point3d* dv3);
  double HessGrad(Point3d* Suu,Point3d *Su,Point3d* Suv,Point3d *Sv);
  void Normalize3D(Point3d* inout);
  double Point2DSize(double d1x,double d1y,double d2x,double d2y);
  double Point2DSizeSq(double d1x,double d1y,double d2x,double d2y);
  double Distance(Point3d *in1,Point3d *in2);
  double LineDistance(Point3d *cen,Point3d *in1,Point3d *in2);
  
  void makeVector(Point3d *out,Point3d *in1,Point3d *in2);
  void CrossVector(Point3d *out,Point3d *in1,Point3d *in2);
  double InnerProduct(Point3d *in1,Point3d *in2);
  double Point3dSize(Point3d *in);
  void Projection(Point3d *out,Point3d *in,Point3d *normal);
  void setVitrualTangents(Point3d *normal,Point3d *outt1,Point3d *outt2); 
  int getBraycentricC(Point3d *evaluation,Point3d *baryparam,Point3d *v1,Point3d *v2,Point3d *v3, Point3d **bc);
 int getBraycentricCQuad(Point3d *evaluation,Point3d *baryparam,Point3d *v1,Point3d *v2,Point3d *v3, Point3d **bc);
  
  void setCenter(Point3d *out,Point3d *in1,Point3d *in2,Point3d *in3);
  void setBCVec(Point3d *out,Point3d *ev,Point3d *v1,Point3d *v2,Point3d *v3);
  void setTIIntersect(Point3d *out,Point3d *fn,Point3d *fp,Point3d *lt,Point3d *lp);
  void setParametricDs(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double t1,double t2,double t3,double A);
  void setParametricDt(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double s1,double s2,double s3,double A);

  void setParametricDss(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double t1,double t2,double t3);
  void setParametricDtt(Point3d *out,Point3d *q1,Point3d *q2,Point3d *q3,double s1,double s2,double s3);

  double getParametricA(double t1,double t2,double t3,double s1,double s2,double s3);
  void ScalarVector(Point3d *out,double dv,Point3d *in);
  void setFloydC(int sx,int sy,int i,int j,double **tmp,double error);
  
};
