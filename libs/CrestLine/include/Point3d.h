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
class Point3d{
 public:
  double x,y,z;
  Point3d(){x=0.0;y=0.0;z=0.0;}
  Point3d(double dx,double dy,double dz){x=dx;y=dy;z=dz;}
  virtual ~Point3d(){x = 0.0;y=0.0;z=0.0;}
 

private:
  Point3d(const Point3d& rhs);
  const Point3d &operator=(const Point3d& rhs);
};
