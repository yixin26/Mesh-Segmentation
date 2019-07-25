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
class SvdSolve{
  
 public:
  SvdSolve(){
  }
  virtual ~SvdSolve(){
  }
  void svdcmp(double **a, int m, int n, double w[], double **v);
  void svbksb(double **u,double *w,double **v,int m,int n,double b[],double x[]);
  
  double pythag(double a, double b);
 private:
  SvdSolve(const SvdSolve& rhs);
  const SvdSolve &operator=(const SvdSolve& rhs);
  
};

