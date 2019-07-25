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
class Eigens{
 public:
  Eigens(){
  }
  virtual ~Eigens(){
  }
  void tqli(double d[], double e[], int n, double **z);
  void tred2(double **a, int n, double d[], double e[]);
  void jacobi(double **a, int n, double d[], double **v, int *nrot);
  void eigsrt(double d[], double **v, int n);
 private:
  double pythag(double a, double b);
  Eigens(const Eigens& rhs);
  const Eigens &operator=(const Eigens& rhs);

};
