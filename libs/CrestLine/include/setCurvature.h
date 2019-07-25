#ifndef SET_CURVATURE
#define SET_CURVATURE

void getCurvature(int*& vids,unsigned vidNum,int**& vadj,int*& vadjNum,double*& vertices, unsigned vNum, unsigned*& faces, unsigned fNum,unsigned neighboreSize,double*& mg1,double*& mg2,double*& dr1,double*& dr2);

void crestLineGen(int**& vadj,int*& vadjNum,double*& vertices, unsigned vNum, unsigned*& faces, unsigned fNum, unsigned neighboreSize, unsigned ridgeTension,
	unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval);

#endif