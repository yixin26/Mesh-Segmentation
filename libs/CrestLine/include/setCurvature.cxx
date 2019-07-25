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
#include "setCurvature.h"
#include<stdio.h>
#include<math.h>
#include"Point3d.h"
#include"PointTool.h"
#include"IDList.h"
#include"PolarList.h"
#include"IDSet.h"
#include"Polyhedron.h"



void crestLineGen(int**& vadj,int*& vadjNum,double*& vertices, unsigned vNum, unsigned*& faces, unsigned fNum, unsigned neighboreSize, unsigned ridgeTension,
	unsigned &pNum,unsigned &cNum,unsigned &eNum, double* &fPoints,unsigned* &fPointType, double* &fFeature, unsigned* &fEdges, unsigned &interval)
{
	Polyhedron *mymesh = new Polyhedron();

	mymesh->readmesh(vadj,vadjNum,vertices, vNum, faces, fNum, neighboreSize, ridgeTension,
		pNum,cNum,eNum,fPoints,fPointType,fFeature,fEdges,interval);

	delete mymesh;
}

void getCurvature(int*& vids,unsigned vidNum,int**& vadj,int*& vadjNum,double*& vertices, unsigned vNum, unsigned*& faces, unsigned fNum, unsigned neighboreSize,double*& mg1,double*& mg2,double*& dr1,double*& dr2)
{
	Polyhedron *mymesh = new Polyhedron();

	mymesh->getCurvature(vids,vidNum,vadj,vadjNum,vertices, vNum, faces, fNum, neighboreSize,mg1,mg2,dr1,dr2);

	delete mymesh;
}

