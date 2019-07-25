/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#ifndef TENSOR_H
#define TENSOR_H

#define _USE_MATH_DEFINES 
#include <math.h>

#include <vector>
using namespace std;

#include "myMesh.h"

namespace GeoProperty{

Vec3 basicNormalTransport(const Vec3 &e1, const Vec3 &e2, const Vec3 &n);

class Tensor{
public:

	Tensor(){ mag1 = mag2 =0; mat[0][0]=mat[0][1]=mat[1][0]=mat[1][1]=0; }
	~Tensor(){}

	static void makeTensor(const Vec3& d1,const Vec3& d2,Tensor &t);

	static void averageTensor(Tensor &t1,Tensor &t2,Tensor &t);
	static void averageTensor(Tensor &t1,Tensor &t2,Tensor &t, double rat);
	static void averageTensor(std::vector<Tensor> &ts,Tensor &t, double rat);

	static void averageTensor(std::vector<Tensor> &ts,Tensor &t);
	static void computeTriangleCurvature(MyMesh*,std::vector<Tensor>&);
// 	static void computeTriangleTensor(MyMesh* orgMesh, const std::vector<Tensor>& vertexTensor, std::vector<Tensor>& faceTensor);
	
	static double compareTensors(Tensor& t1, Tensor& t2, std::vector<Vec2> vs);

public:
	Vec3 dir1,dir2;
	double mag1,mag2;
	double mat[2][2];
};

double VectorMul(Vec2& p11, Vec2& p12, Vec2& p21, Vec2& p22);
double IsPointInTriangle(Vec3 tp1, Vec3 tp2, Vec3 tp3, Vec3 testp, Vec3& projp = Vec3(0, 0, 0));
Tensor lineSearchHessian(Vec3& v, Tensor& tensor);

inline Vec3 basicNormalTransport(const Vec3 &e1, const Vec3 &e2, const Vec3 &n)
{
	Vec3 transNorm = n;
	Vec3  tangent = (e1 + e2);
	tangent.normalize();
	transNorm -= 2 * (transNorm.dot(tangent))*tangent;
	transNorm.normalize();
	return transNorm;
}
inline double VectorMul(Vec2& p11, Vec2& p12, Vec2& p21, Vec2& p22)
{
	return ((p12.x - p11.x)*(p22.y - p21.y) - (p22.x - p21.x)*(p12.y - p11.y));
}
inline double IsPointInTriangle(Vec3 tp1, Vec3 tp2, Vec3 tp3, Vec3 testp, Vec3& projp)
{
	Vec3 ax = tp2 - tp1; ax.normalize();
	Vec3 ay = tp3 - tp1;
	Vec3 nm = ax.cross(ay); nm.normalize();
	ay = nm.cross(ax); ay.normalize();

	testp -= tp1;
	//set org as orgin, and proj triangle vertices to plane;
	tp2 -= tp1;
	tp3 -= tp1;
	//tp1 -= tp1;
	//Vec2 pv1 = Vec2(tp1.dot(ax),tp1.dot(ay));
	Vec2 pv1(0, 0);
	Vec2 pv2 = Vec2(tp2.dot(ax), tp2.dot(ay));
	Vec2 pv3 = Vec2(tp3.dot(ax), tp3.dot(ay));

	Vec2 testpv = Vec2(testp.dot(ax), testp.dot(ay));

	double TMul1 = VectorMul(pv1, testpv, pv2, testpv); if (abs(TMul1) < 1e-10) TMul1 = 0;
	double TMul2 = VectorMul(pv2, testpv, pv3, testpv); if (abs(TMul2) < 1e-10) TMul2 = 0;
	double TMul3 = VectorMul(pv3, testpv, pv1, testpv); if (abs(TMul3) < 1e-10) TMul3 = 0;

	//cout<<"p2:"<<pv2.x<<" "<<pv2.y<<" p3:"<<pv3.x<<" "<<pv3.y<<" pt:"<<testpv.x<<" "<<testpv.y<<" tm:"<<TMul1<<" "<<TMul2<<" "<<TMul3<<endl;
	if (TMul1 * TMul2 < 0 || TMul1 * TMul3 < 0 || TMul2 * TMul3 < 0)//outside
		return -1;

	projp = tp1 + testpv.x * ax + testpv.y * ay;
	// compute the distance
	return abs(nm.dot(testp));
}
}


#endif