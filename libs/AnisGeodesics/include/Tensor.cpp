#include "Tensor.h"
#include "AnisGeodesic.h"
#include <iostream>

#include <Eigen/Dense>
using namespace Eigen;

using namespace GeoProperty;


void Tensor::makeTensor(const Vec3& d1,const Vec3& d2,Tensor &t)
{
	Matrix2d hesMat = Eigen::Matrix2d::Constant(2, 2, 0.0);
	for (unsigned i = 0; i < 2; i++)
	{
		for (unsigned j = 0; j < 2; j++)
		{
			hesMat(i, j) = t.mat[i][j];
		}
	}

	SelfAdjointEigenSolver<Matrix2d> vd(hesMat);
	Vector2d egval = vd.eigenvalues();
	Matrix2d egvec = vd.eigenvectors();

	unsigned ind1 = 0, ind2 = 1;

	t.mag1 = abs(egval[ind1]); t.mag2 = abs(egval[ind2]);
	t.dir1 = egvec(ind1, 0) * d1 + egvec(ind1, 1) * d2; t.dir1.normalize();
	t.dir2 = egvec(ind2, 0) * d1 + egvec(ind2, 1) * d2; t.dir2.normalize();
}

void Tensor::averageTensor(Tensor &t1,Tensor &t2,Tensor &t)
{
	double rat1= 0.5;
	double rat2= 0.5;

	t1.mat[0][0]=t1.mag1; t1.mat[1][1]=t1.mag2;
	t1.mat[0][1]=t1.mat[1][0]=0;

	Vec2 v1_,v2_;
	v1_.x=t2.dir1.dot(t1.dir1); v1_.y=t2.dir1.dot(t1.dir2);
	v2_.x=t2.dir2.dot(t1.dir1); v2_.y=t2.dir2.dot(t1.dir2);

	t2.mat[0][0]= t2.mag1*v1_.x*v1_.x+t2.mag2*v2_.x*v2_.x;
	t2.mat[0][1]= t2.mat[1][0]= t2.mag1*v1_.x*v1_.y+t2.mag2*v2_.x*v2_.y;
	t2.mat[1][1]= t2.mag1*v1_.y*v1_.y+t2.mag2*v2_.y*v2_.y;

	Matrix2d hesMat = Eigen::Matrix2d::Constant(2, 2, 0.0);
	for (unsigned i = 0; i < 2; i++)
	{
		for (unsigned j = 0; j < 2; j++)
		{
			hesMat(i, j) = rat1*t1.mat[i][j] + rat2*t2.mat[i][j];
		}
	}

	SelfAdjointEigenSolver<Matrix2d> vd(hesMat);
	Vector2d egval = vd.eigenvalues();
	Matrix2d egvec = vd.eigenvectors();

	unsigned ind1 = 0, ind2 = 1;
	if (abs(egval[0]) > abs(egval[1])) //////////////////////////////
		swap(ind1, ind2);

	double nl = (abs(t1.mag1) + abs(t1.mag2) + abs(t2.mag1) + abs(t2.mag2));
	if (nl > 1e-10)
		nl /= (2.0*(abs(egval[ind1]) + abs(egval[ind2])));
	else
		nl = 1;

	t.mag1 = egval[ind1] * nl; t.mag2 = egval[ind2] * nl; //how to normalized it is very important;
	t.dir1 = egvec(ind1, 0) * t1.dir1 + egvec(ind1, 1) * t1.dir2; t.dir1.normalize();
	t.dir2 = egvec(ind2, 0) * t1.dir1 + egvec(ind2, 1) * t1.dir2; t.dir2.normalize();
}

void Tensor::averageTensor(std::vector<Tensor> &ts,Tensor &t, double rat)
{
	std::vector<double> rats(ts.size(),rat/double(ts.size()-1)); rats[0]=1-rat;
	Tensor& t1=ts.front();
	t1.mat[0][0]=t1.mag1; t1.mat[1][1]=t1.mag2;
	t1.mat[0][1]=t1.mat[1][0]=0;
	if(ts.size()==1) {t=t1; return;}

	Vec2 v1_,v2_;
	for(unsigned i=1;i<ts.size();i++)
	{
		Tensor& t2=ts[i];
		v1_.x=t2.dir1.dot(t1.dir1); v1_.y=t2.dir1.dot(t1.dir2);
		v2_.x=t2.dir2.dot(t1.dir1); v2_.y=t2.dir2.dot(t1.dir2);
		t2.mat[0][0]= t2.mag1*v1_.x*v1_.x+t2.mag2*v2_.x*v2_.x;
		t2.mat[0][1]= t2.mat[1][0]= t2.mag1*v1_.x*v1_.y+t2.mag2*v2_.x*v2_.y;
		t2.mat[1][1]= t2.mag1*v1_.y*v1_.y+t2.mag2*v2_.y*v2_.y;
	}

	Matrix2d hesMat = Eigen::Matrix2d::Constant(2, 2, 0.0);
	for(unsigned i=0;i<2;i++)
	{
		for (unsigned j = 0; j < 2; j++)
		{
			for (unsigned k = 0; k < ts.size(); k++)
			{
				hesMat(i, j) += rats[k] * ts[k].mat[i][j];
			}
		}
	}

	SelfAdjointEigenSolver<Matrix2d> vd(hesMat);
	Vector2d egval = vd.eigenvalues();
	Matrix2d egvec = vd.eigenvectors();

	unsigned ind1 = 0, ind2 = 1;
	if (abs(egval[0]) < abs(egval[1]))
		swap(ind1, ind2);

	double nl = 0;
	for (unsigned i = 0; i<ts.size(); i++)
		nl += rats[i] * abs(ts[i].mag1) + abs(ts[i].mag2);
	if (nl>1e-10)
		nl /= ts.size()*(abs(egval[ind1]) + abs(egval[ind2]));
	else
		nl = 1;

	t.mag1 = egval[ind1] * nl; t.mag2 = egval[ind2] * nl; //how to normalized it is very important;
	t.dir1 = egvec(ind1, 0) * t1.dir1 + egvec(ind1, 1) * t1.dir2; t.dir1.normalize();
	t.dir2 = egvec(ind2, 0) * t1.dir1 + egvec(ind2, 1) * t1.dir2; t.dir2.normalize();
}

void Tensor::averageTensor(Tensor &t1,Tensor &t2,Tensor &t, double rat)
{
	double rat1= 1-rat;

	t1.mat[0][0]=t1.mag1; t1.mat[1][1]=t1.mag2;
	t1.mat[0][1]=t1.mat[1][0]=0;

	Vec2 v1_,v2_;
	v1_.x=t2.dir1.dot(t1.dir1); v1_.y=t2.dir1.dot(t1.dir2);
	v2_.x=t2.dir2.dot(t1.dir1); v2_.y=t2.dir2.dot(t1.dir2);

	t2.mat[0][0]= t2.mag1*v1_.x*v1_.x+t2.mag2*v2_.x*v2_.x;
	t2.mat[0][1]= t2.mat[1][0]= t2.mag1*v1_.x*v1_.y+t2.mag2*v2_.x*v2_.y;
	t2.mat[1][1]= t2.mag1*v1_.y*v1_.y+t2.mag2*v2_.y*v2_.y;

	Matrix2d hesMat = Eigen::Matrix2d::Constant(2, 2, 0.0);
	for(unsigned i=0;i<2;i++)
	{
		for (unsigned j = 0; j < 2; j++)
		{
			hesMat(i, j) = rat1*t1.mat[i][j] + rat*t2.mat[i][j];
		}
	}

	SelfAdjointEigenSolver<Matrix2d> vd(hesMat);
	Vector2d egval = vd.eigenvalues();
	Matrix2d egvec = vd.eigenvectors();

	unsigned ind1 = 0, ind2 = 1;
	if ((abs(egval[0]) < abs(egval[1])) && (t1.mag1 >= t1.mag2) || (abs(egval[0]) > abs(egval[1])) && (t1.mag1<t1.mag2))
		swap(ind1, ind2);

	double nl = (abs(t1.mag1) + abs(t1.mag2) + abs(t2.mag1) + abs(t2.mag2));
	if (nl>1e-10)
		nl /= (2.0*(abs(egval[ind1]) + abs(egval[ind2])));
	else
		nl = 1;

	t.mag1 = egval[ind1] * nl; t.mag2 = egval[ind2] * nl; //how to normalized it is very important;
	t.dir1 = egvec(ind1, 0) * t1.dir1 + egvec(ind1, 1) * t1.dir2; t.dir1.normalize();
	t.dir2 = egvec(ind2, 0) * t1.dir1 + egvec(ind2, 1) * t1.dir2; t.dir2.normalize();
}

void Tensor::averageTensor(std::vector<Tensor> &ts,Tensor &t)
{
	if (ts.empty()) return;

	Tensor& t1=ts.front();
	t1.mat[0][0]=t1.mag1; t1.mat[1][1]=t1.mag2;
	t1.mat[0][1]=t1.mat[1][0]=0;

	if(ts.size()==1) 
	{
		t=t1; return;
	}
	
	for(unsigned i=1;i<ts.size();i++)
	{
		Tensor& t2=ts[i];

		Vec2 v1(t2.dir1.dot(t1.dir1), t2.dir1.dot(t1.dir2));
		Vec2 v2(t2.dir2.dot(t1.dir1), t2.dir2.dot(t1.dir2));

		t2.mat[0][0] = t2.mag1*v1.x*v1.x + t2.mag2*v2.x*v2.x;
		t2.mat[0][1] = t2.mat[1][0] = t2.mag1*v1.x*v1.y + t2.mag2*v2.x*v2.y;
		t2.mat[1][1] = t2.mag1*v1.y*v1.y + t2.mag2*v2.y*v2.y;
	}

	Matrix2d hesMat = Eigen::Matrix2d::Constant(2, 2, 0.0);
	for (unsigned k = 0; k < ts.size(); k++)
	{
		hesMat(0, 0) += ts[k].mat[0][0];
		hesMat(0, 1) += ts[k].mat[0][1];
		hesMat(1, 0) += ts[k].mat[1][0];
		hesMat(1, 1) += ts[k].mat[1][1];
	}
	SelfAdjointEigenSolver<Matrix2d> vd(hesMat);
	Vector2d egval = vd.eigenvalues();
	Matrix2d egvec = vd.eigenvectors();

	unsigned ind1 = 0, ind2 = 1;
	if (abs(egval[0]) < abs(egval[1]))
		swap(ind1, ind2);

	double nl = 0;
	for (unsigned i = 0; i<ts.size(); i++)
		nl += abs(ts[i].mag1) + abs(ts[i].mag2);
	if (nl>1e-10)
		nl /= ts.size()*(abs(egval[ind1]) + abs(egval[ind2]));
	else
		nl = 1;

	t.mag1 = egval[ind1] * nl; t.mag2 = egval[ind2] * nl; //how to normalized it is very important;
	t.dir1 = egvec(ind1, 0) * t1.dir1 + egvec(ind1, 1) * t1.dir2; t.dir1.normalize();
	t.dir2 = egvec(ind2, 0) * t1.dir1 + egvec(ind2, 1) * t1.dir2; t.dir2.normalize();
}

void Tensor::computeTriangleCurvature(MyMesh* orgMesh,std::vector<Tensor>& faceTensor)
{
	/*
	int fNum = orgMesh->getFaces().size();
	faceTensor.clear(); faceTensor.resize(fNum);
	const auto& fts = orgMesh->getFIter();
#pragma omp parallel for
	for(int i=0;i<fNum;i++){
		auto f_it = fts[i];
		MyMesh::VertexIter v[3]={f_it->vertex_iter(0),f_it->vertex_iter(1),f_it->vertex_iter(2)};

		// local coord system for single triangle;
		Vec3 axix = v[1]->coordinate()- v[0]->coordinate(); axix.normalize();
		Vec3 axiy = v[2]->coordinate()- v[0]->coordinate();
		Vec3 norm = axix.cross(axiy); norm.normalize();
		axiy = norm.cross(axix); axiy.normalize();

		std::vector<Tensor> ts(3);
		for(unsigned i=0;i<3;i++){
			ts[i].dir1 = basicNormalTransport(v[i]->normal(),norm,v[i]->direction(0));
			ts[i].dir2 = basicNormalTransport(v[i]->normal(),norm,v[i]->direction(1));
			ts[i].mag1 = v[0]->magnitude(0);
			ts[i].mag2 = v[0]->magnitude(1);
		}

		Tensor t;
		averageTensor(ts,t);
		faceTensor[i] =t;
/*  
	//second method;  which seems wrong, but whose result is close to method 1;
		double anis[2];
		Jacobi jcb(2);
		for(unsigned i=0;i<2;i++){
			std::vector<std::vector<double> > mat(2,std::vector<double>(2,0));
			for(unsigned j=0;j<3;j++){
				Vec3 projVec = v_its[j]->direction(i); projVec.normalize();
				projVec -= projVec.dot(norm) * norm; projVec.normalize();

				double x = v_its[j]->direction(i).dot(u);	
				double y = v_its[j]->direction(i).dot(v);

				mat[0][0] += v_its[j]->magnitude(i)*x*x;
				mat[0][1] += v_its[j]->magnitude(i)*x*y;
				mat[1][1] += v_its[j]->magnitude(i)*y*y;
			}
			mat[1][0] = mat[0][1];

			jcb.setMatrix(mat);
			jcb.run();

			std::vector<double> egvalue = jcb.getEigenvalues();

			std::vector<std::vector<double> > egvector = jcb.getEigenvectors();
			unsigned ind=0;
			if(egvalue[0]<egvalue[1]) ind=1 ;
			Vec3 td = u*egvector[ind][0] + v*egvector[ind][1]; td.normalize();

			anis[i]= *std::max_element(egvalue.begin(),egvalue.end());
		}
		faceAnis.push_back(abs(anis[1]-anis[0]));
*/
//	}

	faceTensor.clear();
	Tensor t;
	for(auto f_it = orgMesh->getFaces().begin(); f_it!=orgMesh->getFaces().end(); f_it++)
	{
		MyMesh::VertexIter v[3]={f_it->vertex_iter(0),f_it->vertex_iter(1),f_it->vertex_iter(2)};

		// local coord system for single triangle;
		Vec3 axix = v[1]->coordinate()- v[0]->coordinate(); axix.normalize();
		Vec3 axiy = v[2]->coordinate()- v[0]->coordinate();
		Vec3 norm = axix.cross(axiy); norm.normalize();
		axiy = norm.cross(axix); axiy.normalize();

		std::vector<Tensor> ts(3);
		for(unsigned i=0;i<3;i++)
		{
			ts[i].dir1 = basicNormalTransport(v[i]->normal(),norm,v[i]->direction(0));
			ts[i].dir2 = basicNormalTransport(v[i]->normal(),norm,v[i]->direction(1));
			ts[i].mag1 = v[i]->magnitude(0);
			ts[i].mag2 = v[i]->magnitude(1);
		}

		averageTensor(ts,t);
		faceTensor.push_back(t);
	}
}

Tensor GeoProperty::lineSearchHessian(Vec3& v, Tensor& tensor)
{
	double x = v.dot(tensor.dir1);
	double y = v.dot(tensor.dir2);

	/*
	double denom = pow(x,2)*tensor.mag1+pow(y,2)*tensor.mag2;
	double denomSqrt = sqrt(denom);
	double upscale = 1e5;

	tensor.mat[0][0]=tensor.mag1;
	tensor.mat[1][1]=tensor.mag2;
	tensor.mat[0][1]=tensor.mat[1][0]=0;
	Tensor res;
	for (int i = 0; i < 2; i++){
	for (int j = 0; j < 2; j++){
	double t1 = 0.0 , t2 = 0.0;
	for (int k = 0; k < 2; k++){
	t1 += tensor.mat[i][k] * vecProj[k];
	t2 += tensor.mat[j][k] * vecProj[k];
	}
	res.mat[i][j] = (0.5 * tensor.mat[i][j] * upscale) / (denomSqrt * upscale) -
	(0.5 * t1 * upscale * t2) / (denom * upscale * denomSqrt);
	}
	}
	*/

	//make it as a metric tensor

	double a = tensor.mag1;
	double d = tensor.mag2;

	double dis = a*x*x + d*y*y;
	double dis23 = sqrt(dis) * dis;

	Tensor res;
	if (dis23>1e-8)
	{
		res.mat[0][0] = a / dis - a*a*x*x / dis23;
		res.mat[0][1] = -a*d*x*y / dis23;
		res.mat[1][0] = res.mat[0][1];
		res.mat[1][1] = d / dis - d*d*y*y / dis23;
		Tensor::makeTensor(tensor.dir1, tensor.dir2, res);
	}
	return res;
}

// double Tensor::compareTensors(Tensor& t1, Tensor& t2, std::vector<Vec2> vs)
// {
// 	//compute the difference of two tensors.... it is hard to define the difference..
// /*
// 	U1 = {{a1, b1}, {c1, d1}} ;
// 	U2 = {{a2, b2}, {c2, d2}};
// 	V = {Cos[#], Sin[#]} &;
// 	Sum[Power[Sqrt[V[x].U1.V[x]] - Sqrt[V[x].U2.V[x]], 2], {x, 0, Pi,
// 		Pi/360}]
// */
// 	t1.mat[0][0]=t1.mag1; t1.mat[1][1]=t1.mag2;
// 	t1.mat[0][1]=t1.mat[1][0]=0;
// 
// 	Vec2 v1_,v2_;
// 	v1_.x=t2.dir1.dot(t1.dir1); v1_.y=t2.dir1.dot(t1.dir2);
// 	v2_.x=t2.dir2.dot(t1.dir1); v2_.y=t2.dir2.dot(t1.dir2);
// 	t2.mat[0][0]= t2.mag1*v1_.x*v1_.x+t2.mag2*v2_.x*v2_.x;
// 	t2.mat[0][1]= t2.mat[1][0]= t2.mag1*v1_.x*v1_.y+t2.mag2*v2_.x*v2_.y;
// 	t2.mat[1][1]= t2.mag1*v1_.y*v1_.y+t2.mag2*v2_.y*v2_.y;
// 
// 	//Sum[Power[Sqrt[V[x].U1.V[x]] - Sqrt[V[x].U2.V[x]], 2]
// 	double ds=0;
// 	for(unsigned i=0;i<vs.size();i++)
// 	{
// 		ds+= pow(
// 			sqrt(abs((vs[i].x*t1.mat[0][0] + vs[i].y*t1.mat[1][0])*vs[i].x + (vs[i].x*t1.mat[0][1] + vs[i].y*t1.mat[1][1])*vs[i].y)) -
// 			sqrt(abs((vs[i].x*t2.mat[0][0] + vs[i].y*t2.mat[1][0])*vs[i].x + (vs[i].x*t2.mat[0][1] + vs[i].y*t2.mat[1][1])*vs[i].y))
// 			,2);
// 	}
// /*
// 	for(unsigned i=0;i<vs.size();i++){
// 		ds+= abs(
// 			sqrt(abs((vs[i].x*t1.mat[0][0] + vs[i].y*t1.mat[1][0])*vs[i].x + (vs[i].x*t1.mat[0][1] + vs[i].y*t1.mat[1][1])*vs[i].y)) -
// 			sqrt(abs((vs[i].x*t2.mat[0][0] + vs[i].y*t2.mat[1][0])*vs[i].x + (vs[i].x*t2.mat[0][1] + vs[i].y*t2.mat[1][1])*vs[i].y))
// 			);
// 	}
// */
// 	return ds;
// }
