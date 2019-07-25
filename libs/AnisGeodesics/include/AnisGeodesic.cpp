/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#include "geodesic_mesh.h"
#include "geodesic_algorithm_exact.h"
#include <algorithm>
#include <queue>
#include "AnisGeodesic.h"
#include "Tensor.h"

bool AnisGeodesic::isDoGlobalScale=false;
double AnisGeodesic::globalScaleParameter=10.0;
double AnisGeodesic::maxAnis=1.0;
unsigned AnisGeodesic::subTimesMax=100000;
double AnisGeodesic::worstTriangleCost=0.001;
unsigned AnisGeodesic::modifyNeighboreSize=10;
double AnisGeodesic::modifyFieldRatio=0.5;

using namespace GeoProperty;

double AnisGeodesic::computeAngle(const Vec3 &org1,const Vec3 &org2)
{
	Vec3 ang_str = org1;
	Vec3 ang_end = org2;
	ang_str.normalize();
	ang_end.normalize();
	double res = ang_str.dot(ang_end);
	res = res < -1.0 ? -1.0 : res;
	res = res > 1.0 ? 1.0 : res;
	return acos(res);
}
Vec3 AnisGeodesic::projectVecToPlane(const Vec3 &v,const Vec3 &n)
{
	Vec3 norm= n; norm.normalize();
	Vec3 projVec = v; projVec.normalize();
	projVec -= projVec.dot(norm) * norm;
	projVec.normalize();
	return projVec;
}

void AnisGeodesic::creatMidPoint(MyMesh::Vertex& newVer, MyMesh::EdgeIter& e_it)
{
	MyMesh::VertexIter v[]={e_it->vertex_iter(0),e_it->vertex_iter(1)};
	newVer.coordinate()=(v[0]->coordinate()+v[1]->coordinate())/2.;
	newVer.normal()=(v[0]->normal()+v[1]->normal())/2.;
	newVer.normal().normalize();

	Tensor t1,t2;
	t1.dir1 = basicNormalTransport(v[0]->normal(),newVer.normal(),v[0]->direction(0));
	t1.dir2 = basicNormalTransport(v[0]->normal(),newVer.normal(),v[0]->direction(1));
	t1.mag1 = v[0]->magnitude(0);
	t1.mag2 = v[0]->magnitude(1);

	t2.dir1 = basicNormalTransport(v[1]->normal(),newVer.normal(),v[1]->direction(0));
	t2.dir2 = basicNormalTransport(v[1]->normal(),newVer.normal(),v[1]->direction(1));
	t2.mag1 = v[1]->magnitude(0);
	t2.mag2 = v[1]->magnitude(1);

	Tensor t; 
	t.averageTensor(t1,t2,t);
	newVer.direction(0) = t.dir1;
	newVer.direction(1) = t.dir2;
	newVer.magnitude(0) = t.mag1;
	newVer.magnitude(1) = t.mag2;
}
void AnisGeodesic::newPointDirection(MyMesh::VertexIter& v_it, std::vector<MyMesh::VertexIter>& sv_it, double rat)
{
	std::vector<Tensor> ts(sv_it.size()+1);
	Tensor& t1 = ts[0];
	t1.dir1 = v_it->direction(0);
	t1.dir2 = v_it->direction(1);
	t1.mag1 = v_it->magnitude(0);
	t1.mag2 = v_it->magnitude(1);
	for(unsigned i=0;i<ts.size()-1;i++)
	{
		Tensor& tempt=ts[i+1];

		tempt.dir1 = basicNormalTransport(sv_it[i]->normal(),v_it->normal(),sv_it[i]->direction(0));
		if(computeAngle(tempt.dir1,t1.dir1)>M_PI_2) 
			tempt.dir1 = -tempt.dir1;

		tempt.dir2 = basicNormalTransport(sv_it[i]->normal(),v_it->normal(),sv_it[i]->direction(1));
		if(computeAngle(tempt.dir2,t1.dir2)>M_PI_2) 
			tempt.dir2= -tempt.dir2;

		tempt.mag1 = sv_it[i]->magnitude(0);
		tempt.mag2 = sv_it[i]->magnitude(1);
	}

	Tensor t; 
	t.averageTensor(ts,t,rat);
	v_it->direction(0) = t.dir1;
	v_it->direction(1) = t.dir2;
	v_it->magnitude(0) = t.mag1;
	v_it->magnitude(1) = t.mag2;
}
void AnisGeodesic::computeEdgeCost(MyMesh::EdgeIter e_it)
{
	MyMesh::VertexIter v[]={e_it->vertex_iter(0),e_it->vertex_iter(1)};
	Vec3 edgeVec= v[0]->coordinate()- v[1]->coordinate();
	double edgeCost=0;
	for(unsigned k=0;k<2;k++)
	{
		Vec2 edgeVec2D = Vec2(edgeVec.dot(v[k]->direction(0)),edgeVec.dot(v[k]->direction(1)));
		edgeCost+= sqrt(pow(edgeVec2D.x,2)*v[k]->magnitude(0)+pow(edgeVec2D.y,2)*v[k]->magnitude(1));
	}	
	e_it->cost()=edgeCost/2.;
}
MyMesh::EdgeIter AnisGeodesic::worstEdgeInTriangle(MyMesh::FaceIter f_it)
{
	MyMesh::EdgeIter e_its[] = {f_it->edge_iter(0),f_it->edge_iter(1),f_it->edge_iter(2)};
	double& l0 = e_its[0]->cost();
	double& l1 = e_its[1]->cost();
	double& l2 = e_its[2]->cost();
	if(l0>=l1&&l0>=l2)
		return e_its[0];
	else if(l1>=l2)
		return e_its[1];
	else
		return e_its[2];
}

void AnisGeodesic::computeParameter()
{
	//computer scale parameter from global scale parameter; global scale enforces that all metric has same max anisotropy;
// 	if(!isDoGlobalScale || anisMesh==NULL) return;

	if(anisMetric==OUR_MAX || anisMetric==OUR_MIN)
	{
		scaleParameter = (globalScaleParameter - 1) / maxAnis;
	}

	if ((anisMetric == KOVACS_MAX || anisMetric == KOVACS_MIN) && anisMesh != NULL)
	{
		double b = globalScaleParameter*globalScaleParameter;
		double a,newa=1;
		for(MyMesh::VertexIter v_it=orgMesh->getVertices().begin(); v_it!=orgMesh->getVertices().end();v_it++)
		{
			double& k1=v_it->magnitude(0);
			double& k2=v_it->magnitude(1);

			a=(k1*k1-k2*k2*b)/(b-1);
			if(a>=0 && newa < a)
				newa = a;
		}
		scaleParameter=sqrt(newa);
	}

	if(anisMetric==CAMPEN_MAX || anisMetric==CAMPEN_MIN)
	{
		scaleParameter = globalScaleParameter;
	}
}

void AnisGeodesic::anisInit_Euclidean()
{
	if(anisMesh!=NULL) return;

	anisMesh=new MyMesh(orgMesh);

	anisMetric = EUCLIDEAN;
}
// void AnisGeodesic::anisInit_OurMin(unsigned mv)
// {
// 	if(anisMesh!=NULL) delete anisMesh;
// 	anisMesh=new MyMesh(orgMesh);
// 
// 	anisMetric = OUR_MIN;
// 
// 	computeParameter();
// 
// 	metricVariation=mv;
// 	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 	{
// 		Vec3 v1=v_it->direction(0);
// 		v_it->direction(0)=v_it->direction(1);
// 		v_it->direction(1)=v1;
// 	}
// 
// 	if(mv==1)
// 	{
// 		for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 		{
// 			double &k1=v_it->magnitude(0);
// 			double &k2=v_it->magnitude(1);
// 			double s=1+scaleParameter*fabs(k1-k2);
// 			s*=s;
// 			v_it->magnitude(0)=1;
// 			v_it->magnitude(1)=s;
// 		}
// 	}
// 	else if(mv==0)
// 	{
// 		for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 		{
// 			double &k1=v_it->magnitude(0);
// 			double &k2=v_it->magnitude(1);
// 			double s;			
// 			s=1+scaleParameter*fabs(k1-k2);
// 			s*=s;
// 			v_it->magnitude(0)=1/s;
// 			v_it->magnitude(1)=1;
// 		}
// 	}
// 	else if(mv==2)
// 	{
// 		for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 		{
// 			double &k1=v_it->magnitude(0);
// 			double &k2=v_it->magnitude(1);
// 			double s;			
// 			s=1+scaleParameter*fabs(k1-k2);
// 			v_it->magnitude(0)=1/s;
// 			v_it->magnitude(1)=s;
// 		}
// 	}
// }
// void AnisGeodesic::anisInit_OurMax(unsigned mv)
// {
// 	if(anisMesh!=NULL) delete anisMesh;
// 	anisMesh=new MyMesh(orgMesh);
// 
// 	anisMetric = OUR_MAX;
// 
// 	computeParameter();
// 
// 	metricVariation=mv;
// 	if(mv==1)
// 	{
// 		for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 		{
// 			double &k1=v_it->magnitude(0);
// 			double &k2=v_it->magnitude(1);
// 			double s=1+scaleParameter*fabs(k1-k2);
// 			s*=s;
// 			v_it->magnitude(0)=1;
// 			v_it->magnitude(1)=s;
// 		}
// 	}
// 	else if(mv==0)
// 	{
// 		for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 		{
// 			double &k1=v_it->magnitude(0);
// 			double &k2=v_it->magnitude(1);
// 			double s;			
// 			s=1+scaleParameter*fabs(k1-k2);
// 			s*=s;
// 			v_it->magnitude(0)=1/s;
// 			v_it->magnitude(1)=1;
// 		}
// 	}
// 	else if(mv==2)
// 	{
// 		for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++)
// 		{
// 			double &k1=v_it->magnitude(0);
// 			double &k2=v_it->magnitude(1);
// 			double s;			
// 			s=1+scaleParameter*fabs(k1-k2);
// 			v_it->magnitude(0)=1/s;
// 			v_it->magnitude(1)=s;
// 		}
// 	}
// }
void AnisGeodesic::anisInit_Our(MetricType mt, unsigned mv)
{
	if (anisMesh != NULL) delete anisMesh;
	anisMesh = new MyMesh(orgMesh);

	anisMetric = mt;

	computeParameter();

	metricVariation = mv;

	if (mt == OUR_MIN)
	{
		//reflect axis;
		for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
		{
			Vec3 v1 = v_it->direction(0);
			v_it->direction(0) = v_it->direction(1);
			v_it->direction(1) = v1;
		}
	}

	if (mv == 1)
	{
		for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
		{
			double k1 = v_it->magnitude(0);
			double k2 = v_it->magnitude(1);
			double s = 1 + scaleParameter*fabs(k1 - k2);
			s *= s;
			v_it->magnitude(0) = 1;
			v_it->magnitude(1) = s;
		}
	}
	else if (mv == 0)
	{
		for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
		{
			double &k1 = v_it->magnitude(0);
			double &k2 = v_it->magnitude(1);
			double s;
			s = 1 + scaleParameter*fabs(k1 - k2);
			s *= s;
			v_it->magnitude(0) = 1 / s;
			v_it->magnitude(1) = 1;
		}
	}
	else if (mv == 2)
	{
		for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
		{
			double &k1 = v_it->magnitude(0);
			double &k2 = v_it->magnitude(1);
			double s;
			s = 1 + scaleParameter*fabs(k1 - k2);
			v_it->magnitude(0) = 1 / s;
			v_it->magnitude(1) = s;
		}
	}

}
void AnisGeodesic::anisInit_OurMin(unsigned mv)
{
	anisInit_Our(OUR_MIN, mv);
}
void AnisGeodesic::anisInit_OurMax(unsigned mv)
{
	anisInit_Our(OUR_MAX, mv);

}
void AnisGeodesic::anisInit_Pott()
{
	//pottman metric is not presented by tensor.
	if(anisMesh!=NULL) delete anisMesh;
	anisMesh=new MyMesh(orgMesh);

	anisMetric = POTTMANN;
}
// void AnisGeodesic::anisInit_KovsMin()
// {
// 	if(anisMesh!=NULL) delete anisMesh;
// 	anisMesh=new MyMesh(orgMesh);
// 
// 	anisMetric = KOVACS_MIN;
// 
// 	computeParameter();
// 
// 	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++){
// 		Vec3 v1=v_it->direction(0);
// 		v_it->direction(0)=v_it->direction(1);
// 		v_it->direction(1)=v1;
// 	}
// 
// 	double a=scaleParameter*scaleParameter;
// 	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++){
// 		double k1=v_it->magnitude(0);
// 		double k2=v_it->magnitude(1);
// 		v_it->magnitude(0)=k2*k2+a;
// 		v_it->magnitude(1)=k1*k1+a;
// 	}
// }
// void AnisGeodesic::anisInit_KovsMax()
// {
// 	if(anisMesh!=NULL) delete anisMesh;
// 	anisMesh=new MyMesh(orgMesh);
// 
// 	anisMetric = KOVACS_MAX;
// 
// 	computeParameter();
// 
// 	double a=scaleParameter*scaleParameter;
// 	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++){
// 		double k1=v_it->magnitude(0);
// 		double k2=v_it->magnitude(1);
// 
// 		v_it->magnitude(0)=k2*k2+a;
// 		v_it->magnitude(1)=k1*k1+a;
// 	}
// }
void AnisGeodesic::anisInit_Kovs(MetricType mt)
{
	if (anisMesh != NULL) delete anisMesh;
	anisMesh = new MyMesh(orgMesh);

	anisMetric = mt;

	computeParameter();

	if (mt == KOVACS_MIN)
	{
		for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
		{
			Vec3 v1 = v_it->direction(0);
			v_it->direction(0) = v_it->direction(1);
			v_it->direction(1) = v1;
		}
	}

	double a = scaleParameter*scaleParameter;
	for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
	{
		double k1 = v_it->magnitude(0);
		double k2 = v_it->magnitude(1);
		v_it->magnitude(0) = k2*k2 + a;
		v_it->magnitude(1) = k1*k1 + a;
	}

}
void AnisGeodesic::anisInit_KovsMin()
{
	anisInit_Kovs(KOVACS_MIN);
}
void AnisGeodesic::anisInit_KovsMax()
{
	anisInit_Kovs(KOVACS_MAX);
}
// void AnisGeodesic::anisInit_CamMin()
// {
// 	if(anisMesh!=NULL) delete anisMesh;
// 	anisMesh=new MyMesh(orgMesh);
// 
// 	anisMetric = CAMPEN_MIN;
// 
// 	computeParameter();
// 
// 	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++){
// 		Vec3 v1=v_it->direction(0);
// 		v_it->direction(0)=v_it->direction(1);
// 		v_it->direction(1)=v1;
// 	}
// 
// 	double a=scaleParameter*scaleParameter;
// 	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin(); v_it!=anisMesh->getVertices().end();v_it++){
// 		v_it->magnitude(0)=1;
// 		v_it->magnitude(1)=a;
// 	}
// }
// void AnisGeodesic::anisInit_CamMax()
// {
// 	if (anisMesh != NULL) delete anisMesh;
// 	anisMesh = new MyMesh(orgMesh);
// 
// 	anisMetric = CAMPEN_MAX;
// 
// 	computeParameter();
// 
// 	double a = scaleParameter*scaleParameter;
// 	for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++){
// 		v_it->magnitude(0) = 1;
// 		v_it->magnitude(1) = a;
// 	}
// }
void AnisGeodesic::anisInit_Cam(MetricType mt)
{
	if (anisMesh != NULL) delete anisMesh;
	anisMesh = new MyMesh(orgMesh);

	anisMetric = mt;

	if (mt == CAMPEN_MIN)
	{
		for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++)
		{
			Vec3 v1 = v_it->direction(0);
			v_it->direction(0) = v_it->direction(1);
			v_it->direction(1) = v1;
		}
	}
	computeParameter();

	double a = scaleParameter*scaleParameter;
	for (MyMesh::VertexIter v_it = anisMesh->getVertices().begin(); v_it != anisMesh->getVertices().end(); v_it++){
		v_it->magnitude(0) = 1;
		v_it->magnitude(1) = a;
	}
}
void AnisGeodesic::anisInit_CamMin()
{
	anisInit_Cam(CAMPEN_MIN);
}
void AnisGeodesic::anisInit_CamMax()
{
	anisInit_Cam(CAMPEN_MAX);
}

void AnisGeodesic::meshEmbedding_Euclidean()
{
	for(MyMesh::EdgeIter e_it=anisMesh->getEdges().begin();e_it!=anisMesh->getEdges().end();e_it++)
	{
		MyMesh::VertexIter v[]={e_it->vertex_iter(0),e_it->vertex_iter(1)};
		Vec3 edgeVec= v[0]->coordinate()- v[1]->coordinate();
		e_it->cost()=edgeVec.length();
	}
	modifyInvalidTriangles();
}

void AnisGeodesic::meshEmbedding_Tensor(bool isModify)
{
	//tensor based anisotropy metric embedding;
	//Refer to papers: Dual loops meshing: Quality quad layouts on manifolds, Anisotropic quadrangulation, and Our paper;
	for(MyMesh::EdgeIter e_it=anisMesh->getEdges().begin();e_it!=anisMesh->getEdges().end();e_it++)
	{
		MyMesh::VertexIter v[]={e_it->vertex_iter(0),e_it->vertex_iter(1)};
		Vec3 edgeVec= v[0]->coordinate()- v[1]->coordinate();
		double edgeCost=0;
		for(unsigned k=0;k<2;k++)
		{
			Vec2 edgeVec2D = Vec2(edgeVec.dot(v[k]->direction(0)),edgeVec.dot(v[k]->direction(1)));
			edgeCost+= sqrt(pow(edgeVec2D.x,2)*v[k]->magnitude(0)+pow(edgeVec2D.y,2)*v[k]->magnitude(1));
		}
		e_it->cost()=edgeCost/2.;
	}

	if(isModify)
		modifyInvalidTriangles();
}
void AnisGeodesic::meshEmbedding_Pott()
{
	//pottman's(isophotic) metric is actually equal to kovacs, however not descriped by tensor. Practically, they are different.
	//Refer to paper: The isophotic metric and its application to feature sensitive morphology on surfaces
	double a = scaleParameter;
	for(MyMesh::EdgeIter e_it=anisMesh->getEdges().begin();e_it!=anisMesh->getEdges().end();e_it++)
	{
		MyMesh::VertexIter v[]={e_it->vertex_iter(0),e_it->vertex_iter(1)};
		double normalDifference;
		if(true)//directly normal difference
			normalDifference = computeAngle(v[0]->normal(),v[1]->normal());
		else
		{//theoretically more precise by doing normal transport. But find it is not neccessary;
			Vec3 edgeVec= v[0]->coordinate()- v[1]->coordinate();
			Vec3 edgeVecProj[]={projectVecToPlane(edgeVec,v[0]->normal()),projectVecToPlane(edgeVec,v[1]->normal())};
			Vec3 transNorm1 = basicNormalTransport(edgeVecProj[0],edgeVecProj[1],v[0]->normal());
			normalDifference= computeAngle(transNorm1,v[1]->normal());
		}
		e_it->cost() = sqrt(pow(e_it->length(),2)*a + normalDifference*normalDifference);
	}
	modifyInvalidTriangles();
}
void AnisGeodesic::meshEmbedding_Scribble(std::vector<unsigned>& verIds)
{
	//user interaction; user scribble on the mesh, and change the field and embedded mesh around the scribbles;
	if(verIds.empty()) return;

	if(anisMetric!=EUCLIDEAN && anisMetric!=POTTMANN)
	{
		unsigned bigVal = 100000000;
		std::vector<unsigned> edgeIds(1, bigVal);
		for(unsigned i=0;i<verIds.size();i++)
		{
			auto v_it = anisMesh->getVIter()[verIds[i]];
			for(unsigned j=0;j<v_it->edge_iter().size();j++)
			{
				unsigned eid=v_it->edge_iter()[j]->id();
				if(std::find(edgeIds.begin(),edgeIds.end(),eid)==edgeIds.end()) edgeIds.push_back(eid);
			}			
		}

		edgeIds.push_back(0);
		std::sort(edgeIds.begin(),edgeIds.end()); 
		edgeIds.pop_back();
		MyMesh::EdgeIter e_it=anisMesh->getEdges().begin();

		for(unsigned i=1; i<edgeIds.size();i++)
		{
			std::advance(e_it,edgeIds[i]-edgeIds[i-1]);
			MyMesh::VertexIter v[]={e_it->vertex_iter(0),e_it->vertex_iter(1)};
			Vec3 edgeVec= v[0]->coordinate()- v[1]->coordinate();
			double edgeCost=0;
			for(unsigned k=0;k<2;k++)
			{
				Vec2 edgeVec2D = Vec2(edgeVec.dot(v[k]->direction(0)),edgeVec.dot(v[k]->direction(1)));
				edgeCost+= sqrt(pow(edgeVec2D.x,2)*v[k]->magnitude(0)+pow(edgeVec2D.y,2)*v[k]->magnitude(1));
			}
			e_it->cost()=edgeCost/2.;
		}

		modifyInvalidTriangles();
	}
}

bool AnisGeodesic::modifyInvalidTriangles()
{
	//find all invalid triangles, and subdivided along the longest edge;
	//by 1. insert a mid point on that edge, (compute the tensor on the new vertex)
	//   2. and connect to it's opposite vertices on the two opposite triangles,
	//   3. then update datas;

	std::priority_queue<FaceCost> badFaceSet;
	std::vector<MyMesh::FaceIter> faceIts;
	unsigned prioQueueType=1;
	for(auto f_it=anisMesh->getFaces().begin(); f_it!=anisMesh->getFaces().end();f_it++)
	{
		if(f_it->triangleInequality()==false)
		{
			f_it->triangleCost(prioQueueType);
			badFaceSet.push(FaceCost(f_it->id(),f_it->cost()));
		}
		faceIts.push_back(f_it);
	}
	faceIts.reserve(anisMesh->getFaces().size()+10*badFaceSet.size());

	std::vector<int> &vertexOnEdge = meshMap.vertexOnEdge; 
	std::vector<int> &edgeAccestor = meshMap.edgeAccestor;
	if(anisMesh->getVertices().size()==orgMesh->getVertices().size())
	{
		vertexOnEdge.resize(anisMesh->getVertices().size(),-1);
		edgeAccestor.resize(anisMesh->getEdges().size());
		for(unsigned i=0;i<edgeAccestor.size();i++) edgeAccestor[i]=i; //if <0, means inside triangle;
	}

	std::list<MyMesh::Edge> emptyEdges(1);
	std::list<MyMesh::Face> emptyFaces(1);
	unsigned badInit=badFaceSet.size();
	unsigned badLowest=badInit;
	unsigned badHighest=badInit;
	unsigned times=0;
	unsigned faceSize=anisMesh->getFaces().size();
	while(!badFaceSet.empty())
	{

		FaceCost fc =badFaceSet.top(); badFaceSet.pop();
		auto f1= faceIts[fc.id()];
		if(fc.cost()!=f1->cost()) continue;

		if(badLowest>badFaceSet.size()) badLowest=badFaceSet.size();
		if(badHighest<badFaceSet.size()) badHighest=badFaceSet.size();
		if(times>=subTimesMax || fc.cost()<worstTriangleCost)	break;
		times++;

		//find edges,faces,vertices that will be used;
		MyMesh::EdgeIter e0=worstEdgeInTriangle(f1);
		MyMesh::FaceIter f2;
		MyMesh::VertexIter v1=f1->opposite_vertex(e0);
		MyMesh::VertexIter v2;
		MyMesh::VertexIter v3=e0->vertex_iter(0);
		MyMesh::VertexIter v4=e0->vertex_iter(1);
		MyMesh::EdgeIter e13=f1->opposite_edge(v4);
		MyMesh::EdgeIter e14=f1->opposite_edge(v3);
		MyMesh::EdgeIter e23;
		MyMesh::EdgeIter e24;

		auto fop13 = e13->opposite_face(f1);
		auto fop14 = e14->opposite_face(f1);
		MyMesh::FaceIter fop23;
		MyMesh::FaceIter fop24;

		bool isManifold = e0->manifold();
		if(isManifold==true){
			f2= e0->opposite_face(f1);
			v2=f2->opposite_vertex(e0);
			e23=f2->opposite_edge(v4);
			e24=f2->opposite_edge(v3);
			fop23 = e23->opposite_face(f2);
			fop24 = e24->opposite_face(f2);
		}
		else{// a boundary edge;

		}
/*
		how to do subdivision?
		1. add midpoint on the bad edge (add a vertex, two edges, and remove one edge)
		2. connect midpoint to two opposite vertices (add two edges, four faces, and remove two faces)
		3. update references
*/
		//1. create and add a new vertex
		MyMesh::Vertex newVer;
		creatMidPoint(newVer,e0);
		if(anisMetric==OUR_MIN||anisMetric==OUR_MAX)
		{ //tensor need to be normalized
			if(metricVariation==1)
			{
				newVer.magnitude(1) /=newVer.magnitude(0); 
				newVer.magnitude(0)=1; 
			}
			else if(metricVariation==0)
			{
				newVer.magnitude(0) /=newVer.magnitude(1);
				newVer.magnitude(1)=1; 
			}
			else if(metricVariation==2)
			{
				double ta= sqrt(newVer.magnitude(0)*newVer.magnitude(1));
				newVer.magnitude(0) /=ta;
				newVer.magnitude(1) /=ta; 
			}
		}
		MyMesh::VertexIter v0;
		newVer.id()=anisMesh->getVertices().size();
		anisMesh->getVertices().push_back(newVer);
		v0=anisMesh->getVertices().end(); v0--; *v0=newVer;

		//this will be used on tracing geodesy;
		if(edgeAccestor[e0->id()]>=0) 
			vertexOnEdge.push_back(edgeAccestor[e0->id()]);
		else
			vertexOnEdge.push_back(-1);

		//2. create and add 4 new edges;
		MyMesh::Edge newe01,newe02,newe03,newe04;

		newe01.vertex_iter(0)=v0;newe01.vertex_iter(1)=v1;
		newe01.face_iter(0)=emptyFaces.begin();newe01.face_iter(1)=emptyFaces.begin(); //init adjface iter, emptyFaces means nothing
		newe01.length()=(v0->coordinate()-v1->coordinate()).length();

		newe03.vertex_iter(0)=v0;newe03.vertex_iter(1)=v3;
		newe03.face_iter(0)=emptyFaces.begin();newe03.face_iter(1)=emptyFaces.begin();
		newe03.length()=(v0->coordinate()-v3->coordinate()).length();

		newe04.vertex_iter(0)=v0;newe04.vertex_iter(1)=v4;
		newe04.face_iter(0)=emptyFaces.begin();newe04.face_iter(1)=emptyFaces.begin();
		newe04.length()=(v0->coordinate()-v4->coordinate()).length();

		if(isManifold==true)
		{
			newe02.vertex_iter(0)=v0;newe02.vertex_iter(1)=v2;
			newe02.face_iter(0)=emptyFaces.begin();newe02.face_iter(1)=emptyFaces.begin();
			newe02.length()=(v0->coordinate()-v2->coordinate()).length();

			newe01.manifold()=true;newe02.manifold()=true;newe03.manifold()=true;newe04.manifold()=true;
		}
		else
		{
			newe01.manifold()=true;newe03.manifold()=false;newe04.manifold()=false;
		}

		MyMesh::EdgeIter e01,e02,e03,e04;
		newe01.id()=e0->id();*e0=newe01;
		e01=e0;
		newe03.id()=anisMesh->getEdges().size();	anisMesh->getEdges().push_back(newe03);
		e03=anisMesh->getEdges().end(); e03--; *e03=newe03;
		newe04.id()=anisMesh->getEdges().size();	anisMesh->getEdges().push_back(newe04);
		e04=anisMesh->getEdges().end(); e04--; *e04=newe04;
		computeEdgeCost(e01);
		computeEdgeCost(e03);
		computeEdgeCost(e04);

		//update datas that use the new edges;
		edgeAccestor.push_back(-1);		
		edgeAccestor.push_back(-1);

		if(isManifold==true)
		{
			newe02.id()=anisMesh->getEdges().size();	
			anisMesh->getEdges().push_back(newe02);
			e02=anisMesh->getEdges().end(); e02--; *e02=newe02;
			edgeAccestor.push_back(-1);
			computeEdgeCost(e02);
		}

		if(vertexOnEdge.back()>=0)
		{
			edgeAccestor[e03->id()]=edgeAccestor[e0->id()];
			edgeAccestor[e04->id()]=edgeAccestor[e0->id()];
		}

		edgeAccestor[e0->id()]=-1;

		v0->vertex_iter().push_back(v1); v0->edge_iter().push_back(e01);
		v0->vertex_iter().push_back(v3); v0->edge_iter().push_back(e03);
		v0->vertex_iter().push_back(v4); v0->edge_iter().push_back(e04);
		v1->vertex_iter().push_back(v0); v1->edge_iter().push_back(e01);
		v3->vertex_iter().push_back(v0); v3->edge_iter().push_back(e03);
		v4->vertex_iter().push_back(v0); v4->edge_iter().push_back(e04);
		if(isManifold==true)
		{
			v0->vertex_iter().push_back(v2); v0->edge_iter().push_back(e02);
			v2->vertex_iter().push_back(v0); v2->edge_iter().push_back(e02);
		}

		//3. create and add 4 new faces;
		MyMesh::Face newf13,newf14,newf23,newf24;

		newf13.vertex_iter(0)=v0;newf13.vertex_iter(1)=v1;newf13.vertex_iter(2)=v3;
		if(fop13->next_vertex(v1)->id()==v3->id()){newf13.vertex_iter(1)=v3;newf13.vertex_iter(2)=v1;} //consist orientation
		newf13.center_point()=(v0->coordinate()+v1->coordinate()+v3->coordinate())/3.;
		newf13.edge_iter(0)=e01;newf13.edge_iter(1)=e03;newf13.edge_iter(2)=e13;

		newf14.vertex_iter(0)=v0;newf14.vertex_iter(1)=v1;newf14.vertex_iter(2)=v4;
		if(fop14->next_vertex(v1)->id()==v4->id()){newf14.vertex_iter(1)=v4;newf14.vertex_iter(2)=v1;}
		newf14.center_point()=(v0->coordinate()+v1->coordinate()+v4->coordinate())/3.;
		newf14.edge_iter(0)=e01;newf14.edge_iter(1)=e04;newf14.edge_iter(2)=e14;

		if(isManifold==true)
		{
			newf23.vertex_iter(0)=v0;newf23.vertex_iter(1)=v2;newf23.vertex_iter(2)=v3;
			if(fop23->next_vertex(v2)->id()==v3->id()){newf23.vertex_iter(1)=v3;newf23.vertex_iter(2)=v2;}
			newf23.center_point()=(v0->coordinate()+v2->coordinate()+v3->coordinate())/3.;
			newf23.edge_iter(0)=e02;newf23.edge_iter(1)=e03;newf23.edge_iter(2)=e23;

			newf24.vertex_iter(0)=v0;newf24.vertex_iter(1)=v2;newf24.vertex_iter(2)=v4;
			if(fop24->next_vertex(v2)->id()==v4->id()){newf24.vertex_iter(1)=v4;newf24.vertex_iter(2)=v2;}
			newf24.center_point()=(v0->coordinate()+v2->coordinate()+v4->coordinate())/3.;
			newf24.edge_iter(0)=e02;newf24.edge_iter(1)=e04;newf24.edge_iter(2)=e24;
		}

		MyMesh::FaceIter f13,f14,f23,f24;
		newf13.id()=f1->id();	*f1=newf13;
		f13=f1;
		if(isManifold==true)
		{
			newf23.id()=f2->id();	*f2=newf23;
			f23=f2;
		}
		newf14.id()=anisMesh->getFaces().size();	anisMesh->getFaces().push_back(newf14);
		f14=anisMesh->getFaces().end();f14--; *f14=newf14;
		if(isManifold==true)
		{
			newf24.id()=anisMesh->getFaces().size();	anisMesh->getFaces().push_back(newf24);
			f24=anisMesh->getFaces().end();f24--; *f24=newf24;
		}
		faceIts.push_back(f14);
		if(isManifold==true)faceIts.push_back(f24);

		//update adjacent links, for old edges and news';
		if(isManifold==true)
		{
			e01->face_iter(0)=f13;	e01->face_iter(1)=f14;
			e02->face_iter(0)=f23;	e02->face_iter(1)=f24;
			e03->face_iter(0)=f13;	e03->face_iter(1)=f23;
			e04->face_iter(0)=f14;	e04->face_iter(1)=f24;

			e14->face_iter(0)==f1? e14->face_iter(0)=f14: e14->face_iter(1)=f14;
			e24->face_iter(0)==f2? e24->face_iter(0)=f24: e24->face_iter(1)=f24;
		}
		else
		{
			e01->face_iter(0)=f13;	e01->face_iter(1)=f14;
			e03->face_iter(0)=f13;
			e04->face_iter(0)=f14;

			e14->face_iter(0)==f1? e14->face_iter(0)=f14: e14->face_iter(1)=f14;
		}


		//check validation for new faces
		MyMesh::FaceIter f_it[4]={f13,f14,f23,f24};
		unsigned num=4;
		if(isManifold==false) num=2;
		for(unsigned i=0;i<num;i++)
		{
			if(!f_it[i]->triangleInequality())
			{
				f_it[i]->triangleCost(prioQueueType);
				badFaceSet.push(FaceCost(f_it[i]->id(),f_it[i]->cost()));
			}
		}
	}

	badFaces.clear();
	while(!badFaceSet.empty())
	{
		FaceCost fcst = badFaceSet.top();
		badFaces.push_back(fcst.id());
		badFaceSet.pop();
	}
	if(badFaces.size()!=0) std::sort(badFaces.begin(),badFaces.end());

	faceIts.clear();

	return true;
}

bool mySort1(std::pair<unsigned,Vec3> &a,std::pair<unsigned,Vec3> &b)
{
	return a.first<b.first;
}
void AnisGeodesic::modifyAnisField(std::vector<std::pair<unsigned, Vec3> >& newDirection)
{
	if (newDirection.empty()) return;

	std::sort(newDirection.begin(),newDirection.end(),mySort1);

	std::vector<unsigned> frontVertices;
	std::vector<unsigned> usedVertices;
	std::vector<unsigned> newVertices(1,100000000);
	std::vector<std::vector<unsigned> > adjDirection(1);
	auto f_it = anisMesh->getFaces().begin();
	for(unsigned i=0;i<newDirection.size();i++)
	{
		if(i==0)std::advance(f_it,newDirection[i].first);
		else std::advance(f_it,newDirection[i].first-newDirection[i-1].first);
		for(unsigned j=0;j<3;j++)
		{
			unsigned ind = std::find(newVertices.begin(),newVertices.end(),f_it->vertex_iter(j)->id())- newVertices.begin();
			if(ind>=newVertices.size())
			{
				newVertices.push_back(f_it->vertex_iter(j)->id()); adjDirection.push_back(std::vector<unsigned>(1,i));
			}
			else
				adjDirection[ind].push_back(i);
		}
	}
	newVertices.erase(newVertices.begin()); adjDirection.erase(adjDirection.begin());

	double  newMag1,newMag2;
	double s = 1 + modifyFieldRatio*globalScaleParameter;
	if(anisMetric==OUR_MAX || anisMetric==OUR_MIN)
	{
		if(metricVariation==1)
		{
			newMag1 = 1; newMag2 = s*s;
		}
		else if(metricVariation==0)
		{
			newMag1 = 1/(s*s); newMag2 = 1;
		}
		else if(metricVariation==2)
		{
			newMag1 = 1/s; newMag2 = s;
		}
	}
	else if(anisMetric==KOVACS_MAX || anisMetric==KOVACS_MIN)
	{
		newMag1 = 1; newMag2 = s*s;
	}
	else if(anisMetric==CAMPEN_MAX || anisMetric==CAMPEN_MIN)
	{
		newMag1 = 1; newMag2 = s*s;
	}

	cout << "sketched tensor(min:" << newMag1 << ", max:" << newMag2 << ").\n";
	for(unsigned i=0;i<newVertices.size();i++)
	{
		auto v_it = anisMesh->getVIter()[newVertices[i]];
		unsigned num=adjDirection[i].size();
		std::list<MyMesh::Vertex> sv(num);
		std::vector<MyMesh::VertexIter> sv_it(num);
		MyMesh::VertexIter tv = sv.begin();
		for(unsigned j=0;j<adjDirection[i].size();j++)
		{
			tv->normal()=v_it->normal();
			tv->direction(0)= projectVecToPlane(newDirection[adjDirection[i][j]].second,tv->normal());
			tv->direction(1)=tv->normal().cross(tv->direction(0));
			tv->magnitude(0)=newMag1;
			tv->magnitude(1)=newMag2;
			sv_it[j]=tv; tv++;
		}
		newPointDirection(v_it,sv_it,1.);
		{
			if (metricVariation == 1)
			{
				v_it->magnitude(0) = 1; 
				v_it->magnitude(1) = v_it->magnitude(1) / v_it->magnitude(0);
			}
			else if (metricVariation == 0)
			{
				v_it->magnitude(0) = v_it->magnitude(0) / v_it->magnitude(1);
				v_it->magnitude(1) = 1;
			}
			else if (metricVariation == 2)
			{
				s = sqrt(v_it->magnitude(0) / v_it->magnitude(1));
				v_it->magnitude(0) = 1 / s;
				v_it->magnitude(1) = s;
			}
		}
		usedVertices.push_back(newVertices[i]);
	}
	frontVertices=usedVertices; 
	//propagation, find new vertices to be assign with new tensor;
	unsigned propagationDepth = 2;
	//if(modifyNeighboreSize>10) modifyNeighboreSize=10;
	while (propagationDepth <=modifyNeighboreSize)
	{
		//new vertices from front vertices; we used gaussian regresson
		double wt=exp( pow(double(propagationDepth-1)/2.,2)/-2.);
		//double wt=1-pow( double(propagationDepth-1)/5 ,2);
		//double wt=double(modifyNeighboreSize-double(propagationDepth+1)/2)/double(modifyNeighboreSize);

		newVertices.clear();newVertices.push_back(100000000);//a big unsign number;
		adjDirection.clear();adjDirection.push_back(std::vector<unsigned>(0));
		for(unsigned i=0;i<frontVertices.size();i++)
		{
			auto v_it = anisMesh->getVIter()[frontVertices[i]];

			for(unsigned j=0;j<v_it->vertex_iter().size();j++)
			{
				unsigned vid = v_it->vertex_iter()[j]->id();
				if(std::find(usedVertices.begin(),usedVertices.end(),vid)!=usedVertices.end()) continue;

				unsigned newId = std::find(newVertices.begin(),newVertices.end(),vid) - newVertices.begin();
				if(newId>=newVertices.size())
				{
					newVertices.push_back(vid);  
					adjDirection.push_back(std::vector<unsigned>(1,v_it->id()));
				}
				else
				{
					adjDirection[newId].push_back(v_it->id());
				}
			}
		}
		newVertices.erase(newVertices.begin()); adjDirection.erase(adjDirection.begin());

		for(unsigned i=0;i<newVertices.size();i++)
		{
			auto v_it = anisMesh->getVIter()[newVertices[i]];

			std::vector<MyMesh::VertexIter> sv_it;
			for(unsigned j=0;j<adjDirection[i].size();j++)
			{
				auto tv_it = anisMesh->getVIter()[adjDirection[i][j]];
				sv_it.push_back(tv_it);
			}
			newPointDirection(v_it,sv_it,wt);
			usedVertices.push_back(newVertices[i]);
		}
		frontVertices=newVertices;

		propagationDepth++;
	}

	// re-embed mesh;
	meshEmbedding_Scribble(usedVertices);
}

void AnisGeodesic::updateGeoMesh_EdgeLength()
{
/*
	std::vector<std::vector<double> > te;
	for(MyMesh::EdgeIter e_it=theMesh->getEdges().begin();e_it!=theMesh->getEdges().end();e_it++){
		std::vector<double> a(2);
		a[0]=min(e_it->vertex_iter(0)->id(),e_it->vertex_iter(1)->id());
		a[1]=max(e_it->vertex_iter(0)->id(),e_it->vertex_iter(1)->id());
		te.push_back(a);
	}
	std::sort(te.begin(),te.end());

	std::vector<std::vector<double> > te2;
	unsigned dif=0;
	for(unsigned i=0;i<te.size();i++){
		std::vector<double> a;
		a.push_back(meshForAll->edges()[i].adjacent_vertices()[0]->id());
		a.push_back(meshForAll->edges()[i].adjacent_vertices()[1]->id());
		te2.push_back(a);
	}
	std::sort(te2.begin(),te2.end());
	for(unsigned i=0;i<te.size();i++){
		if(te2[i]!=te[i])
			dif++;
	}
	te.clear();
*/
	std::vector<std::vector<double> > te;
	unsigned ind=0;
	for(MyMesh::EdgeIter e_it=anisMesh->getEdges().begin();e_it!=anisMesh->getEdges().end();e_it++)
	{
		std::vector<double> a(4);
		a[0]=min(e_it->vertex_iter(0)->id(),e_it->vertex_iter(1)->id());
		a[1]=max(e_it->vertex_iter(0)->id(),e_it->vertex_iter(1)->id());
		a[2]=e_it->cost();
		a[3]=ind;
		te.push_back(a);
		ind++;
	}
	std::sort(te.begin(),te.end());

	std::vector<unsigned> edgeMap(geoMesh->edges().size());
	ind=0;
	for(unsigned i=0;i<te.size();i++)
	{
		unsigned &id1=geoMesh->edges()[ind].adjacent_vertices()[0]->id();
		unsigned &id2=geoMesh->edges()[ind].adjacent_vertices()[1]->id();
		if(id1==te[i][0] && id2==te[i][1])
		{
			geoMesh->edges()[ind].length() = te[i][2];
			edgeMap[ind]=unsigned(te[i][3]); ind++;
		}
	}

	// update corner angle and saddle points...
	//compute angles for the faces
	unsigned errCount=0;
	std::vector<unsigned> errTri;
	unsigned faceNum=geoMesh->faces().size();
	for(unsigned i=0; i<faceNum; ++i)
	{
		geodesic::Face& f = geoMesh->faces()[i];
		double abc[3];		
		double sum = 0;
		for(unsigned j=0; j<3; ++j)		//compute angle adjacent to the vertex j
		{
			for(unsigned k=0; k<3; ++k)
			{
				geodesic::vertex_pointer v = f.adjacent_vertices()[(j + k)%3];
				abc[k] = f.opposite_edge(v)->length();
			}
			double angle = geodesic::angle_from_edges(abc[0], abc[1], abc[2]);

			if(!((abc[0]+abc[1])>=abc[2] && (fabs(abc[0]-abc[1])<=abc[2]))&&j==0)
			{
				errCount++;errTri.push_back(i);
			}
			//assert(angle>1e-5);						//algorithm works well with non-degenerate meshes only 
			f.corner_angles()[j] = angle;
			sum += angle;
		}
		//assert(std::abs(sum - M_PI) < 1e-5);		//algorithm works well with non-degenerate meshes only 
	}

	//define m_turn_around_flag for vertices
	std::vector<double> total_vertex_angle(geoMesh->vertices().size());
	for(unsigned i=0; i<geoMesh->faces().size(); ++i)
	{
		geodesic::Face& f = geoMesh->faces()[i];
		for(unsigned j=0; j<3; ++j)
		{
			geodesic::vertex_pointer v = f.adjacent_vertices()[j];
			total_vertex_angle[v->id()] += f.corner_angles()[j];
		}
	}

	for(unsigned i=0; i<geoMesh->vertices().size(); ++i)
	{
		geodesic::Vertex& v = geoMesh->vertices()[i];
		v.saddle_or_boundary() = (total_vertex_angle[v.id()] > 2.0*M_PI - 1e-5); 
	}

	for(unsigned i=0; i<geoMesh->edges().size(); ++i)
	{
		geodesic::Edge& e = geoMesh->edges()[i];
		if(e.is_boundary())
		{
			e.adjacent_vertices()[0]->saddle_or_boundary() = true;
			e.adjacent_vertices()[1]->saddle_or_boundary() = true;
		}
	}

	geoAlgo = new  geodesic::GeodesicAlgorithmExact(geoMesh);
	meshMap.edgeMap=edgeMap;

	updateGeoMesh_EdgeList();
}
void AnisGeodesic::updateGeoMesh_EdgeList()
{
	std::vector<unsigned>& edgeMap = meshMap.edgeMap;
	std::vector<int>& edgeAccestor = meshMap.edgeAccestor;

	std::vector<std::vector<unsigned> > &edgeList = meshMap.edgeList;
	edgeList.clear(); edgeList.resize(geoMesh->edges().size());
	for(unsigned i=0;i<geoMesh->edges().size();i++)
	{
		int veid = edgeAccestor[edgeMap[i]];// geoMesh -> ourMetricMesh -> originalMesh
		if(veid<0) continue;
		edgeList[veid].push_back(i);
	}
}

void AnisGeodesic::geodesyInit() 
{
	std::vector<double> tpoints;
	for(MyMesh::VertexIter v_it=anisMesh->getVertices().begin();v_it!=anisMesh->getVertices().end();v_it++)
	{
		tpoints.push_back(v_it->coordinate().x);
		tpoints.push_back(v_it->coordinate().y);
		tpoints.push_back(v_it->coordinate().z);
	}
	std::vector<unsigned> tfaces;
	for(auto f_it=anisMesh->getFaces().begin();f_it!=anisMesh->getFaces().end();f_it++)
	{
		if(f_it->valid_triangle()==false) continue;
		tfaces.push_back(f_it->vertex_iter(0)->id());
		tfaces.push_back(f_it->vertex_iter(1)->id());
		tfaces.push_back(f_it->vertex_iter(2)->id());
	}

	if(geoMesh!=NULL) delete geoMesh;
	geoMesh = new geodesic::Mesh;
	geoMesh->initialize_mesh_data(tpoints, tfaces);		//create internal mesh data structure including edges
	if (geoAlgo != NULL) delete geoAlgo;
	geoAlgo=new  geodesic::GeodesicAlgorithmExact(geoMesh);	//create exact algorithm for the mesh

	updateGeoMesh_EdgeLength();
}
void AnisGeodesic::geodesyPropagation(const IntersectionType stype, const unsigned source_index, const Vec3 &scoord)
{
	if(geoMesh==NULL || geoAlgo==NULL) return;
	if(stype==VERTEX)
	{//vertex
		if(source_index>=geoMesh->vertices().size()) return;
		geodesic::SurfacePoint source(&geoMesh->vertices()[source_index]);		//create source 
		std::vector< geodesic::SurfacePoint> sources(1,source);
		geoAlgo->propagate(sources);
	}
	else if(stype==EDGE)
	{//edge
		if(source_index>=geoMesh->edges().size()) return;
		unsigned sid = getEdgeId(source_index, scoord);
		geodesic::SurfacePoint source(&geoMesh->edges()[sid], scoord.x, scoord.y, scoord.z);
		std::vector< geodesic::SurfacePoint> sources(1,source);
		geoAlgo->propagate(sources);
	}
}	
unsigned AnisGeodesic::getEdgeId(const unsigned src,const Vec3& srcp)
{
	std::vector<unsigned> &edgelt = meshMap.edgeList[src];

	double tdis=100000.; unsigned ti=0;
	geodesic::Point3D* ps=new geodesic::Point3D; ps->set(srcp.x,srcp.y,srcp.z);
	for(unsigned i=0;i<edgelt.size();i++)
	{
		geodesic::vertex_pointer v1=geoMesh->edges()[edgelt[i]].v0();
		geodesic::vertex_pointer v2=geoMesh->edges()[edgelt[i]].v1();
		double td=v1->distance(ps);
		if(tdis>td){
			tdis = td; 
			ti=edgelt[i];
		}
		td=v2->distance(ps);
		if(tdis>td){ 
			tdis = td; 
			ti=edgelt[i];
		}
	}
	return ti;
}
void AnisGeodesic::geodesyTraceBack(const IntersectionType ttype, const unsigned targe_index, const Vec3 &tcoord)
{
	if (geoMesh == NULL || geoAlgo == NULL) return;
	geoPath.clear();	geoPathVE.clear();
	std::vector< geodesic::SurfacePoint> path;	//geodesic path is a sequence of SurfacePoints
	if(ttype==VERTEX)
	{//vertex
		if(targe_index>=geoMesh->vertices().size()) return;
		geodesic::SurfacePoint target(&geoMesh->vertices()[targe_index]);		//create source 
		geoPathCost=geoAlgo->trace_back(target, path);
	}
	else if(ttype==EDGE)
	{//edge
		if(targe_index>=geoMesh->edges().size()) return;
		unsigned tid=getEdgeId(targe_index,tcoord);
		geodesic::SurfacePoint target(&geoMesh->edges()[tid],tcoord.x,tcoord.y,tcoord.z);
		geoPathCost=geoAlgo->trace_back(target, path);
	}

	for(unsigned i = 0; i<path.size(); ++i)
	{
		geodesic::SurfacePoint& s = path[i];
		int veid; 
		IntersectionType veType;
		if(s.type() == geodesic::VERTEX)
		{
			veType=VERTEX;
			geodesic::vertex_pointer v = static_cast<geodesic::vertex_pointer>(s.base_element());	
			if(v->id()<geoMesh->vertices().size())// vertex is actual the vertex of original mesh;
			{
				veid = v->id();
			}
			else
			{
				veid = meshMap.vertexOnEdge[v->id()];
				if(veid<0)//vertex is inside triangle of original mesh
					continue;
				else
				{ // vertex is on the edge of original mesh
					veType=EDGE;
				}
			}
		}
		if(s.type() == geodesic::EDGE)
		{
			veType=EDGE;
			geodesic::edge_pointer e = static_cast<geodesic::edge_pointer>(s.base_element());	
			veid = meshMap.edgeAccestor[meshMap.edgeMap[e->id()]];// geoMesh -> ourMetricMesh -> originalMesh
			if(veid<0)
				continue;
		}
		geoPath.push_back(Vec3(s.x(), s.y(), s.z()));
		geoPathVE.push_back(std::pair<IntersectionType,unsigned>(veType,unsigned(veid)));
	}
	reverse(geoPath.begin(),geoPath.end());
	reverse(geoPathVE.begin(),geoPathVE.end());

	if(geoPath.size()<=1 || geoPathCost<=0)
	{
		geoPath.clear(); geoPathCost=0.;
	}
}
