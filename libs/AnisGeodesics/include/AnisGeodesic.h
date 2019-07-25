/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#ifndef ANIS_GEODESIC_H
#define ANIS_GEODESIC_H

#include "MyMesh.h"
#include "amlVec.h"

enum MetricType
{
	EUCLIDEAN=0,
	OUR_MAX,OUR_MIN,
	POTTMANN,
	KOVACS_MAX,KOVACS_MIN,
	CAMPEN_MAX,CAMPEN_MIN,
	ALL
};

namespace geodesic{
	class Mesh;
	class GeodesicAlgorithmExact;
};

class AnisGeodesic{
public: 

	class FaceCost{ //bad triangle priority queue
	public: FaceCost(unsigned id, double cost): m_id(id),m_cost(cost){}
			~FaceCost(){}
			bool operator<(const FaceCost& other) const{return m_cost<other.m_cost;}
			bool operator>(const FaceCost& other) const{return m_cost>other.m_cost;}
			unsigned& id(){return m_id;}
			double& cost(){return m_cost;}
	private:
		unsigned m_id;
		double m_cost;
	};

	enum IntersectionType
	{
		VERTEX,
		EDGE,
	};

	struct MeshToMeshMap{ //original mesh <-> embedded mesh <-> mesh for anis geodesy
		std::vector<unsigned> edgeMap;
		std::vector<int> vertexOnEdge;
		std::vector<int> edgeAccestor;

		std::vector<std::vector<unsigned> > edgeList;
	};

	AnisGeodesic(){ init(); }
	~AnisGeodesic(){}

	void init(){orgMesh=anisMesh=NULL; geoMesh=NULL; geoAlgo=NULL; anisMetric=EUCLIDEAN;scaleParameter=1.; metricVariation=0;}

	void computeParameter();

	void anisInit_Euclidean();
	void anisInit_OurMin(unsigned mv=0); //our Metric has 3 Variations(0,1,2), 0 - mags>=1; 1 - mags=<1; 2 - mag1<1 & mag2>1 
	void anisInit_OurMax(unsigned mv=0);
	void anisInit_Pott();
	void anisInit_KovsMin();
	void anisInit_KovsMax();
	void anisInit_CamMin();
	void anisInit_CamMax();

	void meshEmbedding_Euclidean();
	void meshEmbedding_Tensor(bool isModify=true); //edges of mesh are rescaled by anis field;
	void meshEmbedding_Pott();
	void meshEmbedding_Scribble(std::vector<unsigned>& verIds);

	void modifyAnisField(std::vector<std::pair<unsigned, Vec3> >& newDirection);

	void updateGeoMesh_EdgeLength();
	void updateGeoMesh_EdgeList();
	void geodesyInit();
	void geodesyPropagation(const IntersectionType stype, const unsigned source_index, const Vec3 &scoord = Vec3(0, 0, 0));
	void geodesyTraceBack(const IntersectionType ttype, const unsigned targe_index, const Vec3 &tcoord = Vec3(0, 0, 0));

public:
	MyMesh *orgMesh,*anisMesh;

	MetricType anisMetric;
	unsigned metricVariation;
	double scaleParameter;
	static bool isDoGlobalScale;
	static double globalScaleParameter;
	static double maxAnis;

	std::vector<unsigned> badFaces;
	static unsigned subTimesMax;
	static double worstTriangleCost;

	geodesic::Mesh *geoMesh;
	geodesic::GeodesicAlgorithmExact *geoAlgo;
	MeshToMeshMap meshMap; //map between orgMesh and geoMesh;

	unsigned seedPoint,freePoint;
	double geoPathCost;
	std::vector<Vec3> geoPath;
	std::vector<std::pair<IntersectionType, unsigned> > geoPathVE; //every point on the path lies on the Vertices or Edges of the mesh;

	//update an exist anis-embedded mesh. This takes in scribbles,which may affect in a certain region;
	static unsigned modifyNeighboreSize;
	static double modifyFieldRatio;


private:

	void anisInit_Our(MetricType mt, unsigned mv);
	void anisInit_Kovs(MetricType mt);
	void anisInit_Cam(MetricType mt);

	Vec3 projectVecToPlane(const Vec3 &v, const Vec3 &n);
	double computeAngle(const Vec3 &org1,const Vec3 &org2);

	//scale mesh;
	void creatMidPoint(MyMesh::Vertex& v, MyMesh::EdgeIter& e_it);
	void newPointDirection(MyMesh::VertexIter& v_it, std::vector<MyMesh::VertexIter>& sv_it, double rat);
	
	void computeEdgeCost(MyMesh::EdgeIter e_it);
	MyMesh::EdgeIter worstEdgeInTriangle(MyMesh::FaceIter f_it);
	bool checkInvalidTriangle(MyMesh*inMesh);

	bool modifyInvalidTriangles();

	unsigned getEdgeId(const unsigned src,const Vec3& srcp);

};

#endif
