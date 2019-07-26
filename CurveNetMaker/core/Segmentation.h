#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "AnisGeodesic.h"
#include "AnisVisual.h"
#include "core/Interaction.h"
#include <GL/glu.h>
#include "Tensor.h"
#include <set>
#include <QtCore/QString>


using namespace GeoProperty;

class TriMesh;

class IncreaseOrder {
public: IncreaseOrder(int id, double val) : m_id(id), m_val(val) {}
		~IncreaseOrder() {}
		bool operator==(const IncreaseOrder& other) const { return (m_val == other.m_val) && (m_id == other.m_id); }
		bool operator>(const IncreaseOrder& other) const { return m_val == other.m_val ? m_id < other.m_id : m_val > other.m_val; }
		bool operator<(const IncreaseOrder& other) const { return m_val == other.m_val ? m_id < other.m_id : m_val < other.m_val; }
public:
	int		m_id;
	double	m_val;
};
class DecreaseOrder {
public: DecreaseOrder(int id, double val) : m_id(id), m_val(val) {}
		~DecreaseOrder() {}
		bool operator==(const DecreaseOrder& other) const { return (m_val == other.m_val) && (m_id == other.m_id); }
		bool operator>(const DecreaseOrder& other) const { return m_val == other.m_val ? m_id < other.m_id : m_val < other.m_val; }
		bool operator<(const DecreaseOrder& other) const { return m_val == other.m_val ? m_id < other.m_id : m_val > other.m_val; }
public:
	int		m_id;
	double	m_val;
};

typedef DecreaseOrder HighFunction;
typedef IncreaseOrder ClusterArea; //std::set.. use inc
class DecreaseOrderExt
{
public:
	DecreaseOrderExt(int id, double val, double val2) : m_id(id), m_val(val), m_val_ext(val2) {}
	DecreaseOrderExt(const DecreaseOrderExt& r) : m_id(r.m_id), m_val(r.m_val), m_val_ext(r.m_val_ext) {}
	~DecreaseOrderExt() {}
	bool operator=(const DecreaseOrderExt& other)
	{
		this->m_id = other.m_id;
		this->m_val = other.m_id;
		this->m_val_ext = other.m_val_ext;
	}
	bool operator<(const DecreaseOrderExt& other) const
	{
		return m_val_ext == other.m_val_ext ? m_id < other.m_id : m_val_ext > other.m_val_ext;
	}
	bool operator>(const DecreaseOrderExt& other) const
	{
		return m_val_ext == other.m_val_ext ? m_id < other.m_id : m_val_ext < other.m_val_ext;
	}
public:
	int		m_id;
	double	m_val;
	double	m_val_ext;
};

class MGraph {
public:
	struct Edge {
	public:
		Edge(int nid1 = 0, int nid2 = 0) :n1(nid1), n2(nid2) {}
		int n1, n2;
		double w, award;
		bool ort;
		bool isCut;
	};

	struct Vertex {
		int id;
		Vec3 coord;
		std::list<int> vadjs;
		std::list<int> eadjs;
	};

	MGraph() {}
	~MGraph() {}

public:

	int vNum;
	int eNum;
	std::vector<Edge> es;
	std::vector<Vertex> vs;
};

class Mitani_WaterShed
{
	// implement the paper "Making Papercraft Toys from Meshes
	//	using Strip - based Approximate Unfolding" by Jun Mitani  Hiromasa Suzuki. siggraph 2004.
	// and paper "Least Squares Conformal Maps for Automatic Texture Atlas Generation"
	// by Bruno L¨¦vy Sylvain Petitjean Nicolas Ray J¨¦rome Maillot.



public:
	Mitani_WaterShed() {};
	~Mitani_WaterShed() {};

private:


public:

	static void distanceToFeature_ShortestPath(std::list<MyMesh::Vertex>& vers, std::list<MyMesh::Edge>& edges, std::set<int>& sources, std::vector<double>& res);
	static void distanceToFeature_ShortestPath(MGraph& pg, std::set<int>& sources, std::vector<double>& res);

	// 	static void harmonicDistance(MyMesh *inMesh,std::vector<double>& weights,std::vector<int>& sources_sinks,std::vector<double>& res);

	static int grow(const std::vector<std::list<MyMesh::Vertex>::iterator>& vers, const std::vector<double>& distanceToFeature, std::vector<bool>& isExtrema, std::vector<unsigned>& vertexLabels);
	static int grow2(MGraph& pg, const std::vector<double>& distanceToFeature, const std::vector<std::set<unsigned> >& neighbours, std::vector<unsigned>& vertexLabels);

	static void thresholdClusters(MyMesh*, std::vector<unsigned>&, unsigned& segNum, std::vector<bool>&, double thres = 0);

public:
	std::vector<bool>		isFeature;
	std::vector<double>		distanceToFeature;
};

class MeshProperty {
public:

	MeshProperty() {
		myMesh = NULL; d_fitting = 4; d_monge = 4; crestline_scale = 2; //modify crestline scale
		sharpnessRangeLow = 1.0; sharpnessRangeHigh = 0.5; sharpnessRangeInterval = 1.0;
		metricVariation = 0; metricParameter = 0.05; //modify metric parameter
		showSDF = showGradSDF = false; showConcavity = false; gradLineWidth = 5; showSDFByText = showSDFLpp = false;
	}
	~MeshProperty() {}

	MyMesh** getMyMesh() { return &myMesh; }

protected:

	MyMesh * myMesh;

	unsigned	d_fitting;
	unsigned	d_monge;

	std::vector<double> sdf_faces, sdf_vts, sdf_grads;

	std::vector<double> con_vals;

	std::vector<Vec3> func_grad_vecs;
	std::vector<Tensor> func_tensors;
	std::vector<Vec3> vLps;
	std::vector<double> Lpps;

	std::vector<double> v_sharpness;
	std::vector<double> v_ridgeness;
	std::vector<double> v_anisotropy;

	std::vector<std::pair<double, double> > curMag;
	std::vector<std::vector<Vec3> > curDir;

	unsigned	metricVariation;
	double		metricParameter;
	std::vector<Tensor> faceAnis;

public:

	void computeVertexProperty_TRIMESH(double scale);
	void computeFaceTensor(std::vector<Tensor>&);

	void faceTensorUpdate(unsigned mt, double mp) { metricVariation = mt; metricParameter = mp; }
public:

	bool		showSDFByText, showSDFLpp;

	unsigned	crestline_scale;
	double		sharpnessRangeLow, sharpnessRangeHigh, sharpnessRangeInterval;
	bool		showSDF, showGradSDF, showConcavity;

	double gradLineWidth;
};
class FeatureLine : public MeshProperty
{

	struct CellFeature
	{
		Vec3 mid;
		Vec3 rep;
		bool isBoundary;
	};
	struct EdgeFeature
	{
		EdgeFeature() :isCut(false) {}
		Vec3 mid;
		Vec3 rep;
		int lab;

		bool	isOriented;
		double	dualLength;
		double	strength;

		bool	isCut;
	};
	///representation of features on a mesh.
	// it duals to a mesh edge, with each end duals to a triangle of mesh,
	// includes coordinates, thickness, length, the way it connected(orientation)...
	struct GraphFeature
	{
		//AwardFunction awardType;
		double alpha;

		std::vector<EdgeFeature> ef;
		std::vector<CellFeature> cf;
	};

public:
	FeatureLine()
	{
		showCrestLine = false;  graphFeature.alpha = 0.5;
		//lost function variations 
		//("1. Cost - a*Length"))
		//("2. Cost - a*Length*Anisotropic Strength"))
		//("3. Cost - a*Length*Ridgeness"))
		//("4. Cost - a*Length*Sharpness"))
		//("5. L- a*L\'"))
		//("6. L-(a+ L\')"))

		featureStrengthThres = 0;
		featureWidth = 2.0;
	}
	~FeatureLine() {}

protected:

	std::vector<Vec3> crestPoints;
	std::vector<unsigned> crestPointsAttriID;
	std::vector<std::vector<double> > crestAttris;
	std::vector<std::vector<unsigned> > crestEdges;
	std::vector<bool> crestEdgesVisible;
	unsigned crestLineoffset;

	std::vector<double> f_sharpness;
	std::vector<double> f_ridgeness;
	std::vector<double> f_anisotropy;

	//gc(graph cut), group features for graph cut
	std::vector<std::vector<unsigned> > gcFeatures;
	std::vector<std::vector<unsigned> > gcFeaturesOrig;
	std::vector<double > gcFeaturesStrength;


	//temp;
	std::vector<bool> isfeatureEdges;
	std::vector<std::vector<Vec3> > featureLines;

public:

	//void ComputeFeature(double curThres, double dirGradThres);
	//void GroupFeature();

	void computeMultiScaleFeatures(std::vector<bool>& boundaryEdge = std::vector<bool>(0));
	void computeCrestLine(std::vector<bool>& boundaryEdge = std::vector<bool>(0));
	void tuningCrestLine(double anisStrength, double sharpness, double ridgeness);

protected:

	void drawCrestLine();
	void featureGraphInitial();

public:

	bool			showCrestLine;
	double			featureStrengthThres;

	//features on mesh;
	GraphFeature	graphFeature;

	//for debug
	std::vector < bool > unusedCrestline;
	double featureWidth;

}; //feature lines
class MeshSegment : public FeatureLine
{
public:

	struct BoundaryVertex
	{
	public:
		BoundaryVertex() { ngbs.reserve(3);/*only valid in our algorithm; for comparison, this is not guaranteed;*/ }
		BoundaryVertex(const BoundaryVertex& v) { fid = v.fid; pos = v.pos; ngbs = v.ngbs; }

	public:

		int fid;
		Vec3 pos;
		std::vector<int> ngbs;
	};

	struct MGTriple {
		//int i;
		MGTriple(int a = -1, double b = 0.0) :j(a), val(b) {}
		int j;
		double val;
	};


	MeshSegment() { init(); }
	~MeshSegment() {}

	void init();

private:

	//initial
	std::vector<Vec3>					points;
	std::vector<std::vector<unsigned> > faces;
	//MyMesh*								myMesh;

	//compute curvatures;
	TriMesh * triMesh;

	//represented by a graph, rather than a graph.
	MGraph		mg;
	unsigned				featureType;

	std::vector<unsigned> vertexLabel;
	std::vector<unsigned> vertexLabelBeforeMerge;
	std::vector<unsigned> vertexLabelInit;
	unsigned segNumber;

	std::vector<double> weightGraphEdgeThickness;
	std::vector<ColorEngine::color> weightGraphColor;
	std::vector<ColorEngine::color> costColor;
	std::vector<ColorEngine::color> awardColor;

	std::vector<bool> subgraphVers;

	std::vector<MyInteraction::SketchFace> userSketches;

	//boundary smoothness
	std::vector<std::vector<Vec3> > boundaryCurves;
	std::vector<Vec3> boundaryJoints;


	//unsigned patchSmoothTimes;
	std::vector<double> v_difference;

	std::vector<ColorEngine::color> patchColors;

public:

	//in-out
	TriMesh * getTrimesh() { return triMesh; }
	//MyMesh** getMyMesh() { return &myMesh; }
	MyInteraction** getInteractionUtils() { return &m_interaction; }

	void rescaleMesh();
	void readTrimesh(const char* fileName);

	//compute weighted graph
	void computeCurvature();

	void updateGlabalAwardAlpha();
	void updateLocalAwardAlpha(std::vector<unsigned>& ToEnhance);
	void getWeights(std::vector<double> &, std::vector<double>&, std::vector<double>& t = std::vector<double>(0), std::vector<bool>& = std::vector<bool>(0));

	//clustering on weighted graph
	void overSegmentation(unsigned featureType = 0);
	//merge
	void mergePartition();
	//tiny patches remove
	void mergeSmallPatches();

	//for interaction. operation performs in sub-graph
	void modifySegmentBySketch(Interaction::SCRIBBLE_TYPE type);
	void smoothScribble(std::vector<unsigned>& sface, std::vector<Vec3>& scur);

	//post-process: smoothing on boundary.
	bool vertexMoveToLocalMinima(const std::vector<BoundaryVertex>& Vbs, unsigned i, BoundaryVertex&);
	void boundarySmooth(bool tmp = true);

	//opengl visualization
	void DrawGraph();
	void DrawCurve();

	//save result
	void saveSegmentation();

private:

	//opengl
	void genSegPatchColor(bool val);
	void drawGroupFeatures();
	void drawSegmentation();
	void drawSegBoundary();
	void drawUserSketch();
	void drawInteraction();

	/// mark graph edge with feature tag from input
	void featureFromCrestline();
	void featureFromSketches();

	//metric! weighting the edge
	void computeDualEdgeCostAndLength();
	void computeDualEdgeStrength();

	void featureToCurves();
	void mgraphInit(MGraph&);

	//clustering via wathershed under our anisotrpic metric on pri/dual graph;
	void Mitani_Watershed();
	void Mitani_Watershed_Dual();

	void convertFaceToVertexLabelling(std::vector<unsigned>&, std::vector<unsigned>& faceLabel);
	void convertVertexToFaceLabelling(std::vector<unsigned>&, std::vector<unsigned>& faceLabel);

	//merge on clusters
	void regionMerging(MGraph& gph, std::vector<unsigned>& vertexLabel, unsigned& labelId);
	void LMP_Merging(const std::vector<unsigned>& vertexLabel, std::vector< std::list< std::list< MGTriple > >::iterator >&pCIter,
		unsigned labelId, std::vector<unsigned>& newLabels);

	//partial segmentation; restrict the whole pipeline perform on a specific region/ sub-graph
	void Mitani_Watershed_Dual_Partial(std::vector<unsigned>&);
	void LMP_Partitioning();

	void getBoundaryOfClusters();

	//interactive operation and partial clustering...
	void findNeighboreOfTriangles(unsigned, std::vector<unsigned>&);
	void eraseFeatures(std::vector<unsigned>& fs);
	void eraseCluster(std::vector<unsigned>& es);
	void addSkethFeatures(std::vector<MyInteraction::SketchFace>& fs);

	void OverSegmentationOnPartialMesh(std::vector<unsigned>& fs, std::vector<unsigned>& es);
	void mergePartitionOnPartialMesh();
	void regionMergingOnPartialMesh(MGraph& gph, std::vector<bool>&, std::vector<unsigned>& vertexLabel, unsigned& labelId);
	void mergeLocalClusterByFeatureEnhancing(std::vector<unsigned>&);

	void genPatchColor(std::vector<ColorEngine::color>& patchColors, int segNumber, bool reGenerate);
	void patchColorMinChanged(std::vector<unsigned>& tlabels, int patNum);

public:

	QString m_filename;

	char								curvatureScale[100];
	char								curvatureSmooth[100];

	bool showSegmentation, showSegmentationBoundary, showWeightGraph;
	bool showMesh, showAnistropy, showSharpness, showRidgeness;

	double featureExtParam;

	unsigned gcFeaturePriority;
	bool gcIsMerge;

	std::vector<double> edgeWeightParameter;

	bool awardNormalized;

	bool isMergeSmallPatch;
	double smallPatchThres;

	bool autoBoundarySmooth;
	unsigned boundarySmoothTimes;

	//user interactive
	MyInteraction* m_interaction;

	std::set<unsigned> realTimeSelectedTriangles;
	std::vector<unsigned> sketchFaces;
	std::vector<Vec3> sketchCurves;
	double strokeStrength, paintParam;
	unsigned paintParamType, strokeSize;
	bool isSmoothStroke;

	//parameters for other methods;
	bool showWatershedField;
	Mitani_WaterShed ws;
	bool isWatershedMerge;
	double waterShedMergeParam;
	unsigned mitani_Watershed_Ringsize;
	bool isDualGraph_Watershed;
	bool isAnisGeodesics_Watershed;

	//for debug
	bool graphCutLocally;
	bool useAllFeatureLine;
	double boundaryWidth;
};

#endif