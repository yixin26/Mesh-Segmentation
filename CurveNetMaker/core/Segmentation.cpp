/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 //2014-11-02
 **************************************************************/
#include "time.h"
#include "core/Segmentation.h"
#include "setCurvature.h"

#include <Eigen/Dense>
#include <set>
#include "TriMesh.h"
#include "TriMesh_algo.h"

#include <unordered_map>
using namespace std;
using namespace GeoProperty;

#include <QtCore/QString>

//There are many libs for computing differential property of discrete surface, we includes CGAL and TriMesh;
//In default setting, we use Trimesh, which seems has smoother field over the surface;

void MeshSegment::computeCurvature()
{
	TriMesh *tmesh = new TriMesh(*triMesh);
	unsigned vNum = tmesh->vertices.size();
	tmesh->colors.resize(triMesh->vertices.size());
	colorbycurv(tmesh, curvatureScale, curvatureSmooth);

	unsigned i = 0; double maxAnis = 0;
	for (auto v_it = myMesh->getVertices().begin(); v_it != myMesh->getVertices().end(); v_it++)
	{
		v_it->normal() = Vec3(tmesh->normals[i][0], tmesh->normals[i][1], tmesh->normals[i][2]);
		v_it->color() = Vec3(tmesh->colors[i][0], tmesh->colors[i][1], tmesh->colors[i][2]);
		v_it->direction(0) = Vec3(tmesh->pdir1[i][0], tmesh->pdir1[i][1], tmesh->pdir1[i][2]);
		v_it->direction(1) = Vec3(tmesh->pdir2[i][0], tmesh->pdir2[i][1], tmesh->pdir2[i][2]);
		v_it->magnitude(0) = fabs(tmesh->curv1[i]);
		v_it->magnitude(1) = fabs(tmesh->curv2[i]);

		double& k1 = v_it->magnitude(0);
		double& k2 = v_it->magnitude(1);

		double diff = fabs(k1 - k2);

		if (maxAnis < diff)
		{
			maxAnis = diff;
		}
		i++;
	}
	AnisGeodesic::maxAnis = maxAnis;
	delete tmesh;
}
void MeshProperty::computeVertexProperty_TRIMESH(double scale)
{
	QString  mystring = QString::number(scale, 10, 3);

	//init trimesh;
	TriMesh *triMesh = new TriMesh;
	int vNum = myMesh->getVertices().size();
	auto vInds = myMesh->getVIter();
	triMesh->vertices.resize(vNum);
	for (int i = 0; i < vNum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			triMesh->vertices[i][j] = vInds[i]->coordinate()[j];
		}
	}
	int fNum = myMesh->getFaces().size();
	auto fInds = myMesh->getFIter();
	triMesh->faces.resize(fNum);
	for (int i = 0; i < fNum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			triMesh->faces[i][j] = fInds[i]->vertex_iter(j)->id();
		}
	}

	//compute curvature
	triMesh->need_normals();
	triMesh->colors.resize(triMesh->vertices.size());
	colorbycurv(triMesh, "0.0", mystring.toStdString().data());

	//write to...
	v_anisotropy.clear(); v_anisotropy.resize(myMesh->getVertices().size(), 0);
	curMag.clear(); curMag.resize(myMesh->getVertices().size());
	curDir.clear(); curDir.resize(myMesh->getVertices().size(), std::vector<Vec3>(3));

	unsigned i = 0;
	for (auto v_it = myMesh->getVertices().begin(); v_it != myMesh->getVertices().end(); v_it++, i++)
	{
		curDir[i][0] = Vec3(triMesh->normals[i][0], triMesh->normals[i][1], triMesh->normals[i][2]);
		curDir[i][1] = Vec3(triMesh->pdir1[i][0], triMesh->pdir1[i][1], triMesh->pdir1[i][2]);
		curDir[i][2] = Vec3(triMesh->pdir2[i][0], triMesh->pdir2[i][1], triMesh->pdir2[i][2]);
		curMag[i].first = fabs(triMesh->curv1[i]);
		curMag[i].second = fabs(triMesh->curv2[i]);
		v_anisotropy[i] = fabs(fabs(triMesh->curv1[i]) - fabs(triMesh->curv2[i]));
	}

	delete triMesh;
}
void MeshProperty::computeFaceTensor(std::vector<Tensor>& faceAnis)
{
	//avg of vertex tensor

	if (myMesh == NULL) return;

	cout << "compute tensor field on triangles(metric type:" << metricVariation << " metric parameter:" << metricParameter << ")" << endl;

	Tensor::computeTriangleCurvature(myMesh, faceAnis);

	if (metricVariation == 0)
	{
		for (unsigned i = 0; i < faceAnis.size(); i++)
		{
			Tensor& ft = faceAnis[i];
			double s = 1 + metricParameter*abs(ft.mag2 - ft.mag1);
			s *= s;
			ft.mag2 = 1 / s;
			ft.mag1 = 1;
		}
	}
	else if (metricVariation == 1)
	{
		for (unsigned i = 0; i < faceAnis.size(); i++)
		{
			Tensor& ft = faceAnis[i];
			double s = 1 + metricParameter*abs(ft.mag2 - ft.mag1);
			s *= s;
			ft.mag2 = 1;
			ft.mag1 = s;
		}
	}
	else
	{
		for (unsigned i = 0; i < faceAnis.size(); i++)
		{
			Tensor& ft = faceAnis[i];
			double s = 1 + metricParameter*abs(ft.mag2 - ft.mag1);
			s *= s;
			ft.mag2 = 1 / s;
			ft.mag1 = s;
		}
	}
}
void FeatureLine::computeMultiScaleFeatures(std::vector<bool>& boundaryEdge)
{
	/* compute multiscale local geometric features both for vertex tensors and feature prunning*/
	/*for vertex tensor, need to be avaraged by multi scale tensors*/

	unsigned numV = myMesh->getVertices().size();
	std::vector<double> sharpness(numV, 1);
	std::vector<double> ridgeness(numV, 1);
	std::vector<double> anisotropy(numV, 1);

	unsigned numIter = 0;
	for (double i = sharpnessRangeLow; i <= sharpnessRangeHigh; i += sharpnessRangeInterval, numIter++){}
	cout << numIter << " scales." << endl;

	numIter = 0;
	for (double i = sharpnessRangeLow; i <= sharpnessRangeHigh; i += sharpnessRangeInterval)
	{
		computeVertexProperty_TRIMESH(i);
		for (unsigned j = 0; j < numV; j++)
			anisotropy[j] += v_anisotropy[j];
	}
	numIter++;


	//avg of vers; f_xxx is for pruning of features;
	v_sharpness.swap(sharpness);
	v_ridgeness.swap(ridgeness);
	v_anisotropy.swap(anisotropy);
	f_sharpness.clear(); f_sharpness.resize(crestAttris[0].size(), 0);
	f_anisotropy.clear(); f_anisotropy.resize(crestAttris[0].size(), 0);
	f_ridgeness.clear(); f_ridgeness.resize(crestAttris[0].size(), 0);
	const auto& fts = myMesh->getFIter();
	for (unsigned i = 0; i < crestEdges.size(); i++)
	{
		if (crestEdges[i][2] > myMesh->getFaces().size()) continue;
		double esp = 0;
		double eani = 0;
		double erd = 0;
		for (unsigned j = 0; j < 3; j++)
		{
			esp += v_sharpness[fts[crestEdges[i][2]]->vertex_iter(j)->id()];
			eani += v_anisotropy[fts[crestEdges[i][2]]->vertex_iter(j)->id()];
			erd += v_ridgeness[fts[crestEdges[i][2]]->vertex_iter(j)->id()];
		}
		f_anisotropy[crestPointsAttriID[crestEdges[i][0]]] += eani;
		f_ridgeness[crestPointsAttriID[crestEdges[i][0]]] += erd;
		f_sharpness[crestPointsAttriID[crestEdges[i][0]]] += esp;
	}
}

//We use the lib by Shin Yoshizawa for computing Crestline(ridges&varines), theoretically detailed in Paper: Fast and Robust Detection of Crest Lines on Mesh;
void FeatureLine::computeCrestLine(std::vector<bool>& boundaryEdge)
{
	//we modify the interface to carry in our input, such as the neighbores of each vertex;
	//When user want to compute crestline inside a patch whose frontier are salient features globally,
	//she wants those salient features become less/none salient. This happens when we want to refine the segmentation iteratively, which requiring dynamic multi-sacle features.
	//This is simply done with boundaryEdge, that is constraints, that neighbores of a vertices can only lie on one side of which. 

	cout << "compute crest line(scale:" << crestline_scale << ")" << endl;

	if (boundaryEdge.empty()) boundaryEdge.resize(myMesh->getEdges().size(), false);

	unsigned vNum = myMesh->getVertices().size();
	unsigned fNum = myMesh->getFaces().size();

	//write vertices, faces, and adjcency of vers;
	double* vertices = new double[vNum * 3];
	unsigned i = 0;
	for (auto v_it = myMesh->getVertices().begin(); v_it != myMesh->getVertices().end(); v_it++)
	{
		vertices[i * 3 + 0] = v_it->coordinate().x;
		vertices[i * 3 + 1] = v_it->coordinate().y;
		vertices[i * 3 + 2] = v_it->coordinate().z;
		i++;
	}

	unsigned* faces = new unsigned[fNum * 3];
	i = 0;
	std::vector<unsigned> faceMap;
	for (auto if_it = myMesh->getFaces().begin(); if_it != myMesh->getFaces().end(); if_it++)
	{
		unsigned cnt = 0;
		for (unsigned j = 0; j<3; j++)
		{
			if (boundaryEdge[if_it->edge_iter(j)->id()])
			{
				cnt++;
			}
		}
		if (cnt>1) { fNum--; continue; }
		faceMap.push_back(if_it->id());

		faces[i * 3 + 0] = if_it->vertex_iter(0)->id();
		faces[i * 3 + 1] = if_it->vertex_iter(1)->id();
		faces[i * 3 + 2] = if_it->vertex_iter(2)->id();
		i++;
	}

	const auto& vts = myMesh->getVIter();
	std::vector<std::vector<int> > vadjs(myMesh->getVertices().size());
	for (auto v_it = myMesh->getVertices().begin(); v_it != myMesh->getVertices().end(); v_it++)
	{
		std::vector<int> usedVers;
		usedVers.push_back(v_it->id());
		std::vector<int> frontVerts = usedVers;
		unsigned count = 0;
		while (count<crestline_scale)
		{
			count++;
			std::vector<int> nextVerts(1, -1);
			while (!frontVerts.empty())
			{
				auto vft = vts[frontVerts.back()];
				frontVerts.pop_back();

				for (unsigned i = 0; i<vft->vertex_iter().size(); i++)
				{
					unsigned tid = vft->vertex_iter()[i]->id();
					if (std::find(usedVers.begin(), usedVers.end(), tid) == usedVers.end() && std::find(nextVerts.begin(), nextVerts.end(), tid) == nextVerts.end())
					{
						if (!boundaryEdge[vft->edge_iter()[i]->id()])
						{
							usedVers.push_back(tid); nextVerts.push_back(tid);
						}
					}
				}
			}
			nextVerts.erase(nextVerts.begin());
			frontVerts.swap(nextVerts);
		}
		usedVers.erase(usedVers.begin());
		vadjs[v_it->id()] = usedVers;
	}

	int** vadj = NULL;
	int* vadjNum = NULL;

	unsigned vn = vadjs.size();
	vadj = new int*[vn];
	vadjNum = new int[vn];
	for (int i = 0; i<vn; i++)
	{
		vadj[i] = new int[vadjs[i].size()];
		vadjNum[i] = vadjs[i].size();
		for (int j = 0; j<vadjs[i].size(); j++)		
		{
			vadj[i][j] = vadjs[i][j];
		}
	}

	unsigned pNum, cNum, eNum;
	double* fPoints;
	unsigned* fPointType;
	double* fFeature;
	unsigned* fEdges;
	unsigned offset;
	crestLineGen(vadj, vadjNum, vertices, vNum, faces, fNum, crestline_scale, 1,
		pNum, cNum, eNum, fPoints, fPointType, fFeature, fEdges, crestLineoffset);

	//read crestlines;
	Vec3 tp; unsigned tid;
	crestPoints.clear(); crestPointsAttriID.clear();
	for (unsigned int i = 0; i<pNum; ++i)
	{
		tp.x = fPoints[i * 3 + 0];
		tp.y = fPoints[i * 3 + 1];
		tp.z = fPoints[i * 3 + 2];
		crestPoints.push_back(tp);
		crestPointsAttriID.push_back(fPointType[i]);
	}
	crestAttris.clear(); crestAttris.resize(3);
	std::vector<double> &ridge = crestAttris[0];
	std::vector<double> &sphere = crestAttris[1];
	std::vector<double> &cynlinder = crestAttris[2];
	for (unsigned int i = 0; i<cNum; ++i)
	{
		ridge.push_back(fFeature[i * 3 + 0]);
		sphere.push_back(fFeature[i * 3 + 1]);
		cynlinder.push_back(fFeature[i * 3 + 2]);
	}
	std::vector<unsigned> te(3);
	crestEdges.clear();
	crestEdgesVisible.clear();
	for (unsigned i = 0; i<eNum; i++)
	{
		if (fEdges[i * 3 + 2] >= fNum) continue;

		te[0] = fEdges[i * 3 + 0];
		te[1] = fEdges[i * 3 + 1];
		te[2] = faceMap[fEdges[i * 3 + 2]];
		crestEdges.push_back(te);
		crestEdgesVisible.push_back(true);
	}

	delete vertices;
	delete faces;
	for (int i = 0; i<vn; i++)
		delete vadj[i];
	delete vadjNum;
	delete fPointType;
	delete fFeature;
	delete fEdges;
}
void FeatureLine::tuningCrestLine(double anisStrength, double sharpness, double ridgeness)
{
	//datas are kept, filtered is active by crestEdgesVisible. 
	//invisible crestEdge will not be shown and used for segmentation.

	std::vector<double> &anis = f_anisotropy;
	std::vector<double> &ridge = f_ridgeness;
	std::vector<double> &sharp = f_sharpness;

	std::vector<bool> visibleAttri(ridge.size(), true);

	double maxVal, minVal;
	maxVal = *std::max_element(sharp.begin(), sharp.end());
	minVal = *std::min_element(sharp.begin(), sharp.end());
	double rang = maxVal - minVal;
	for (unsigned i = 0; i<ridge.size(); i++)
	{
		if (rang>0 && ((sharp[i] - minVal) / rang)<sharpness)
		{
			visibleAttri[i] = false;
		}
	}

	maxVal = *std::max_element(anis.begin(), anis.end());
	minVal = *std::min_element(anis.begin(), anis.end());
	rang = maxVal - minVal;
	for (unsigned i = 0; i<anis.size(); i++)
	{
		if (rang>0 && ((anis[i] - minVal) / rang)<anisStrength)
		{
			visibleAttri[i] = false;
		}
	}

	maxVal = *std::max_element(ridge.begin(), ridge.end());
	minVal = *std::min_element(ridge.begin(), ridge.end());
	rang = maxVal - minVal;
	for (unsigned i = 0; i<ridge.size(); i++)
	{
		if (rang>0 && ((ridge[i] - minVal) / rang)<ridgeness)
		{
			visibleAttri[i] = false;
		}
	}

	for (unsigned i = 0; i<crestEdgesVisible.size(); i++)
	{
		unsigned &v1 = crestEdges[i][0];
		unsigned &v2 = crestEdges[i][1];
		if (visibleAttri[crestPointsAttriID[v1]] == false || visibleAttri[crestPointsAttriID[v2]] == false)
			crestEdgesVisible[i] = false;
		else
			crestEdgesVisible[i] = true;
	}

	cout << "crestline updated" << endl;
}

void FeatureLine::featureGraphInitial()
{
	//initialization of graph, representing nodes,edges,cells and the features.
	if (myMesh == NULL) return;

	std::vector<EdgeFeature>& tef = graphFeature.ef;
	std::vector<CellFeature>& tcf = graphFeature.cf;
	tef.clear(); tcf.clear();
	tef.resize(myMesh->getEdges().size());
	tcf.resize(myMesh->getFaces().size());
	for (auto f_it = myMesh->getFaces().begin(); f_it != myMesh->getFaces().end(); f_it++)
	{
		tcf[f_it->id()].mid = f_it->center_point();
		tcf[f_it->id()].rep = f_it->center_point();
		tcf[f_it->id()].isBoundary = false;
	}
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		tef[e_it->id()].mid = (e_it->vertex_iter(0)->coordinate() + e_it->vertex_iter(1)->coordinate())*0.5;
		tef[e_it->id()].rep = tef[e_it->id()].mid;
		tef[e_it->id()].lab = -1;
		tef[e_it->id()].isOriented = true;
		tef[e_it->id()].dualLength = -1;
		tef[e_it->id()].strength = 1;
// 		tef[e_it->id()].isCut = false;
	}
}

void MeshSegment::init()
{
	showMesh=true;
	showAnistropy=showRidgeness=showSharpness=false;
	crestLineoffset=0;

	crestPoints.clear(); crestPointsAttriID.clear();crestAttris.clear();
	crestEdges.clear();crestEdgesVisible.clear();

	showWeightGraph=false; showSegmentation=showSegmentationBoundary=true;
	featureExtParam = 0.0;

	featureType = 1;
	awardNormalized = false;

	gcFeaturePriority=1;
	gcIsMerge=true;

	segNumber = 1;
	vertexLabel.clear();

	isWatershedMerge = true;
	isMergeSmallPatch = true;
	smallPatchThres= 5.0;
	waterShedMergeParam = 3;

	autoBoundarySmooth=true;
	boundarySmoothTimes=3;
	v_difference.clear();

	boundaryCurves.clear();
	boundaryJoints.clear();

	strokeStrength=0.5;
	paintParam = 0.05;
	paintParamType = 0;

	showWatershedField = false;
	mitani_Watershed_Ringsize  = 1;
	isAnisGeodesics_Watershed = true;
	isDualGraph_Watershed = true;

	useAllFeatureLine = false;

	boundaryWidth = 2.0;

	m_interaction = NULL;
	strokeSize = 1;
	isSmoothStroke = true;

	graphCutLocally = false;

	points.clear();	faces.clear(); myMesh = NULL;

	triMesh = NULL; strcpy(curvatureScale, "0.");	strcpy(curvatureSmooth, "0.1");

}

void MeshSegment::featureFromCrestline()
{
	//crestline are not clean! may includes invalid/incorret indices;
	//only visible and valid feature edges are wrote to feature graph.
	if(myMesh == NULL) return;

	unsigned cnt = 0;	unusedCrestline.clear(); unusedCrestline.resize(crestEdgesVisible.size(), false);//crestlines are shown in blue and red, invalid ones are shown in yellow. this is for debug.
	for(unsigned i=0;i<crestEdgesVisible.size();i++)
	{
		unusedCrestline[i] = true;	cnt++;

		if (!crestEdgesVisible[i])
			continue;

		unsigned faceID = crestEdges[i][2];
		if(faceID>myMesh->getFaces().size())
		{
			crestEdgesVisible[i] = false;
			continue; // crest line data is not alway clean
		}

		unusedCrestline[i] = false;	cnt--;

		auto f_it = myMesh->getFIter()[faceID];

		std::vector<unsigned> edgeUsed;
		for (unsigned j=0;j<2;j++)
		{
			unsigned vid = crestEdges[i][j];
			Vec3 voe = crestPoints[crestEdges[i][j]];
			//simply check if the crestPoint is on the edge;
			for (unsigned k=0;k<3;k++)
			{
				double dis=0;
				Vec3 v1 = f_it->edge_iter(k)->vertex_iter(0)->coordinate();
				Vec3 v2 = f_it->edge_iter(k)->vertex_iter(1)->coordinate();
				Vec3 v12 = v1-v2;
				v1 = voe-v1; v2=voe-v2;
				if(abs(v1.length()+v2.length()-v12.length())< 1e-7) //when very close, then it is
				{
					graphFeature.ef[f_it->edge_iter(k)->id()].lab = 1;
					edgeUsed.push_back(f_it->edge_iter(k)->id());
				}
			}
		}
		if(edgeUsed.size()==4)
		{//feature edge is on the edge;
			unsigned tid;
			for (unsigned j=0;j<3;j++)
			{
				for(unsigned k=j+1;k<4;k++)
				{
					if(edgeUsed[j]==edgeUsed[k])
					{
						graphFeature.ef[edgeUsed[k]].lab = -1;
						k=4;j=4;
					}
				}
			}
		}
		else if (edgeUsed.size() == 3)
		{//feature edge go through one vertex;
			unsigned vind;
			double dis = DBL_MAX;
			unsigned tind;
			for (unsigned k = 0; k < 3; k++)
			{
				Vec3 v1 = f_it->vertex_iter(0)->coordinate();
				double tdis = min(v1.dot(crestPoints[crestEdges[i][0]]), v1.dot(crestPoints[crestEdges[i][1]]));
				if (tdis < dis)
				{
					dis = tdis; vind = k; 
					tind = v1.dot(crestPoints[crestEdges[i][0]]) < v1.dot(crestPoints[crestEdges[i][1]]) ? 1 : 0;
				}
			}
			unsigned te = f_it->opposite_edge(f_it->vertex_iter(vind))->id();
			graphFeature.cf[faceID].rep = (crestPoints[crestEdges[i][0]] + crestPoints[crestEdges[i][1]])*0.5;
			graphFeature.ef[te].rep = crestPoints[crestEdges[i][tind]];

			te = f_it->opposite_edge(f_it->vertex_iter((vind+1)%3))->id();
			graphFeature.ef[te].rep = crestPoints[crestEdges[i][(tind+1)%2]];
		}
		else if (edgeUsed.size() == 2)
		{
			//general case;
			graphFeature.cf[faceID].rep = (crestPoints[crestEdges[i][0]] + crestPoints[crestEdges[i][1]])*0.5;
			graphFeature.ef[edgeUsed.front()].rep = crestPoints[crestEdges[i][0]];
			graphFeature.ef[edgeUsed.back()].rep = crestPoints[crestEdges[i][1]];
		}
		else
		{
			cnt++;	unusedCrestline[i] = true;
		}
	}

	if (cnt != 0)	cout << "bad crestline edges:" << cnt << endl;
}
void MeshSegment::featureFromSketches()
{
	//user input sketches as new features;
	for(unsigned i=0;i<userSketches.size();i++)
	{
		unsigned faceID = userSketches[i].fid;
		if (faceID>myMesh->getFaces().size()) continue; //double checks

		graphFeature.ef[userSketches[i].e1].lab = 2;
		graphFeature.ef[userSketches[i].e2].lab = 2;

		graphFeature.cf[faceID].rep=userSketches[i].fpos;
		graphFeature.ef[userSketches[i].e1].rep = userSketches[i].ep1;
		graphFeature.ef[userSketches[i].e2].rep = userSketches[i].ep2;
	}
}

void MeshSegment::computeDualEdgeStrength()
{
	double maxCost = 0;
	//the dual length is not the euclidean length of vector, but the length under our metric.
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		graphFeature.ef[e_it->id()].strength = 1;
		if (graphFeature.ef[e_it->id()].lab == -1)//a normal edge
		{
			graphFeature.ef[e_it->id()].dualLength = 0;
		}
		else//a feature edge or sketched edge
		{
			if (graphFeature.ef[e_it->id()].lab == 1)
			{
				double edgeCost = 0;
				if (e_it->manifold())
				{
					MyMesh::FaceIter f[] = { e_it->face_iter(0), e_it->face_iter(1) };
					unsigned fid[] = { f[0]->id(), f[1]->id() };
					Vec3 edgeVec = graphFeature.cf[f[0]->id()].rep - graphFeature.cf[f[1]->id()].rep;
					for (unsigned k = 0; k < 2; k++)
					{
						Vec2 edgeVec2D = Vec2(edgeVec.dot(faceAnis[fid[k]].dir1), edgeVec.dot(faceAnis[fid[k]].dir2));
						edgeCost += sqrt(pow(edgeVec2D.x, 2)*faceAnis[fid[k]].mag2 + pow(edgeVec2D.y, 2)*faceAnis[fid[k]].mag1);// two directions of tensor are swapped, since we computed dual edge's cost; mag1>mag2;
					}
					edgeCost = edgeCost*0.5;
				}
				else
				{
					MyMesh::FaceIter f[] = { e_it->face_iter(0) };
					unsigned fid[] = { f[0]->id(), f[1]->id() };
					Vec3 edgeVec = graphFeature.cf[f[0]->id()].rep - graphFeature.ef[e_it->id()].rep;
					for (unsigned k = 0; k < 1; k++)
					{
						Vec2 edgeVec2D = Vec2(edgeVec.dot(faceAnis[fid[k]].dir1), edgeVec.dot(faceAnis[fid[k]].dir2));// two direction were swapped, since we computed dual edge's cost; and in faceAnis, mag1>mag2;
						edgeCost += sqrt(pow(edgeVec2D.x, 2)*faceAnis[fid[k]].mag2 + pow(edgeVec2D.y, 2)*faceAnis[fid[k]].mag1);
					}
				}
				graphFeature.ef[e_it->id()].dualLength = edgeCost;
				maxCost = max(maxCost, edgeCost);
			}
		}
		maxCost = max(maxCost, e_it->cost());
	}
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		if (graphFeature.ef[e_it->id()].lab == 2)
		{
			graphFeature.ef[e_it->id()].dualLength = 10 * maxCost;
		}
	}
}
void MeshSegment::computeDualEdgeCostAndLength()
{
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		if (e_it->manifold())
		{
			MyMesh::FaceIter f[] = { e_it->face_iter(0), e_it->face_iter(1) };
			unsigned fid[] = { f[0]->id(), f[1]->id() };
			Vec3 edgeVec = graphFeature.cf[f[0]->id()].rep - graphFeature.cf[f[1]->id()].rep;
			double edgeCost = 0;
			for (unsigned k = 0; k < 2; k++)
			{
				Vec2 edgeVec2D = Vec2(edgeVec.dot(faceAnis[fid[k]].dir1), edgeVec.dot(faceAnis[fid[k]].dir2));
				edgeCost += sqrt(pow(edgeVec2D.x, 2)*faceAnis[fid[k]].mag1 + pow(edgeVec2D.y, 2)*faceAnis[fid[k]].mag2);
			}
			e_it->cost() = edgeCost*0.5;
			//feature award;
			if (graphFeature.ef[e_it->id()].lab != -1)
			{
				graphFeature.ef[e_it->id()].dualLength = edgeVec.length();
			}
		}
		else
		{
			MyMesh::FaceIter f[] = { e_it->face_iter(0) };
			unsigned fid[] = { f[0]->id(), f[1]->id() };
			Vec3 edgeVec = graphFeature.cf[f[0]->id()].rep - graphFeature.ef[e_it->id()].rep;
			double edgeCost = 0;
			for (unsigned k = 0; k < 1; k++)
			{
				Vec2 edgeVec2D = Vec2(edgeVec.dot(faceAnis[fid[k]].dir1), edgeVec.dot(faceAnis[fid[k]].dir2));// two direction were swapped, since we computed dual edge's cost; and in faceAnis, mag1>mag2;
				edgeCost += sqrt(pow(edgeVec2D.x, 2)*faceAnis[fid[k]].mag1 + pow(edgeVec2D.y, 2)*faceAnis[fid[k]].mag2);
			}
			e_it->cost() = edgeCost + 1e-6; //make sure weight is bigger than 0;
			//feature award;
			if (graphFeature.ef[e_it->id()].lab != -1)
			{
				graphFeature.ef[e_it->id()].dualLength = edgeVec.length();
			}
		}
	}
}
void MeshSegment::updateGlabalAwardAlpha()
{
	if (!edgeWeightParameter.empty())
	{
		double brushParameter = paintParam;
		if (paintParamType == 0)
		{
			if (m_interaction->scribbleType == Interaction::CONCEAL)
				brushParameter = -paintParam;
			for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
			{
				if (graphFeature.ef[e_it->id()].lab != -1)//a normal edge
				{
					edgeWeightParameter[e_it->id()] += brushParameter;
				}
			}
			graphFeature.alpha += brushParameter;
		}
		else if (paintParamType == 1)
		{
			if (m_interaction->scribbleType == Interaction::CONCEAL)
				brushParameter = 1.0 / paintParam;
			for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
			{
				if (graphFeature.ef[e_it->id()].lab != -1)//a normal edge
				{
					edgeWeightParameter[e_it->id()] *= brushParameter;
				}
			}
			graphFeature.alpha *= brushParameter;
		}
	}
}
void MeshSegment::updateLocalAwardAlpha(std::vector<unsigned>& fs)
{
	cout << "Boost/Decrease Alpha by ";
	if (paintParamType == 0)
	{
		if (m_interaction->scribbleType == Interaction::REVEAL)
			cout << "+";
		else
			cout << "-";
	}
	else if (paintParamType == 1)
	{
		if (m_interaction->scribbleType == Interaction::REVEAL)
			cout << "*";
		else
			cout << "/";
	}
	if (paintParamType != 2)
		cout << " " << paintParam << endl;
	else
		cout << "a fix number:" << 1e-20 << endl;

	const auto& faceIters = myMesh->getFIter();
	std::set<int> enhancedEdge;
	for (int i = 0; i < fs.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			auto eit = faceIters[fs[i]]->edge_iter(j);
			if (graphFeature.ef[eit->id()].lab != -1)
			{
				enhancedEdge.insert(eit->id());
			}
		}
	}
	double brushParameter = paintParam;
	if (paintParamType == 0)
	{
		if (m_interaction->scribbleType == Interaction::CONCEAL)
			brushParameter = -paintParam;
		for (auto i = enhancedEdge.begin(); i != enhancedEdge.end(); i++)
		{
			edgeWeightParameter[*i] += brushParameter;
		}
	}
	else if (paintParamType == 1)
	{
		if (m_interaction->scribbleType == Interaction::CONCEAL)
			brushParameter = 1.0 / paintParam;
		for (auto i = enhancedEdge.begin(); i != enhancedEdge.end(); i++)
		{
			edgeWeightParameter[*i] *= brushParameter;
		}
	}
	else
	{
		for (auto i = enhancedEdge.begin(); i != enhancedEdge.end(); i++)
		{
			edgeWeightParameter[*i] = 1e-20;
		}
	}
	cout << "max weight param:" << *std::max_element(edgeWeightParameter.begin(), edgeWeightParameter.end()) << endl;
	cout << "min weight param:" << *std::min_element(edgeWeightParameter.begin(), edgeWeightParameter.end()) << endl;
}
void MeshSegment::getWeights(std::vector<double>& param, std::vector<double>& edgeWeight, std::vector<double>& edgeAward, std::vector<bool>& featureEdges)
{
	if (param.empty())
	{
		edgeWeightParameter.resize(myMesh->getEdges().size(), graphFeature.alpha);
	}

	edgeWeight.clear(); edgeWeight.resize(myMesh->getEdges().size());
	edgeAward.clear(); edgeAward.resize(myMesh->getEdges().size(), 0);
	featureEdges.clear(); featureEdges.resize(myMesh->getEdges().size(), false);
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		auto tc = e_it->cost();

		double edgeParam = param[e_it->id()];
		if (graphFeature.ef[e_it->id()].lab != -1)
		{
			if (awardNormalized)
			{
				double normalizationParam = edgeParam;
				double anisVal = graphFeature.ef[e_it->id()].strength;
				double K = 1 - exp(-0.5 * normalizationParam*(anisVal));
				edgeAward[e_it->id()] = tc * (1 + K);
			}
			else
			{
				edgeAward[e_it->id()] = edgeParam * graphFeature.ef[e_it->id()].dualLength *graphFeature.ef[e_it->id()].strength;
			}

			edgeWeight[e_it->id()] = tc - edgeAward[e_it->id()];
			featureEdges[e_it->id()] = true;
		}
		else
		{
			edgeWeight[e_it->id()] = tc;
		}
	}
}

void MeshSegment::featureToCurves()
{
	//group features to curves;

	std::vector<std::vector<unsigned> > featureEdges;
	std::vector<std::vector<unsigned> > edgeInFaces(myMesh->getFaces().size());
	//this is only correct in an oriented mesh;
	for(auto e_it=myMesh->getEdges().begin();e_it!=myMesh->getEdges().end();e_it++)
	{
		if(graphFeature.ef[e_it->id()].lab ==-1 || !e_it->manifold()) continue;

		std::vector<unsigned> newEdge;
		if(e_it->face_iter(0)->next_vertex(e_it->vertex_iter(0))==e_it->vertex_iter(1))
		{
			newEdge.push_back(e_it->face_iter(0)->id());
			newEdge.push_back(e_it->face_iter(1)->id());
		}
		else
		{
			newEdge.push_back(e_it->face_iter(1)->id());
			newEdge.push_back(e_it->face_iter(0)->id());
		}
		newEdge.push_back(e_it->id());
		featureEdges.push_back(newEdge);

		edgeInFaces[newEdge[0]].push_back(featureEdges.size()-1);
		edgeInFaces[newEdge[1]].push_back(featureEdges.size()-1);
	}

	//feature lines;
	featureLines.clear();

	gcFeatures.clear();
	std::vector<bool> edgeUsed(featureEdges.size(),false);
	unsigned groupNum=0;
	while(true)
	{

		unsigned eid = std::find(edgeUsed.begin(),edgeUsed.end(),false) - edgeUsed.begin();
		if(eid>=edgeUsed.size()) break;

		edgeUsed[eid]=true;
		//edgeFeatureLabels[featureEdges[eid][2]] = groupNum;

		std::vector<unsigned> gcFeature;
		gcFeature.push_back(featureEdges[eid][2]);
		std::vector<Vec3> featureLine;
		featureLine.push_back(graphFeature.ef[featureEdges[eid][2]].rep);
		//find group feature
		unsigned frontEdge = eid;
		unsigned frontFace = featureEdges[eid][1];
		while(true)
		{
			if(edgeInFaces[frontFace].size()==1) break;
			unsigned nextEdge;
			for(unsigned i=0;i<edgeInFaces[frontFace].size();i++)
			{
				if(edgeInFaces[frontFace][i]!=frontEdge)
				{
					nextEdge = edgeInFaces[frontFace][i]; 
					break;
				}
			}
			if(edgeUsed[nextEdge]) break;

			edgeUsed[nextEdge] = true;
			//edgeFeatureLabels[featureEdges[nextEdge][2]] = groupNum;
			gcFeature.push_back(featureEdges[nextEdge][2]);
			featureLine.push_back(graphFeature.ef[featureEdges[nextEdge][2]].rep);

			if(featureEdges[nextEdge][0]==frontFace)
			{
				frontFace = featureEdges[nextEdge][1];
			}
			else
			{
				frontFace = featureEdges[nextEdge][0];
				graphFeature.ef[featureEdges[nextEdge][2]].isOriented = false;
			}
			frontEdge=nextEdge;
		}
		reverse(gcFeature.begin(), gcFeature.end());
		reverse(featureLine.begin(), featureLine.end());

		frontEdge = eid;
		frontFace = featureEdges[eid][0];
		while(true)
		{
			if(edgeInFaces[frontFace].size()==1) break;
			unsigned nextEdge;
			for(unsigned i=0;i<edgeInFaces[frontFace].size();i++)
			{
				if(edgeInFaces[frontFace][i]!=frontEdge)
				{
					nextEdge = edgeInFaces[frontFace][i]; break;
				}
			}
			if(edgeUsed[nextEdge]) break;

			edgeUsed[nextEdge] = true;
			//edgeFeatureLabels[featureEdges[nextEdge][2]] = groupNum;
			gcFeature.push_back(featureEdges[nextEdge][2]);
			featureLine.push_back(graphFeature.ef[featureEdges[nextEdge][2]].rep);

			if(featureEdges[nextEdge][1]==frontFace)
			{
				frontFace = featureEdges[nextEdge][0];
			}
			else
			{
				frontFace = featureEdges[nextEdge][1];
				graphFeature.ef[featureEdges[nextEdge][2]].isOriented = false;
			}
			frontEdge=nextEdge;
		}
		if (!useAllFeatureLine && gcFeature.size() < 3)
		{
			for (int f = 0; f < gcFeature.size(); f++)
			{
				graphFeature.ef[gcFeature[f]].lab = -1;
			}

		}
		else
		{
			groupNum++;
			gcFeatures.push_back(gcFeature);

			featureLines.push_back(featureLine);
		}

	}
}

void MeshSegment::mgraphInit(MGraph& mg)
{
	//a data sturcture for graph cut. 

	std::vector<double> edgeWeight(myMesh->getEdges().size());
	std::vector<double> edgeAward(myMesh->getEdges().size(), 0);
	std::vector<bool> featureEdges(myMesh->getEdges().size(), false);
	getWeights(edgeWeightParameter,edgeWeight, edgeAward, featureEdges);

	mg.vNum = myMesh->getVertices().size();
	mg.vs.clear(); mg.vs.resize(mg.vNum);
	mg.eNum = myMesh->getEdges().size();
	mg.es.clear(); mg.es.resize(mg.eNum);

	int i=0;
	for(auto vit = myMesh->getVertices().begin(); vit!=myMesh->getVertices().end(); vit++)
	{
		mg.vs[i].id = i;
		for(int j=0;j<vit->vertex_iter().size();j++)
		{
			mg.vs[i].vadjs.push_back(vit->vertex_iter()[j]->id());
			mg.vs[i].eadjs.push_back(vit->edge_iter()[j]->id());
		}
		i++;
	}
	i=0;
	try
	{
		for (auto eit = myMesh->getEdges().begin(); eit != myMesh->getEdges().end(); eit++)
		{
			mg.es[i].isCut = featureEdges[i];
			mg.es[i].n1 = eit->vertex_iter(0)->id();
			mg.es[i].n2 = eit->vertex_iter(1)->id();
			if (std::isnan(edgeWeight[i]) || std::isnan(edgeAward[i]))
			{
				throw ;
			}
			mg.es[i].w = edgeWeight[i];
			mg.es[i].award = edgeAward[i];
			mg.es[i].ort = graphFeature.ef[i].isOriented;
			i++;
		}
	}
	catch (double e)
	{
		std::cout << "The input graph weight has NAN...\n";
	}
}

//////////////////////////////////////////////////////////////////////////
//compute over segmentation of mesh
//requires everything of features, initialization of graph, partitioning algorithms, merging algorithms, and final smoothing.
void MeshSegment::overSegmentation(unsigned ftype)
{
	cout << "-------------------- GRAPH PARTITION -------------" << endl;
	//1. init GraphFeature data structure;
	{
		clock_t tstr = clock();
		featureGraphInitial();
		featureFromCrestline();
		featureFromSketches();

		computeFaceTensor(faceAnis);
		computeDualEdgeCostAndLength();
		computeDualEdgeStrength();

		featureToCurves();

		mgraphInit(mg);

		clock_t tinit = clock() - tstr;
		cout << "alpha:" << graphFeature.alpha << endl;
		cout << endl << "graph initialization:" << tinit / 1000 << "sec" << tinit % 1000 << "mm" << endl;
	}


	if (isDualGraph_Watershed)
		Mitani_Watershed_Dual();
	else
		Mitani_Watershed();

	genPatchColor(patchColors, segNumber,false);
	getBoundaryOfClusters();

	vertexLabelInit = vertexLabel;
}
void MeshSegment::regionMerging(MGraph& gph,
	std::vector<unsigned>& vertexLabel,
	unsigned& labelId)
{
	clock_t tstr = clock();

	labelId = *std::max_element(vertexLabel.begin(), vertexLabel.end()) + 1;

	std::vector<std::set<int> > adjVers(labelId);
	for (int i = 0; i < (int)gph.es.size(); i++)
	{
		unsigned vid0 = min(vertexLabel[gph.es[i].n1], vertexLabel[gph.es[i].n2]);
		unsigned vid1 = max(vertexLabel[gph.es[i].n1], vertexLabel[gph.es[i].n2]);
		if (vid0 == vid1)continue;
		adjVers[vid0].insert(vid1);
	}


	std::list< std::list< MGTriple > > pairClusterCost(labelId);
	typedef std::list< std::list< MGTriple > >::iterator PCIter;
	std::vector< PCIter > pCIter;

	int count = 0;
	MGTriple tmg;
	for (PCIter pct = pairClusterCost.begin(); pct != pairClusterCost.end(); pct++, count++)
	{
		if (adjVers[count].empty())
		{
			tmg.j = count;
			tmg.val = 0;
			pct->push_back(tmg);
		}
		else
		{
			for (std::set<int>::iterator ita = adjVers[count].begin(); ita != adjVers[count].end(); ita++)
			{
				tmg.j = *ita;
				tmg.val = 0;
				pct->push_back(tmg);
			}
		}
		pCIter.push_back(pct);
	}

	std::list< MGTriple >::iterator tpc;

	for (int i = 0; i < (int)gph.es.size(); i++)
	{
		unsigned vid0 = min(vertexLabel[gph.es[i].n1], vertexLabel[gph.es[i].n2]);
		unsigned vid1 = max(vertexLabel[gph.es[i].n1], vertexLabel[gph.es[i].n2]);
		if (vid0 == vid1)continue;

		for (tpc = pCIter[vid0]->begin(); tpc != pCIter[vid0]->end(); tpc++)
		{
			if (tpc->j == vid1) break;
		}

		tpc->val += gph.es[i].w;
	}

	//algorithms;
	std::vector<unsigned> newLabels;
	LMP_Merging(vertexLabel, pCIter, labelId, newLabels);

	vertexLabel.clear(); vertexLabel.resize(gph.vNum);
	labelId = 0;
	std::vector<bool> visitedVer(vertexLabel.size(), false);
	for (size_t v = 0; v < gph.vNum; v++)
	{
		if (visitedVer[v] == true) continue;

		unsigned currentLabel = unsigned(newLabels[v]);
		std::vector<MGraph::Vertex*> frontVer(1, &(gph.vs[v]));
		while (!frontVer.empty())
		{
			MGraph::Vertex* fv = frontVer.back(); frontVer.pop_back();
			vertexLabel[fv->id] = labelId;
			for (std::list<int>::iterator i = fv->vadjs.begin(); i != fv->vadjs.end(); i++)
			{
				MGraph::Vertex* adjv = &(gph.vs[*i]);
				if (visitedVer[adjv->id] == true) continue;
				if (unsigned(newLabels[adjv->id]) == currentLabel)
				{
					visitedVer[adjv->id] = true;
					frontVer.push_back(adjv);
				}
			}
		}
		labelId++;
	}

	cout << endl << "cluster number after merging:" << labelId << endl;

	//get boundary;
	double cutCost = 0;
	for (int i = 0; i < (int)gph.es.size(); i++)
	{
		unsigned vs[] = { gph.es[i].n1, gph.es[i].n2 };
		if (vertexLabel[vs[0]] != vertexLabel[vs[1]])
		{
			gph.es[i].isCut = true;
			cutCost += gph.es[i].w;
		}
		else
		{
			gph.es[i].isCut = false;
		}
	}

	clock_t tinit = clock() - tstr;
	cout << "merging takes:" << tinit / 1000 << "sec" << tinit % 1000 << "mm" << endl;
	cout << "Cost:" << cutCost << endl;
}
void MeshSegment::mergePartition()
{
	boundaryCurves.clear();
	if(!gcIsMerge) return;
	cout << "*   *   *   *   *   *  MERGE  *   *   *   *   *   *" << endl;

	//just for color 
	std::vector<unsigned> tlabels = vertexLabel;
	unsigned patNum = segNumber;

	//merging
	vertexLabel=vertexLabelBeforeMerge;

	cout << "alpha:" << graphFeature.alpha << endl;
	mgraphInit(mg);

	double cutCost;
	regionMerging(mg, vertexLabel, segNumber);

	if (isMergeSmallPatch)
	{
		mergeSmallPatches();
	}
	else
	{
		for(int i=0;i<mg.eNum;i++)
		{
			graphFeature.ef[i].isCut = mg.es[i].isCut;
		}
	}

	if (autoBoundarySmooth)
	{
		cout << endl;
		boundarySmooth(true);
	}

	cout << "\n-----------------------------------------------------" << endl;

	//label id for most of vertices should remain the same, causing least changes;
	patchColorMinChanged(tlabels, patNum);
}

void MeshSegment::mergeSmallPatches()
{
	//prioritized by patch size, similar as greedy merging, by weight between patches.

	auto& mesh = myMesh;
	std::vector<unsigned>& vertexLabels = vertexLabel;
	double thres = smallPatchThres;
	unsigned& segNum = segNumber;

	double averageArea = 0;
	std::vector<double> faceAreas(mesh->getFaces().size());
	for (auto i = mesh->getFaces().begin(); i != mesh->getFaces().end(); i++)
	{
		faceAreas[i->id()] = i->triangleCost(2);
		averageArea += faceAreas[i->id()];
	}

	{//in this case, we threshold by triangle number, this is for our algorithm.
		averageArea /= double(faceAreas.size());
		thres *= averageArea;
	}

	//initialize area of all clusters;
	std::vector<double> scores(segNum, 0.0);
	for (auto i = mesh->getFaces().begin(); i != mesh->getFaces().end(); i++)
	{
		unsigned v1 = i->vertex_iter(0)->id();
		unsigned v2 = i->vertex_iter(1)->id();
		unsigned v3 = i->vertex_iter(2)->id();
		if (vertexLabels[v1] == vertexLabels[v2] && vertexLabels[v1] == vertexLabels[v3])
		{
			scores[vertexLabels[v1]] += faceAreas[i->id()];
		}
	}

	//shared boundary length between clusters;
	std::vector<double> sharesRow(segNum, 0.0);
	std::vector<std::vector<double> > sharesMatrix(segNum, sharesRow);
	for (auto i = mesh->getEdges().begin(); i != mesh->getEdges().end(); i++)
	{
		int vid0 = vertexLabels[i->vertex_iter(0)->id()];
		int vid1 = vertexLabels[i->vertex_iter(1)->id()];
		if (vid0 != vid1)
		{
			sharesMatrix[vid0][vid1] = sharesMatrix[vid1][vid0] += graphFeature.ef[i->id()].dualLength; //in dual graph case, use i->dualLength;
		}
	}

	//build adjacency and priority queue;
	std::multiset<ClusterArea> clusterQueue;
	for (int i = 0; i < scores.size(); i++)
	{
		clusterQueue.insert(ClusterArea(i, scores[i]));
	}

	//adjcents between clusters;
	std::list<std::multiset<DecreaseOrderExt> > clusterAdjcencyList(segNum);
	std::vector<std::list<std::multiset<DecreaseOrderExt> >::iterator> clusterAdjcency;
	for (auto i = clusterAdjcencyList.begin(); i != clusterAdjcencyList.end(); i++)
	{
		clusterAdjcency.push_back(i);
	}

	for (auto i = mesh->getEdges().begin(); i != mesh->getEdges().end(); i++)
	{
		int vid0 = vertexLabels[i->vertex_iter(0)->id()];
		int vid1 = vertexLabels[i->vertex_iter(1)->id()];
		if (vid0 != vid1)
		{
			if (clusterAdjcency[vid0]->find(DecreaseOrderExt(vid1, scores[vid1], sharesMatrix[vid0][vid1])) == clusterAdjcency[vid0]->end())
				clusterAdjcency[vid0]->insert(DecreaseOrderExt(vid1, scores[vid1], sharesMatrix[vid0][vid1]));

			if (clusterAdjcency[vid1]->find(DecreaseOrderExt(vid0, scores[vid0], sharesMatrix[vid0][vid1])) == clusterAdjcency[vid1]->end())
				clusterAdjcency[vid1]->insert(DecreaseOrderExt(vid0, scores[vid0], sharesMatrix[vid0][vid1]));
		}
	}

	//clusters
	std::list<std::set<int>> clusters(segNum);
	std::vector<std::list<std::set<int>>::iterator> iclusters;
	for (auto i = clusters.begin(); i != clusters.end(); i++)
	{
		iclusters.push_back(i);
	}
	for (int i = 0; i < vertexLabels.size(); i++)
	{
		iclusters[vertexLabels[i]]->insert(i);
	}

	while (!clusterQueue.empty())
	{
		ClusterArea ci = *clusterQueue.begin(); clusterQueue.erase(clusterQueue.begin());
		if (ci.m_val > thres)
		{
			break;
		}
		if (clusterAdjcency[ci.m_id]->empty())
		{
			continue;
		}

		DecreaseOrderExt cj_ = *clusterAdjcency[ci.m_id]->begin();
		ClusterArea cj(cj_.m_id, cj_.m_val);

		clusterAdjcency[ci.m_id]->erase(clusterAdjcency[ci.m_id]->begin()); //remove cj from ci's

		auto tp = clusterAdjcency[cj.m_id]->begin();
		for (; tp != clusterAdjcency[cj.m_id]->end(); tp++)
		{
			if (tp->m_id == ci.m_id)
				break;
		}
		if (tp == clusterAdjcency[cj.m_id]->end())
		{
			cout << "invalid merge" << endl;
		}
		else
		{
			clusterAdjcency[cj.m_id]->erase(tp);
		}

		auto i = std::find(clusterQueue.begin(), clusterQueue.end(), cj);
		if (i == clusterQueue.end())
		{
			cout << "invalid reference" << endl;
		}
		else
		{
			clusterQueue.erase(i);
		}

		//0. merge ci to cj.
		//1. cj is still in clusterQueue, need to update it, the new area is not just sum of ci and cj.
		//2. clear ci's clusterAdjacency, and add adj to cj, if not existed in cj. add cj to ci's adjacency, if not existed.

		//0. merge ci to cj.
		for (auto i = iclusters[ci.m_id]->begin(); i != iclusters[ci.m_id]->end(); i++)
		{
			vertexLabels[*i] = cj.m_id;
		}
		iclusters[cj.m_id]->insert(iclusters[ci.m_id]->begin(), iclusters[ci.m_id]->end());
		iclusters[ci.m_id]->clear();


		//1. cj is still in clusterQueue, need to update it, the new area is not just sum of ci and cj.
		double newScore;
		if (false)
		{//simple way, but incorrect, for instance, if ci,cj's area are both zero, then the sum is zero, however, the actual size is more than zero.
			newScore = ci.m_val + cj.m_val;
		}
			{
				newScore = 0;

				//I'll rewrite the part later. It doesn't need to revisit all triangles;
				for (auto i = mesh->getFaces().begin(); i != mesh->getFaces().end(); i++)
				{
					unsigned v1 = i->vertex_iter(0)->id();
					unsigned v2 = i->vertex_iter(1)->id();
					unsigned v3 = i->vertex_iter(2)->id();
					if (vertexLabels[v1] == vertexLabels[v2] && vertexLabels[v1] == vertexLabels[v3] && vertexLabels[v1] == cj.m_id)
					{
						newScore += faceAreas[i->id()];
					}
				}
			}
		clusterQueue.insert(ClusterArea(cj.m_id, newScore));

		//2. clear ci's clusterAdjacency, and add adj to cj, if not existed in cj.
		for (auto i = clusterAdjcency[ci.m_id]->begin(); i != clusterAdjcency[ci.m_id]->end(); i++)
		{
			int ci_a = i->m_id; // for each cluster

			//if ci_a not in cj, then add to cj, if exist, then update their shared boundary, m_val_ext;
			//if cj not in ci_a, then add to ci_a, if exist, then update their shared boundary, m_val_ext;
			auto ip = clusterAdjcency[ci_a]->begin();
			for (; ip != clusterAdjcency[ci_a]->end(); ip++)
			{
				if (ip->m_id == cj.m_id && ip->m_val == cj.m_val)
					break;
			}
			if (ip == clusterAdjcency[ci_a]->end())
			{
				clusterAdjcency[ci_a]->insert(DecreaseOrderExt(cj.m_id, newScore, i->m_val_ext));
			}
			else
			{
				double temp = ip->m_val_ext;
				clusterAdjcency[ci_a]->erase(ip);
				clusterAdjcency[ci_a]->insert(DecreaseOrderExt(cj.m_id, newScore, i->m_val_ext + temp));
			}

			ip = clusterAdjcency[cj.m_id]->begin();
			for (; ip != clusterAdjcency[cj.m_id]->end(); ip++)
			{
				if (ip->m_id == i->m_id && ip->m_val == i->m_val)
					break;
			}
			if (ip == clusterAdjcency[cj.m_id]->end())
			{
				clusterAdjcency[cj.m_id]->insert(DecreaseOrderExt(i->m_id, i->m_val, i->m_val_ext));
			}
			else
			{
				double temp = ip->m_val_ext;
				clusterAdjcency[cj.m_id]->erase(ip);
				clusterAdjcency[cj.m_id]->insert(DecreaseOrderExt(i->m_id, i->m_val, i->m_val_ext + temp));
			}


			//remove ci from ci_a
			ip = clusterAdjcency[ci_a]->begin();
			for (; ip != clusterAdjcency[ci_a]->end(); ip++)
			{
				if (ip->m_id == ci.m_id && ip->m_val == ci.m_val)
					break;
			}
			if (ip == clusterAdjcency[ci_a]->end())
			{
				cout << endl << ci_a << " error " << ci.m_id << " " << ci.m_val << endl;
			}
			else
			{
				clusterAdjcency[ci_a]->erase(ip);
			}
		}

		clusterAdjcency[ci.m_id]->clear();

		//for all cj_a in cj, even though not appeared in ci, need to update cj_a's element cj...
		for (auto i = clusterAdjcency[cj.m_id]->begin(); i != clusterAdjcency[cj.m_id]->end(); i++)
		{
			auto ip = clusterAdjcency[i->m_id]->begin();
			for (; ip != clusterAdjcency[i->m_id]->end(); ip++)
			{
				if (ip->m_id == cj.m_id)
					break;
			}
			if (ip == clusterAdjcency[i->m_id]->end())
			{
				cout << "bad reference" << endl;
			}
			else
			{
				clusterAdjcency[i->m_id]->erase(ip);
				clusterAdjcency[i->m_id]->insert(DecreaseOrderExt(cj.m_id, newScore, i->m_val_ext));
			}
		}
	}

	//update label;
	int ind = 0;
	for (auto i = clusters.begin(); i != clusters.end(); i++)
	{
		if (i->empty())
			continue;

		for (auto v = i->begin(); v != i->end(); v++)
		{
			vertexLabels[*v] = ind;
		}

		ind++;
	}
	segNum = ind;

	cout << endl << "cluster number after merging(by removing small patches):" << segNumber << endl;
	double cutCost = 0;
	for (auto i = myMesh->getEdges().begin(); i != myMesh->getEdges().end(); i++)
	{
		unsigned vs[] = { i->vertex_iter(0)->id(), i->vertex_iter(1)->id() };
		if (vertexLabel[vs[0]] != vertexLabel[vs[1]])
		{
			graphFeature.ef[i->id()].isCut = true;
			cutCost += mg.es[i->id()].w;
		}
		else
		{
			graphFeature.ef[i->id()].isCut = false;
		}
	}
	cout << "Cost:" << cutCost << endl;
}

void MeshSegment::getBoundaryOfClusters()
{
	for(auto e_it=myMesh->getEdges().begin();e_it!=myMesh->getEdges().end();e_it++)
	{
		unsigned vs[]={e_it->vertex_iter(0)->id(),e_it->vertex_iter(1)->id()};
		if(vertexLabel[vs[0]] != vertexLabel[vs[1]])
			graphFeature.ef[e_it->id()].isCut = true;
		else
			graphFeature.ef[e_it->id()].isCut = false;
	}
}

//smooth the curve network on 2d manifold domain
bool MeshSegment::vertexMoveToLocalMinima(const std::vector<BoundaryVertex>& Vbs, unsigned i, BoundaryVertex& newBv)
{
	//1. find the direction; 2. find the step size;
	//the energy function is sum of length of all vectors from source to neighbores, and the length is computed by our anis metric;
	//the direction is negative of gradient direction of energy function;
	//step is simply half of length from source to center of neighbores, or computed by lineSearch;

	const BoundaryVertex& vb = Vbs[i];
	Tensor& tvb = faceAnis[vb.fid];

	Vec3 grad(0, 0, 0); //first order derivative
	std::vector<Tensor> hsTensors; //second order derivative matrix
	//energy function is linear, process gradient of each neighbore and sum them up;
	for (unsigned j = 0; j < vb.ngbs.size(); j++)
	{
		const BoundaryVertex& nvb = Vbs[vb.ngbs[j]];
		Tensor& tnvb = faceAnis[nvb.fid];

		Vec3 v;
		v = vb.pos - nvb.pos;
		if (v.length() < 1e-4)
		{
// 			cout << "vertices are too close";
			continue;
		}

		if (false)
		{
			hsTensors.push_back(lineSearchHessian(v, tvb));
			hsTensors.push_back(lineSearchHessian(v, tnvb));
		}

		Vec2 vecProjL = Vec2(v.dot(tvb.dir1), v.dot(tvb.dir2));
		double lenL = sqrt(pow(vecProjL.x, 2)*tvb.mag1 + pow(vecProjL.y, 2)*tvb.mag2);
		vecProjL.x *= tvb.mag1; vecProjL.y *= tvb.mag2;
		lenL *= 2;
		if (lenL < 1e-4)
		{
			vecProjL = Vec2(0, 0);
		}
		else
		{
			vecProjL.x /= lenL;
			vecProjL.y /= lenL;
		}

		Vec2 vecProjR = Vec2(v.dot(tnvb.dir1), v.dot(tnvb.dir2));
		double lenR = sqrt(pow(vecProjR.x, 2)*tnvb.mag1 + pow(vecProjR.y, 2)*tnvb.mag2);
		vecProjR.x *= tnvb.mag1; vecProjR.y *= tnvb.mag2;
		lenR *= 2;
		if (lenR < 1e-4)
		{
			vecProjR = Vec2(0, 0);
		}
		else{
			vecProjR.x /= lenR;
			vecProjR.y /= lenR;
		}

		v = vecProjL.x*tvb.dir1 + vecProjL.y*tvb.dir2;
		grad += v;
		v = vecProjR.x*tnvb.dir1 + vecProjR.y*tnvb.dir2;
		grad += v;
	}

	Vec3 grad0 = grad;
	grad.normalize();
	grad = -grad;

	//step length
	Vec3  p(0, 0, 0);
	for (unsigned j = 0; j < vb.ngbs.size(); j++)
	{
		p += Vbs[vb.ngbs[j]].pos;
	}
	if (vb.ngbs.size() == 2)
	{
		p = p*0.5 - vb.pos;
	}
	else if (vb.ngbs.size() == 3)
	{
		p = p / 3.0 - vb.pos;
	}
	else
	{
		cout << "smooth boundary error" << endl;
	}
	double stepLen = p.length() /** 0.5*/;

	//step length by line search, all transport to tvb.dir1&tvb.dir2 frame;
	if (false)
	{
		Vec3 norm = tvb.dir1.cross(tvb.dir2); norm.normalize();
		Tensor hsTensor;
		for (int j = 0; j < hsTensors.size(); j++)
		{
			Tensor& f2 = hsTensors[j];
			Vec3 tNorm = f2.dir1.cross(f2.dir2); tNorm.normalize();
			f2.dir1 = basicNormalTransport(tNorm, norm, f2.dir1);
			f2.dir2 = basicNormalTransport(tNorm, norm, f2.dir2);

			Vec2 v1_, v2_;
			v1_.x = f2.dir1.dot(tvb.dir1); v1_.y = f2.dir1.dot(tvb.dir2);
			v2_.x = f2.dir2.dot(tvb.dir1); v2_.y = f2.dir2.dot(tvb.dir2);

			f2.mat[0][0] = f2.mag1*v1_.x*v1_.x + f2.mag2*v2_.x*v2_.x;
			f2.mat[0][1] = f2.mat[1][0] = f2.mag1*v1_.x*v1_.y + f2.mag2*v2_.x*v2_.y;
			f2.mat[1][1] = f2.mag1*v1_.y*v1_.y + f2.mag2*v2_.y*v2_.y;

			for (unsigned i = 0; i < 2; i++)
			{
				for (unsigned j = 0; j < 2; j++)
				{
					hsTensor.mat[i][j] += f2.mat[i][j];
				}
			}
		}
		Tensor::makeTensor(tvb.dir1, tvb.dir2, hsTensor);

		Vec2 f1 = Vec2(grad0.dot(hsTensor.dir1), grad0.dot(hsTensor.dir2));
		double denom = sqrt(pow(f1.x, 2)*hsTensor.mag1 + pow(f1.y, 2)*hsTensor.mag2);
		double stepLen0 = grad0.length() / denom;

		//if(stepLen0 > stepLen) cout<<"yes";
		stepLen = min(stepLen0, stepLen);
		//stepLen = stepLen0;
	}

	//new position
	p = vb.pos + grad * stepLen * 0.5;

	auto f_it = myMesh->getFIter()[vb.fid];
	std::set<unsigned> usedFaces;
	for (unsigned fn = 0; fn < 3; fn++)
	{
		auto vt = f_it->vertex_iter(fn);
		for (unsigned fv = 0; fv < vt->edge_iter().size(); fv++)
		{
			auto et = vt->edge_iter()[fv];
			usedFaces.insert(et->face_iter(0)->id());
			if (et->manifold())
			{
				usedFaces.insert(et->face_iter(1)->id());
			}
		}
	}

	//project new position back to mesh; we only look at a small neighbores;
	std::vector<double> disToTris;
	std::vector<Vec3> projPoints;
	for (auto fn = usedFaces.begin(); fn != usedFaces.end(); fn++)
	{
		auto f = myMesh->getFIter()[*fn];
		Vec3 projp;
		double d = IsPointInTriangle(f->vertex_iter(0)->coordinate(),
			f->vertex_iter(1)->coordinate(), f->vertex_iter(2)->coordinate(), p, projp);
		disToTris.push_back(d);
		projPoints.push_back(projp);
	}

	unsigned leastId = 0;;
	double leastDis = DBL_MAX;
	for (unsigned d = 0; d < disToTris.size(); d++)
	{
		if (disToTris[d] != -1 && disToTris[d] < leastDis)
		{
			leastDis = disToTris[d];
			leastId = d;
		}
	}

	auto it_face = usedFaces.begin();
	std::advance(it_face, leastId);
	if (disToTris[leastId] == -1 || projPoints[leastId] == Vec3(0, 0, 0) || *it_face > myMesh->getFaces().size() /*|| leastDis > 1e-3*/)
	{
		return false;
	}
	else	//if it is in the neighborhood
	{
		newBv.fid = *it_face;
		newBv.pos = projPoints[leastId];
		return true;
	}
}
void MeshSegment::boundarySmooth(bool doSmooth)
{
	if(vertexLabel.empty()) return;
	cout << endl << "smooth...";

	clock_t tstr = clock();

	std::vector<BoundaryVertex> Vbs;

	//init boundary network
	int newId = 0;
	std::vector<int> VbMap(myMesh->getFaces().size(),-1);
	for(auto e_it=myMesh->getEdges().begin();e_it!=myMesh->getEdges().end();e_it++)
	{
		if(graphFeature.ef[e_it->id()].isCut == false)continue;

		int vid1,vid2;
		if(e_it->manifold())
		{
			int fs[]={e_it->face_iter(0)->id(),e_it->face_iter(1)->id()};
			if(VbMap[fs[0]]==-1) 
			{
				VbMap[fs[0]] = newId; newId++;
				Vbs.push_back(BoundaryVertex());
				Vbs.back().fid = fs[0];
				Vbs.back().pos = graphFeature.cf[fs[0]].rep;
			}
			if(VbMap[fs[1]]==-1) 
			{
				VbMap[fs[1]] = newId; newId++;
				Vbs.push_back(BoundaryVertex());
				Vbs.back().fid = fs[1];
				Vbs.back().pos = graphFeature.cf[fs[1]].rep;
			}
			vid1 = VbMap[fs[0]];
			vid2 = VbMap[fs[1]];

			if (std::find(Vbs[vid1].ngbs.begin(), Vbs[vid1].ngbs.end(), vid2) == Vbs[vid1].ngbs.end())
				Vbs[vid1].ngbs.push_back(vid2);
			if (std::find(Vbs[vid2].ngbs.begin(), Vbs[vid2].ngbs.end(), vid1) == Vbs[vid2].ngbs.end())
				Vbs[vid2].ngbs.push_back(vid1);
		}
	}

	//smooth by anis tensor for several times.
	int smt = doSmooth ? boundarySmoothTimes : 0;
	for (unsigned cnt = smt; cnt > 0; cnt--)
	{
		std::vector<BoundaryVertex> newVbs = Vbs;
		for (unsigned i = 0; i < Vbs.size(); i++) //smooth for each vertex
		{
			if (Vbs[i].ngbs.size() == 1) continue; //boundary fixed;

			BoundaryVertex newBv;
			if (vertexMoveToLocalMinima(Vbs, i, newBv))
			{
				newVbs[i].fid = newBv.fid;
				newVbs[i].pos = newBv.pos;
			}
		}
		Vbs = newVbs;
	}

	//trace boundary curves from vertices and their adjacency.
	boundaryJoints.clear();
	for(int i=0;i<Vbs.size();i++)
	{
		if(Vbs[i].ngbs.size()>2) boundaryJoints.push_back(Vbs[i].pos);
	}
	std::vector<bool> vbUsed(Vbs.size(),false);
	boundaryCurves.clear();
	while(true)
	{
		unsigned vid = std::find(vbUsed.begin(),vbUsed.end(),false) - vbUsed.begin();
		if(vid>=vbUsed.size()) break;

		vbUsed[vid]=true;
		if(Vbs[vid].ngbs.size()!=2) continue;

		std::vector<int> boundFace;
		boundFace.push_back(vid);

		unsigned frontFace = Vbs[vid].ngbs.front();
		while(Vbs[frontFace].ngbs.size()==2 && !vbUsed[frontFace])
		{
			boundFace.push_back(frontFace);
			vbUsed[frontFace]=true;
			if(!vbUsed[Vbs[frontFace].ngbs.front()])
			{
				frontFace = Vbs[frontFace].ngbs.front();
			}
			else if(!vbUsed[Vbs[frontFace].ngbs.back()])
			{
				frontFace = Vbs[frontFace].ngbs.back();
			}
			else
			{
				if(Vbs[Vbs[frontFace].ngbs.front()].ngbs.size()==2)
					frontFace = Vbs[frontFace].ngbs.back();
				else
					frontFace = Vbs[frontFace].ngbs.front();
			}
		}
		boundFace.push_back(frontFace);
		vbUsed[frontFace]=true;

		reverse(boundFace.begin(),boundFace.end());

		frontFace = Vbs[vid].ngbs.back();
		while(Vbs[frontFace].ngbs.size()==2 && !vbUsed[frontFace])
		{
			boundFace.push_back(frontFace);
			vbUsed[frontFace]=true;
			if(!vbUsed[Vbs[frontFace].ngbs.front()])
			{
				frontFace = Vbs[frontFace].ngbs.front();
			}
			else if(!vbUsed[Vbs[frontFace].ngbs.back()])
			{
				frontFace = Vbs[frontFace].ngbs.back();
			}
			else{
				if(Vbs[Vbs[frontFace].ngbs.front()].ngbs.size()==2)
					frontFace = Vbs[frontFace].ngbs.back();
				else
					frontFace = Vbs[frontFace].ngbs.front();
			}
		}
		boundFace.push_back(frontFace);
		vbUsed[frontFace]=true;

		std::vector<Vec3> bds(boundFace.size());
		for( int i=0;i<bds.size();i++)
		{
			bds[i] = Vbs[boundFace[i]].pos;
		}
		boundaryCurves.push_back(bds);
	}

	//special cases: degenerate cases;
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		if (graphFeature.ef[e_it->id()].isCut == false)continue;
		int vid1, vid2;
		try{
			if (e_it->manifold()){
				int fs[] = { e_it->face_iter(0)->id(), e_it->face_iter(1)->id() };
				vid1 = VbMap[fs[0]];
				vid2 = VbMap[fs[1]];

				if (vid1>Vbs.size())
					throw vid1;
				else if (vid2 > Vbs.size())
					throw vid2;

				if (Vbs[vid1].ngbs.size() == 3 && Vbs[vid2].ngbs.size() == 3)
				{
// 					cout << "trigger special case.";
					std::vector<Vec3> bds(2);
					bds[0] = Vbs[vid1].pos;
					bds[1] = Vbs[vid2].pos;
					boundaryCurves.push_back(bds);
				}
			}
		}
		catch (int e)
		{
			cout << "Numeric issue" << e << endl;
		}
	}
	cout<<" done in "<<clock()-tstr<<"mm"<<endl;
}

//////////////////////////////////////////////////////////////////////////
//User Interaction of Segmentation
//Locally update of segmentation.
void MeshSegment::eraseFeatures(std::vector<unsigned>& fs)
{
	std::vector<bool> visFace(myMesh->getFaces().size(),true);
	for(int i=0;i<fs.size();i++)
	{
		visFace[fs[i]]=false;
	}
	for(int i=0;i<crestEdgesVisible.size();i++)
	{
		if(crestEdges[i][2]>myMesh->getFaces().size()) continue;
		if(!visFace[crestEdges[i][2]]) crestEdgesVisible[i] = false;
	}
	for(int i=0;i<userSketches.size();i++)
	{
		if(!visFace[userSketches[i].fid])
		{
			userSketches.erase(userSketches.begin()+i); i--;
		}
	}

// 	{
// 		featureGraphInitial();
// 		featureFromSketches();
// 		featureFromCrestline();
// 		computeFaceTensor(faceAnis);
// 		computeDualEdgeCostAndLength();
// 		computeDualEdgeStrength();
// 	}
}
void MeshSegment::eraseCluster(std::vector<unsigned>& fs)
{
	if(!vertexLabel.empty() && fs.size()>1)
	{
		cout << "erase cluster...";

		unsigned labelId = *std::max_element(vertexLabel.begin(), vertexLabel.end()) + 1;

		std::vector<std::set<int> > vadjs(labelId);
		for (unsigned i = 0; i < fs.size(); i++)
		{
			for (unsigned j = 0; j < 3; j++)
			{
				auto eit = myMesh->getFIter()[fs[i]]->edge_iter(j);
				int l1 = min(vertexLabel[eit->vertex_iter(0)->id()],vertexLabel[eit->vertex_iter(1)->id()]);
				int l2 = max(vertexLabel[eit->vertex_iter(0)->id()],vertexLabel[eit->vertex_iter(1)->id()]);
				if(l1 != l2)
				{
					vadjs[l1].insert(l2);
				}
			}			
		}

		{
			//reduce the strength of feature edge lying between merged patches.
			for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
			{
				int l1 = min(vertexLabel[e_it->vertex_iter(0)->id()], vertexLabel[e_it->vertex_iter(1)->id()]);
				int l2 = max(vertexLabel[e_it->vertex_iter(0)->id()], vertexLabel[e_it->vertex_iter(1)->id()]);

				if (!vadjs[l1].empty() && vadjs[l1].find(l2) != vadjs[l1].end())
				{
					if (graphFeature.ef[e_it->id()].lab != -1)
						edgeWeightParameter[e_it->id()] = 1e-20;//1e-20 a small number;
				}
			}
		}

		std::vector<std::list<int> > labels(labelId);
		for (unsigned i = 0; i < labelId; i++)
		{
			labels[i].push_back(i);
		}

		for (unsigned i = 0; i < vadjs.size(); i++)
		{
			std::vector<int> frontLabel;
			for(std::set<int>::iterator it = vadjs[i].begin();it != vadjs[i].end(); it++)
			{
				frontLabel.push_back(*it);
			}
			std::set<int> adjLabel;
			while(!frontLabel.empty())
			{
				int cl = frontLabel.back(); frontLabel.pop_back();
				adjLabel.insert(cl);
				for(std::set<int>::iterator it = vadjs[cl].begin();it != vadjs[cl].end(); it++)
				{
					frontLabel.push_back(*it);
				}
			}
			for(std::set<int>::iterator it = adjLabel.begin();it != adjLabel.end(); it++)
			{
				labels[i].push_back(*it);
				labels[*it].clear();
				vadjs[*it].clear();
			}
		}

		std::vector<int> newLabels(labels.size());
		int newLabelId=0;
		for (unsigned i = 0; i < labels.size(); i++)
		{
			if (labels[i].size() > 0)
			{
				for(std::list<int>::iterator itn = labels[i].begin(); itn!= labels[i].end(); itn++)
				{
					newLabels[*itn] = newLabelId;
				}
				newLabelId++;
			}
		}

		for(unsigned i=0;i<vertexLabel.size();i++)
		{
			vertexLabel[i] = newLabels[vertexLabel[i]];
		}
		
		cout<<"cluster number:"<<newLabelId;
		segNumber = newLabelId;
		getBoundaryOfClusters();
	}
}
void MeshSegment::addSkethFeatures(std::vector<MyInteraction::SketchFace>& fs)
{
	for(int i=0;i<fs.size();i++)
	{
		userSketches.push_back(fs[i]);
	}
}
void MeshSegment::findNeighboreOfTriangles(unsigned regionSize, std::vector<unsigned>& ToErase)
{
	if (regionSize == 0)return;

	std::vector<bool> visFace(myMesh->getFaces().size(), true);
	for (int i = 0; i < ToErase.size(); i++)
	{
		visFace[ToErase[i]] = false;
	}

	for (unsigned count = 1; count <= regionSize; count++)
	{
		unsigned fsize = ToErase.size();
		for (unsigned i = 0; i < fsize; i++)
		{
			auto f_it = myMesh->getFIter()[ToErase[i]];
			for (unsigned j = 0; j < 3; j++)
			{
				for (unsigned t = 0; t < f_it->vertex_iter(j)->edge_iter().size(); t++)
				{
					for (unsigned k = 0; k < 2; k++)
					{
						if (!f_it->vertex_iter(j)->edge_iter()[t]->manifold())
							break;
						unsigned fid = f_it->vertex_iter(j)->edge_iter()[t]->face_iter(k)->id();
						if (visFace[fid])
						{
							visFace[fid] = false; 
							ToErase.push_back(fid);
						}
					}

				}
			}
		}
	}
}
void  MeshSegment::smoothScribble(std::vector<unsigned>& sface,std::vector<Vec3>& scur)
{
// 	if (isSmoothStroke == false)
		return;

	std::map < unsigned, std::pair<unsigned,Vec3> > fmap;
	for (unsigned i = 0; i < sface.size(); i++)
	{
		fmap[sface[i]] = std::pair<unsigned, Vec3>(i,scur[i]);
	}
	for (unsigned i = 0; i < fmap.size(); i++)
	{
		auto tf = fmap.begin(); std::advance(tf, i);
		fmap[tf->first] = std::pair<unsigned, Vec3>(i, tf->second.second);
	}

	std::vector<BoundaryVertex> Vbs(fmap.size());
	std::vector<int> VbMap(myMesh->getFaces().size(), -1);
	for (auto it_f = fmap.begin(); it_f != fmap.end(); it_f++)
	{
		unsigned tind = it_f->second.first;
		VbMap[it_f->first] = tind;
		Vbs[tind].fid = it_f->first;
		Vbs[tind].pos = it_f->second.second;		
	}

	for (unsigned i = 1; i < sface.size(); i++)
	{
		if (sface[i - 1] == sface[i]) continue;

		unsigned vid1 = fmap[sface[i - 1]].first;
		unsigned vid2 = fmap[sface[i]].first;

		Vbs[vid1].ngbs.push_back(vid2);
		Vbs[vid2].ngbs.push_back(vid1);
	}

	//smooth by anis tensor for several times.
	const auto& fts = myMesh->getFIter();
	for (int cnt = boundarySmoothTimes; cnt > 0; cnt--)
	{
		std::vector<BoundaryVertex> newVbs = Vbs;
		for (int i = 0; i < Vbs.size(); i++) //smooth for each vertex
		{
			if (Vbs[i].ngbs.size() == 1) continue;

			BoundaryVertex newBv;
			if (vertexMoveToLocalMinima(Vbs, i, newBv))
			{
				if ((newVbs[i].pos - newBv.pos).length() > 1e-3 || Vbs[i].fid != newBv.fid)
				{
					newVbs[i].fid = newBv.fid;
					newVbs[i].pos = newBv.pos;
				}
			}
		}
		Vbs = newVbs;
	}

	//write to sketches;
	std::vector<unsigned> newfaces;
	std::vector<Vec3> newVers;
	for (unsigned i = 0; i < sface.size(); i++)
	{
		newfaces.push_back(Vbs[VbMap[sface[i]]].fid);
		newVers.push_back(Vbs[VbMap[sface[i]]].pos);
	}
	sface.swap(newfaces);
	scur.swap(newVers);
	cout << "scribble smoothed" << endl;
}
void MeshSegment::OverSegmentationOnPartialMesh(std::vector<unsigned>& fs, std::vector<unsigned>& es)
{
	cout << "------------------ PARTITIAL PARTITION -----------" << endl;
	//1. init GraphFeature data structure;
	{
		clock_t tstr = clock();
		featureGraphInitial();
		featureFromCrestline();
		featureFromSketches();

		computeFaceTensor(faceAnis);
		computeDualEdgeCostAndLength();
		computeDualEdgeStrength();

		featureToCurves();

		mgraphInit(mg);

		clock_t tinit = clock() - tstr;
		cout << "alpha:" << graphFeature.alpha << endl;
		cout << endl << "graph initialization:" << tinit / 1000 << "sec" << tinit % 1000 << "mm" << endl;
	}

	Mitani_Watershed_Dual_Partial(fs);

// 	genPatchColor(patchColors, segNumber, false);
// 	getBoundaryOfClusters();
}
void MeshSegment::regionMergingOnPartialMesh(MGraph& gph,
	std::vector<bool>&subgraphVers,
	std::vector<unsigned>& vertexLabel,
	unsigned& labelId)
{
	clock_t tstr = clock();
	cout << "perform in subgraph" << endl;

	//find patches corresponding to subgraph
	//get map between global index(whole mesh) and local index(subgraph), including vertex index and label index;
	std::set<unsigned> patchLabels;
	std::map<unsigned, unsigned> globalIndToLocalInd;
	unsigned num = 0;
	for (unsigned i = 0; i < subgraphVers.size(); i++)
	{
		if (subgraphVers[i])
		{
			patchLabels.insert(vertexLabel[i]);
			globalIndToLocalInd[i] = num;
			num++;
		}
	}
	std::map<unsigned, unsigned> globalLabelToLocalLabel;
	unsigned ind = 0;
	for (auto i = patchLabels.begin(); i != patchLabels.end(); i++, ind++)
	{
		globalLabelToLocalLabel[*i] = ind;
	}
	std::vector<unsigned> localLabels(num); //local index, and local label
	for (auto i = globalIndToLocalInd.begin(); i != globalIndToLocalInd.end(); i++)
	{
		localLabels[i->second] = globalLabelToLocalLabel[vertexLabel[i->first]];
	}

	//local adjacency of patches, the pair indices with cost;
	labelId = patchLabels.size();
	std::vector<std::set<int> > adjVers(labelId);
	for (int i = 0; i < (int)gph.es.size(); i++)
	{
		if (!subgraphVers[gph.es[i].n1] || !subgraphVers[gph.es[i].n2]) continue;
		
		unsigned vid0 = min(globalLabelToLocalLabel[vertexLabel[gph.es[i].n1]], globalLabelToLocalLabel[vertexLabel[gph.es[i].n2]]);
		unsigned vid1 = max(globalLabelToLocalLabel[vertexLabel[gph.es[i].n1]], globalLabelToLocalLabel[vertexLabel[gph.es[i].n2]]);

		if (vid0 == vid1)continue;
		adjVers[vid0].insert(vid1);
	}

	std::list< std::list< MGTriple > > pairClusterCost(labelId);
	typedef std::list< std::list< MGTriple > >::iterator PCIter;
	std::vector< PCIter > pCIter;

	//cost between two patches is the sum of weight of boundary edges in between.
	int count = 0;
	MGTriple tmg;
	for (PCIter pct = pairClusterCost.begin(); pct != pairClusterCost.end(); pct++, count++)
	{
		if (adjVers[count].empty())
		{
			tmg.j = count;
			tmg.val = 0;
			pct->push_back(tmg);
		}
		else
		{
			for (std::set<int>::iterator ita = adjVers[count].begin(); ita != adjVers[count].end(); ita++)
			{
				tmg.j = *ita;
				tmg.val = 0;
				pct->push_back(tmg);
			}
		}
		pCIter.push_back(pct);
	}

	std::list< MGTriple >::iterator tpc;
	for (int i = 0; i < (int)gph.es.size(); i++)
	{
		if (!subgraphVers[gph.es[i].n1] || !subgraphVers[gph.es[i].n2]) continue;

		unsigned vid0 = min(globalLabelToLocalLabel[vertexLabel[gph.es[i].n1]], globalLabelToLocalLabel[vertexLabel[gph.es[i].n2]]);
		unsigned vid1 = max(globalLabelToLocalLabel[vertexLabel[gph.es[i].n1]], globalLabelToLocalLabel[vertexLabel[gph.es[i].n2]]);

		if (vid0 == vid1)continue;

		for (tpc = pCIter[vid0]->begin(); tpc != pCIter[vid0]->end(); tpc++)
		{
			if (tpc->j == vid1) break;
		}
		tpc->val += gph.es[i].w;
	}

	//merging algorithms; Interfaces 
	std::vector<unsigned> newLabels;
	LMP_Merging(localLabels, pCIter, labelId, newLabels);

	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (!e->manifold()) continue;

		unsigned n1 = e->vertex_iter(0)->id();
		unsigned n2 = e->vertex_iter(1)->id();
		if (subgraphVers[n1] && subgraphVers[n2])
		{

			if (newLabels[globalIndToLocalInd[n1]] != newLabels[globalIndToLocalInd[n2]])
			{
				gph.es[e->id()].isCut = true;
			}
			else
			{
				gph.es[e->id()].isCut = false;
			}
		}
		else
		{
			if (vertexLabel[e->vertex_iter(0)->id()] != vertexLabel[e->vertex_iter(1)->id()])
				gph.es[e->id()].isCut = true;
			else
				gph.es[e->id()].isCut = false;
		}
	}

	unsigned vNum = myMesh->getVertices().size();
	unsigned eNum = myMesh->getEdges().size();
	vertexLabel.clear(); vertexLabel.resize(gph.vNum);
	labelId = 0;
	std::vector<bool> visitedVer(vNum, false);
	std::vector<bool> visitedEdge(eNum, false);
	for (size_t v = 0; v < vNum; v++)
	{
		if (visitedVer[v] == true) continue;

		std::vector<int> frontVer(1, v);
		//propagates vertex with same label, and group them into a subgraph
		while (!frontVer.empty())
		{
			int fv = frontVer.back(); frontVer.pop_back();
			vertexLabel[fv] = labelId;

			auto v = gph.vs[fv].vadjs.begin();
			auto e = gph.vs[fv].eadjs.begin();
			for (; v != mg.vs[fv].vadjs.end(); v++, e++)
			{
				if (gph.es[*e].isCut || visitedEdge[*e] == true)
					continue;

				visitedEdge[*e] = true;

				if (visitedVer[*v] == true)
					continue;

				visitedVer[*v] = true;

				frontVer.push_back(*v);
			}
		}
		labelId++;
	}

	cout << endl << "cluster number after merging:" << labelId << endl;

	clock_t tinit = clock() - tstr;
	cout << "merging takes:" << tinit / 1000 << "sec" << tinit % 1000 << "mm" << endl;
}
void MeshSegment::mergePartitionOnPartialMesh()
{
	boundaryCurves.clear();

	cout << "*   *   *   *   *   * PARTIAL MERGE  *   *   *   *   *   *" << endl;
	auto tlabels = vertexLabel;
	auto patNum = segNumber;

	vertexLabel = vertexLabelBeforeMerge;

	mgraphInit(mg);

	regionMergingOnPartialMesh(mg, subgraphVers, vertexLabel, segNumber);

	if (isMergeSmallPatch)
	{
		mergeSmallPatches();
	}
	else
	{
		for (int i = 0; i < mg.eNum; i++)
		{
			graphFeature.ef[i].isCut = mg.es[i].isCut;
		}
	}

	if (autoBoundarySmooth)
	{
		cout << endl;
		boundarySmooth(true);
	}

	cout << "\n-----------------------------------------------------" << endl;

	//label id for most of vertices should remain the same, causing least changes;
// 	patchColorMinChanged(tlabels, patNum);
}
void MeshSegment::mergeLocalClusterByFeatureEnhancing(std::vector<unsigned>& ToEnhance)
{
	const auto& faceIters = myMesh->getFIter();
	const auto& verIters = myMesh->getVIter();

	//1. determine which clusters does paints lie in
	std::set<int> inClusters;
	for (int i = 0; i < ToEnhance.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			auto vit = faceIters[ToEnhance[i]]->vertex_iter(j);
			inClusters.insert(vertexLabel[vit->id()]);
		}
	}

	//2. find out clusters(before merging) in clusters by step 1.
	//for all vertices in clusters by step 1, find their initial clusters
	std::vector<bool> includedClusters(segNumber, false);
	for (auto i = inClusters.begin(); i != inClusters.end(); i++)
	{
		includedClusters[*i] = true;
	}
	std::map<int, int> mapToNewClusters;
	std::vector<bool> includedVers(verIters.size(), false);
	for (int i = 0; i < verIters.size(); i++)
	{
		if (includedClusters[vertexLabel[i]])
		{
			mapToNewClusters[vertexLabelInit[i]] = vertexLabel[i];
			includedVers[i] = true;
		}
	}

	subgraphVers = includedVers;
	//assign initial vertexlabel to seleted regions, and keep non-seleted the same; update vertexLabel, and call merging
	{//merge locally
		int clusterSize = inClusters.size();
		int clusterSizeOrig = mapToNewClusters.size();
		std::set<int> newLabels;
		for (auto i = inClusters.begin(); i != inClusters.end(); i++)
		{
			newLabels.insert(*i);
		}
		for (int i = 0; i < clusterSizeOrig - clusterSize; i++)
		{
			newLabels.insert(segNumber + i);
		}

		for (int i = 0; i < includedVers.size(); i++)
		{
			if (includedVers[i])
			{
				int labelId = *newLabels.begin(); newLabels.erase(newLabels.begin());
				std::vector<MyMesh::VertexIter> frontVer(1, verIters[i]);
				//propagates vertex with same label, and group them into a subgraph
				while (!frontVer.empty())
				{
					auto fv = frontVer.back(); frontVer.pop_back();
					vertexLabel[fv->id()] = labelId;

					for (int j = 0; j < fv->vertex_iter().size(); j++)
					{
						auto av = fv->vertex_iter()[j];
						if (vertexLabelInit[fv->id()] == vertexLabelInit[av->id()] && includedVers[av->id()])
						{
							frontVer.push_back(av);
							includedVers[av->id()] = false;
						}
					}
				}
			}
		}
		if (!newLabels.empty())
			cout << "something is wrong when assigning new labels when painting" << endl;

		vertexLabelBeforeMerge.swap(vertexLabel);
	}
	//3. change local merging parameters by paint setting.. and computer the graph weight, then perform merge in clusters by step 1.
	updateLocalAwardAlpha(ToEnhance);

	mergePartitionOnPartialMesh();
}
void MeshSegment::modifySegmentBySketch(Interaction::SCRIBBLE_TYPE type)
{
	if(sketchFaces.size()<2) return; //return by empty;
	clock_t tstr = clock();

	std::vector<unsigned> tlabels = vertexLabel;
	unsigned patNum = segNumber;

	std::vector<unsigned> ToErase = sketchFaces; //Erase features in 'sketchFaces';

	if (type == Interaction::ADD)
	{
		std::vector<unsigned> edgeSequence;//edge sequence along unfolded triangle strip of scribbles;
		std::vector<MyInteraction::SketchFace> ToAdd; //sequence of triangles;
		m_interaction->userSketchSamplingFromSketchFaces(sketchFaces,sketchCurves,ToAdd, edgeSequence);//find sketch curve from unordered sketch faces;
		
		if (edgeSequence.empty()) return;

		if (!graphCutLocally)
		{
			unsigned ringSize = strokeSize;
			findNeighboreOfTriangles(ringSize, ToErase); //includes the neighbores of sketches, erase a larger region
			// 		auto tp = paintParamType;		paintParamType = 2;
			// 		updateLocalAwardAlpha(ToErase);	paintParamType = tp;
			eraseFeatures(ToErase);

			addSkethFeatures(ToAdd);//add new sketches;
		}

		OverSegmentationOnPartialMesh(sketchFaces,edgeSequence);//do over segmentation;
		mergePartitionOnPartialMesh();//do merging in seleted patches;
	}
	else if (type == Interaction::ERASE)
	{
		eraseCluster(ToErase);
		if (autoBoundarySmooth)	boundarySmooth(true);
	}
	else
	{
		mergeLocalClusterByFeatureEnhancing(ToErase);//do merging in seleted patches;
	}

	patchColorMinChanged(tlabels, patNum);//label id for most of vertices should remain the same;

	sketchFaces.clear(); sketchCurves.clear();

	clock_t totalTime=clock()-tstr;
	cout<<"total time:"<<totalTime/1000<<"sec"<<totalTime%1000<<"mm"<<endl;
}