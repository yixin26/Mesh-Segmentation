/***************************************************************
* Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
**************************************************************/

#include "core/Segmentation.h"
#include <fstream>
#include <iomanip>
#include "TriMesh.h"
#include "TriMesh_algo.h"

using namespace std;

void MeshSegment::rescaleMesh()
{
	cout << "mesh normalized to unit cube" << endl;

	MyMesh::BBox m_bbox;

	m_bbox.bBoxMin = Vec3(1e10, 1e10, 1e10);
	m_bbox.bBoxMax = Vec3(-1e10, -1e10, -1e10);
	m_bbox.centroid = Vec3(0, 0, 0);

	unsigned int numPoints = 0;

	for (unsigned i = 0; i < points.size(); i++)
	{
		auto &p = points[i];
		for (unsigned int d = 0; d < 3; d++)
		{
			if (p[d] > m_bbox.bBoxMax[d])
				m_bbox.bBoxMax[d] = p[d];
			if (p[d] < m_bbox.bBoxMin[d])
				m_bbox.bBoxMin[d] = p[d];
		}
		m_bbox.centroid += p;
		numPoints++;
	}

	m_bbox.centroid *= 1.0 / numPoints;

	//find the sides of bbox
	Vec3 bBox;
	for (unsigned int d = 0; d < 3; d++)
	{
		bBox[d] = abs(m_bbox.bBoxMin[d] - m_bbox.bBoxMax[d]);
	}

	//find the largest side
	double largestSide = (bBox[0] > bBox[1]) ? bBox[0] : bBox[1];
	largestSide = (largestSide > bBox[2]) ? largestSide : bBox[2];

	//translate the points so centroid is at origin
	//scale the points so bbox fits in canonical cube
	for (unsigned i = 0; i < points.size(); i++)
	{
		auto &p = points[i];
		p = p - m_bbox.centroid;
		//scale the points so bbox fits in canonical cube
		p *= 1.0 / largestSide;
	}

	m_bbox.isInit = true;
}
void MeshSegment::readTrimesh(const char* fileName)
{
	triMesh = TriMesh::read(fileName);

	points.clear();	faces.clear();
	int vNum = triMesh->vertices.size();
	points.reserve(vNum);
	for (int i = 0; i < vNum; i++)
	{
		points.push_back(Vec3(triMesh->vertices[i][0], triMesh->vertices[i][1], triMesh->vertices[i][2]));
	}
	int fNum = triMesh->faces.size();
	faces.reserve(fNum);
	std::vector<unsigned> tf(3);
	for (int i = 0; i < fNum; i++)
	{
		tf[0] = triMesh->faces[i][0]; tf[1] = triMesh->faces[i][1]; tf[2] = triMesh->faces[i][2];
		faces.push_back(tf);
	}

	rescaleMesh();

	for (int i = 0; i < vNum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			triMesh->vertices[i][j] = points[i][j];
		}
	}

	triMesh->need_normals();

	myMesh = new MyMesh;
	myMesh->initialize(points, faces);
	computeCurvature();
}

void MeshSegment::saveSegmentation()
{
	QString ClusterFile = m_filename;
	ClusterFile.replace(ClusterFile.lastIndexOf(".") + 1, 3, "label");

	std::ofstream File;
	File.open(ClusterFile.toStdString(), std::ios::out);
	if (!File.good())
		return;

	//store vertex labels;
	File << "[num of vertices] [num of clusters]" << endl;
	File << myMesh->getVertices().size() << ' ' << *std::max_element(vertexLabel.begin(), vertexLabel.end()) + 1 << endl;
	File << "[index of vertex] [label of vertex]" << endl;
	for (MyMesh::VertexIter v_it = myMesh->getVertices().begin(); v_it != myMesh->getVertices().end(); v_it++)
	{
		File << "v" << v_it->id() << ' ' << vertexLabel[v_it->id()] << endl;
	}
	File.close();

	//store original boundary
	boundarySmooth(false);

	QString CurveFile = m_filename;
	CurveFile.replace(CurveFile.lastIndexOf(".") + 1, 3, "curve");

	File.open(CurveFile.toStdString(), std::ios::out);
	if (!File.good())
		return;

	//store vertex labels;
	File << "[num of curves]" << endl;
	File << boundaryCurves.size() << endl;

	for (int i = 0; i < boundaryCurves.size(); i++)
	{
		File << "[length of current curve]" << endl;
		File << boundaryCurves[i].size() << endl;
		for (int j = 0; j < boundaryCurves[i].size(); j++)
		{
			File << boundaryCurves[i][j].x << ' ' << boundaryCurves[i][j].y << ' ' << boundaryCurves[i][j].z << endl;
		}
	}
	File.close();

	//store smoothed boundary
	boundarySmooth(true);

	CurveFile = m_filename;
	CurveFile.replace(CurveFile.lastIndexOf("."), 13, "_smooth.curve");

	File.open(CurveFile.toStdString(), std::ios::out);
	if (!File.good())
		return;

	//store vertex labels;
	File << "[num of curves]" << endl;
	File << boundaryCurves.size() << endl;

	for (int i = 0; i < boundaryCurves.size(); i++)
	{
		File << "[length of current curve]" << endl;
		File << boundaryCurves[i].size() << endl;
		for (int j = 0; j < boundaryCurves[i].size(); j++)
		{
			File << boundaryCurves[i][j].x << ' ' << boundaryCurves[i][j].y << ' ' << boundaryCurves[i][j].z << endl;
		}
	}
	File.close();

}
