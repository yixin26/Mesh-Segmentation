/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#ifndef ANIS_VISUAL_H
#define ANIS_VISUAL_H

#include "myMesh.h"
#include <vector>
#include "amlVec.h"
#include "colorEngine.h"

class AnisGeodesic;

class AnisVisual{
public:
	AnisVisual(){
		showMesh=false;anisGeodesy=NULL;isDoSrcToAll=showDistanceField=showAllPaths=showIsoline=false;
	}

	void genDisplayList(){
		wireMesh=++count;solidMesh=++count; tensorLineMax=++count;tensorLineMin=++count;tensorEllipse=++count;
		distanceField=++count;allPaths=++count;isoline=++count;
	}
	void drawMesh();
	void drawSolidMesh(); void solidMeshList();
	void drawWireMesh(double w = 1.0, Vec3 col = Vec3(0.7, 0.7, 0.7));  void wireMeshList();

	void drawDirection(unsigned dir);	void tensorLineList(unsigned dir);
	void drawEllipse();					void tensorEllipseList();

	bool drawDistanceField();	void distanceFieldList();
	void drawIsoline();			void isoLineList();
	void drawAllPaths();		void allPathsList();

	void genDistanceFieldColor();
	void genIsoLine();
	void genDisFieldRank();


	void drawBadFaces();

public:
	AnisGeodesic* anisGeodesy;

	bool showMesh;

	static unsigned count;
	unsigned wireMesh,solidMesh;
	unsigned tensorLineMax,tensorLineMin,tensorEllipse;
	unsigned distanceField,allPaths,isoline;

	static unsigned lineWidth,shapeSize,tensorShapeType;

	static bool isDoSrcToAll,showDistanceField,showIsoline,showAllPaths;
	std::vector<std::vector<Vec3> > srcToAllPaths;
	std::vector<double> srcToAllPathsLength;
	std::vector<std::vector<Vec3> > isoLines;

	static bool showInvalidTriangle;

private:
	void ellipse(const MyMesh::VertexIter& v_it,std::vector<double>& sin360,std::vector<double>& cos360);
	void tensorShape(const MyMesh::VertexIter& v_it,std::vector<double>& sin360,std::vector<double>& cos360);

	std::vector<ColorEngine::color> distanceFieldColor;
	std::vector<unsigned>  disFieldMagnitudeRank;
};

#endif