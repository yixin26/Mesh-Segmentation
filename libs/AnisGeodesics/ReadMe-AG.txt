Anisotrpic Geodesics

The implementation uses Kirsanov's code(https://code.google.com/p/geodesic/).
 
///////////////////////////////////////////////////

#include "AnisGeodesic.h"
#include "AnisVisual.h"


int main()
{
	//there are different type of anisotropic setting; 
	//use OurMax as example;
	AnisGeodesic* ourGeodesy2= new AnisGeodesic;

	//initiazation of mesh;
	MyMesh myMesh=new MyMesh;
	myMesh->initialize(points,faces);
	ourGeodesy2>orgMesh = myMesh;

	ourGeodesy2>scaleParameter=0.1;

	computeCurvature();

	//init field;
	ourGeodesy2->anisInit_OurMin(ourGeodesy2->metricVariation);ourGeodesy2->meshEmbedding_Tensor();

	//init geodesic
	ourGeodesy2->geodesyInit();

	//select a vertex, and propagate the geodesics
	ourGeodesy2->geodesyPropagation(s.type,s.id,s.pos);

	//trace back the geodesy between two ends;
	ourGeodesy2->geodesyTraceBack(t.type,t.id,t.pos);

	return 0;
}


