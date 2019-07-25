/***************************************************************
* Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
**************************************************************/

/*    Import note
  Mouse     position (0,0) is on the left top;
  GL buffer position (0,0) is on the left bottom;

  glGetIntegerv(GL_VIEWPORT, viewport);
  this func returns the size of framework, not the region of gl widget.!!
  So, use widget.geometry() instead! 

*/


#ifndef MY_INTERACTION
#define MY_INTERACTION

#include <GL/glu.h>
//#include <glut.h>
#include <iostream>
#include <set>
#include <map>
#include "myMesh.h"
#include <QtCore/QString>

namespace Interaction
{
	enum USER_INPUT_TYPE{ GEODESICS, SKETCH };
	enum SCRIBBLE_TYPE{	REVEAL, CONCEAL, ADD, ERASE};
	enum SELECT_REGION{	PIXEL, RECTANGLE};
}

class MyInteraction
{
public:

	struct Hit
	{
		unsigned hittedTriangle;
		unsigned hittedNode; //the node of triangle;
		Vec3 hittedVertex; //intersection on triangle;
		Vec3 hittedNearVertex;
	};

	struct SketchFace
	{
		int fid;
		int e1, e2;
		Vec3 fpos;
		Vec3 ep1, ep2;
	};

	MyInteraction()
	{
		init();
	};
	~MyInteraction(){};
	void clearSketches()
	{
		scribbleCurve.clear(); scribbleCurves.clear(); scribbleTriangleStrip.clear(); scribbleTriangleStrips.clear();
	}
	void init()
	{
		mesh = NULL; gl_buffer = NULL;	rgbBuffer = NULL; sr = Interaction::PIXEL;  user_input_type = Interaction::SKETCH; scribbleType = Interaction::REVEAL;
		glObjColors.clear(); windowSize_x = windowSize_y = 0;	scribleRadius = 0; brushSize = 20;
		hit.hittedVertex = hit.hittedNearVertex = AML::double3(0, 0, 0); hit.hittedNode = 0;	hit.hittedTriangle = 0;
		isControlPress = isAltPress = isShiftPress = isScrolling = false;
		clearSketches();
	}
	MyMesh** getMyMesh(){ return &mesh; }

	void drawTargets();
	void vertexLineIntersection(const std::pair<Vec3, Vec3>& Line, const Vec3& P, Vec3&I, double& dist);
	double IsPointInTriangle(Vec3 tp1, Vec3 tp2, Vec3 tp3, Vec3 testp, Vec3& projp);
	int triangleLineIntersection(const std::pair<Vec3, Vec3>& R, const std::vector<Vec3>& T, Vec3 &I=Vec3(0,0,0));
	void getRay(unsigned mouseX, unsigned mouseY, std::pair<Vec3, Vec3> &Ray);
	bool glSelection1(unsigned mouse_x, unsigned mouse_y, unsigned& hittedTriangle);
	void initSelectionColor();
	void updateBackBuffer(int,int,int,int);
	bool glSelection2(unsigned mouse_x, unsigned mouse_y, unsigned& hittedTriangle);
	bool glSelection2_Rectangle(unsigned mouse_x, unsigned mouse_y, unsigned width, unsigned height, std::set<unsigned>& hittedTriangles);
	unsigned getClosestVertex(unsigned mouse_x, unsigned mouse_y, unsigned hittedTriangle);
	double getClosestVertex(unsigned mouse_x, unsigned mouse_y, const std::vector<std::vector<Vec3>>& curves, std::pair<unsigned, unsigned>& s);
	void genStraightScreenLine(unsigned str_x, unsigned str_y, unsigned end_x, unsigned end_y, std::vector<std::pair<unsigned, unsigned>>& sequence);
	void connectTwoTriangles(Vec3& p, Vec3& norm, Vec3& dir, unsigned f1, unsigned f2, std::vector<Vec3>& pts, std::vector<unsigned>& fs);
	bool scribleOnMesh(unsigned mouseX, unsigned mouseY, Hit& lasthit, Hit& hit, bool isSingleHit = true, std::vector<Hit>& hits = std::vector<Hit>(0));
	bool brushOnMesh(unsigned mouseX, unsigned mouseY, unsigned nbsize, std::set<unsigned>& addToTriangles, Hit& hit); //scribbling with a radius;

	void pushToMultipleScribble();
	void popLastDraw(){ scribbleCurves.pop_back();scribbleTriangleStrips.pop_back(); }

	//post-process;
	void userSketchSamplingFromSketchFaces(std::vector<unsigned>& sketchFaces, std::vector<Vec3>& sketchCurves, std::vector<SketchFace>&, std::vector<unsigned>&);

	void processSketchForSeg(std::vector<unsigned>& sketchFaces, std::vector<Vec3>& sketchCurves);

private:

	double pointToPlane(Vec3& p, Vec3& pl_n, Vec3& pl_p);

private:

	MyMesh* mesh; //the shape where interaction works on;

	GLuint	*gl_buffer; //selection 1

	int windowSize_x, windowSize_y; //frame size
	std::vector<AML::ubyte3> glObjColors;//selection 2
	unsigned char* rgbBuffer; //selection 2, read back buffer;	

public:

	Interaction::USER_INPUT_TYPE user_input_type;
	Interaction::SCRIBBLE_TYPE scribbleType;
	Interaction::SELECT_REGION sr; //type of selection; a vertex or a region;
	
	Hit hit; //location of intersection, the triangle, the vertex in triangle and the vertex from near view point.

	std::vector<Vec3> scribbleCurve; // a sequeence of hits on the mesh;
	std::vector<std::vector<Vec3> > scribbleCurves;
	std::vector<unsigned> scribbleTriangleStrip;
	std::vector<std::vector<unsigned> > scribbleTriangleStrips;

	unsigned brushSize; //radius of pixels
	double scribleRadius;//radius of 3d object, for visualization;

	//hot keys;
	bool isControlPress, isAltPress, isShiftPress, isScrolling;
};
inline void MyInteraction::drawTargets()											// Draws The Targets (Needs To Be Seperate)
{
	for (auto f = mesh->getFaces().begin(); f != mesh->getFaces().end(); f++)
	{
		glLoadName(f->id());										// Assign Object A Name (ID)
		glBegin(GL_TRIANGLES);
		for (int j = 0; j < 3; j++)
		{
			auto& v = f->vertex_iter(j)->coordinate();
			glVertex3d(v.x, v.y, v.z);
		}
		glEnd();
	}

}
inline void MyInteraction::vertexLineIntersection(const std::pair<Vec3, Vec3>& Line, const Vec3& P, Vec3&I, double& dist)
{
	Vec3 vecLine = Line.first - Line.second; vecLine.normalize();

	Vec3 vp = P - Line.first;
	double s1 = abs(vecLine.dot(vp));

	vp = P - Line.second;
	double s2 = abs(vecLine.dot(vp));

	double l = s1 + s2;
	s1 /= l;
	s2 /= l;

	I = Line.first*s2 + Line.second*s1;
	dist = (P - I).length();
}

inline double MyInteraction::IsPointInTriangle(Vec3 tp1, Vec3 tp2, Vec3 tp3, Vec3 testp, Vec3& projp)
{
	Vec3 ax = tp2 - tp1; ax.normalize();
	Vec3 ay = tp3 - tp1;
	Vec3 nm = ax.cross(ay); nm.normalize();
	ay = nm.cross(ax); ay.normalize();

	std::pair<Vec3, Vec3> Ray;
	Ray.first = testp + nm;
	Ray.second = testp - nm;

	std::vector<Vec3> Tri(3);
	Tri[0] = tp1;Tri[1] = tp2;Tri[2] = tp3;
	if (triangleLineIntersection(Ray, Tri, projp) == 1)
		return (testp-projp).length();
	else
		return -1;
}
inline int MyInteraction::triangleLineIntersection(const std::pair<Vec3, Vec3>& R, const std::vector<Vec3>& T, Vec3 &I)
{
	Vec3    u, v, n;              // triangle vectors
	Vec3    dir, w0, w;           // ray vectors
	double     r, a, b;              // params to calc ray-plane intersect

	// get triangle edge vectors and plane normal
	u = T[1] - T[0];
	v = T[2] - T[0];
	n = u.cross(v);              // cross product
	if (n == Vec3(0, 0, 0))             // triangle is degenerate
		return -1;                  // do not deal with this case

	dir = R.second - R.first;              // ray direction vector
	w0 = R.first - T[0];
	a = -n.dot(w0);
	b = n.dot(dir);
	if (fabs(b) < 0.000001)
	{     // ray is  parallel to triangle plane
		if (a == 0)                 // ray lies in triangle plane
			return 2;
		else return 0;              // ray disjoint from plane
	}

	// get intersect point of ray with triangle plane
	r = a / b;
	if (r < 0.0)                    // ray goes away from triangle
		return 0;                   // => no intersect
	// for a segment, also test if (r > 1.0) => no intersect

	I = R.first + r * dir;            // intersect point of ray and plane

	// is I inside T?
	double    uu, uv, vv, wu, wv, D;
	uu = u.dot(u);
	uv = u.dot(v);
	vv = v.dot(v);
	w = I - T[0];
	wu = w.dot(u);
	wv = w.dot(v);
	D = uv * uv - uu * vv;

	// get and test parametric coords
	double s, t;
	s = (uv * wv - vv * wu) / D;
	if (s < 0.0 || s > 1.0)         // I is outside T
		return 0;
	t = (uv * wu - uu * wv) / D;
	if (t < 0.0 || (s + t) > 1.0)  // I is outside T
		return 0;

	return 1;                       // I is in T
}
inline void MyInteraction::getRay(unsigned mouseX, unsigned mouseY, std::pair<Vec3,Vec3> &Ray)
{
	GLint aViewport[4];
	glGetIntegerv(GL_VIEWPORT, aViewport);
	aViewport[2] = windowSize_x;
	aViewport[3] = windowSize_y;


	GLdouble matMV[16], matProj[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, matMV);
	GLdouble wx, wy, wz;  //  temp world x, y, z coords
	glGetDoublev(GL_PROJECTION_MATRIX, matProj);
	//  note viewport[3] is height of window in pixels
	mouseY = aViewport[3] - (GLint)mouseY - 1;
	//printf ("Coordinates at cursor are (%4d, %4d)\n", mouseX, mouseY);

	gluUnProject((GLdouble)mouseX, (GLdouble)mouseY, 0,
		matMV, matProj, aViewport, &wx, &wy, &wz);
	Ray.first = Vec3(wx, wy, wz);
	gluUnProject((GLdouble)mouseX, (GLdouble)mouseY, 1,
		matMV, matProj, aViewport, &wx, &wy, &wz);
	Ray.second = Vec3(wx, wy, wz);
}
inline bool MyInteraction::glSelection1(unsigned mouse_x, unsigned mouse_y, unsigned& hittedTriangle)											// This Is Where Selection Is Done
{
	if (mesh == NULL) return false;

	if (gl_buffer == NULL)		gl_buffer = new GLuint[mesh->getFaces().size()];
	GLuint	*buffer = gl_buffer;										// Set Up A Selection Buffer
	//buffer = new GLuint[faces.size()];
	GLint	hits;												// The Number Of Objects That We Selected

	// The Size Of The Viewport. [0] Is <x>, [1] Is <y>, [2] Is <length>, [3] Is <width>
	GLint	viewport[4];

	// This Sets The Array <viewport> To The Size And Location Of The Screen Relative To The Window
	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(mesh->getFaces().size(), buffer);								// Tell OpenGL To Use Our Array For Selection

	// Puts OpenGL In Selection Mode. Nothing Will Be Drawn.  Object ID's and Extents Are Stored In The Buffer.
	(void)glRenderMode(GL_SELECT);

	glInitNames();												// Initializes The Name Stack
	glPushName(0);												// Push 0 (At Least One Entry) Onto The Stack

	glMatrixMode(GL_PROJECTION);								// Selects The Projection Matrix
	glPushMatrix();												// Push The Projection Matrix
	glLoadIdentity();											// Resets The Matrix

	// This Creates A Matrix That Will Zoom Up To A Small Portion Of The Screen, Where The Mouse Is.
	gluPickMatrix((GLdouble)mouse_x, (GLdouble)(viewport[3] - mouse_y), 1.0f, 1.0f, viewport);

	// Apply The Perspective Matrix
	gluPerspective(45.0f, (GLfloat)(viewport[2] - viewport[0]) / (GLfloat)(viewport[3] - viewport[1]), 0.1f, 100.1f);
	glMatrixMode(GL_MODELVIEW);									// Select The Modelview Matrix
	drawTargets();												// Render The Targets To The Selection Buffer
	glMatrixMode(GL_PROJECTION);								// Select The Projection Matrix
	glPopMatrix();												// Pop The Projection Matrix
	glMatrixMode(GL_MODELVIEW);									// Select The Modelview Matrix
	hits = glRenderMode(GL_RENDER);								// Switch To Render Mode, Find Out How Many
	// Objects Were Drawn Where The Mouse Was
	if (hits > 0)												// If There Were More Than 0 Hits
	{
		int	choose = buffer[3];									// Make Our Selection The First Object
		int depth = buffer[1];									// Store How Far Away It Is 

		for (int loop = 1; loop < hits; loop++)					// Loop Through All The Detected Hits
		{
			// If This Object Is Closer To Us Than The One We Have Selected
			if (buffer[loop * 4 + 1] < GLuint(depth))
			{
				choose = buffer[loop * 4 + 3];						// Select The Closer Object
				depth = buffer[loop * 4 + 1];						// Store How Far Away It Is
			}
		}

		hittedTriangle = choose;
		return true;
	}
	else return false;
}
inline void MyInteraction::initSelectionColor()
{
	if (mesh == NULL)return;

	unsigned fsize = mesh->getFaces().size();
	glObjColors.resize(fsize);

	unsigned r, g, b;
	for (unsigned i = 0; i<fsize; i++)
	{
		AML::ubyte3 &col = glObjColors[i];

		r = i / 65536;
		g = (i - r * 65536) / 256;
		b = (i - r * 65536 - g * 256);

		col.x = r;
		col.y = g;
		col.z = b;
	}
	std::cout << "initial color map" << std::endl;
}
inline void MyInteraction::updateBackBuffer(int pos_x,int pos_y,int width,int height)
{
	if (mesh == NULL) return;

	glDrawBuffer(GL_BACK);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);
	for (auto f = mesh->getFaces().begin(); f != mesh->getFaces().end(); f++)
	{
		AML::ubyte3 &col = glObjColors[f->id()];
		glColor3ub(col.x, col.y, col.z);
		for (int j = 0; j<3; j++)
		{
			auto& v = f->vertex_iter(j)->coordinate();
			glVertex3d(v.x, v.y, v.z);
		}
	}
	glEnd();

	windowSize_x = width;
	windowSize_y = height;

	if (rgbBuffer != NULL) delete[]rgbBuffer;
	rgbBuffer = new unsigned char[width*height * 3];

	std::cout << "buffer updated" << std::endl;

	if (rgbBuffer != NULL)
	{
		glReadBuffer(GL_BACK);
		std::cout << "read back buffer..." << std::endl;
		glReadPixels((GLint)0, (GLint)pos_y, (GLint)width, (GLint)height, GL_RGB, GL_UNSIGNED_BYTE, rgbBuffer);
	}


	if (false)
	{
		byte*  bmpBuffer = (byte*)malloc(windowSize_x*windowSize_y * 3);
		char fileName[400];
		if (bmpBuffer != NULL) 
		{
			glReadBuffer(GL_BACK);
			glReadPixels((GLint)0, (GLint)pos_y, (GLint)windowSize_x, (GLint)windowSize_y, GL_RGB, GL_UNSIGNED_BYTE, bmpBuffer);
			for (int i = 0; i < windowSize_x*windowSize_y * 3; i += 3)
			{
				byte tmp = *(bmpBuffer + i);
				*(bmpBuffer + i) = *(bmpBuffer + i + 2);
				*(bmpBuffer + i + 2) = tmp;
			}

			QString imageFile("debug.off");
			imageFile.replace(imageFile.lastIndexOf(".") + 1, 3, "bmp");
			FILE *filePtr = fopen(imageFile.toStdString().data(), "wb");

			if (filePtr != NULL)
			{
				BITMAPFILEHEADER  bitmapFileHeader;
				BITMAPINFOHEADER  bitmapInfoHeader;
				bitmapFileHeader.bfType = ((WORD)('M' << 8) | 'B'); ;  //"BM"
				bitmapFileHeader.bfSize = windowSize_x*windowSize_y * 3 + sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
				bitmapFileHeader.bfReserved1 = 0;
				bitmapFileHeader.bfReserved2 = 0;
				bitmapFileHeader.bfOffBits =
					sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
				bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
				bitmapInfoHeader.biWidth = windowSize_x - 1;
				bitmapInfoHeader.biHeight = windowSize_y - 1;
				bitmapInfoHeader.biPlanes = 1;
				bitmapInfoHeader.biBitCount = 24;
				bitmapInfoHeader.biCompression = BI_RGB;
				bitmapInfoHeader.biSizeImage = 0;
				bitmapInfoHeader.biXPelsPerMeter = 0; // ?
				bitmapInfoHeader.biYPelsPerMeter = 0; //  ?
				bitmapInfoHeader.biClrUsed = 0;
				bitmapInfoHeader.biClrImportant = 0;

				fwrite(&bitmapFileHeader, sizeof(BITMAPFILEHEADER), 1, filePtr);
				fwrite(&bitmapInfoHeader, sizeof(BITMAPINFOHEADER), 1, filePtr);
				for (int i = 0; i < windowSize_x*windowSize_y * 3; i++)
					fputc(bmpBuffer[i], filePtr);
				//	fwrite(bmpBuffer, windowWidth*windowHeight*3, 1,  filePtr);
				fflush(filePtr);
				fclose(filePtr);
			}
			free(bmpBuffer);
		}

	}

}
inline bool MyInteraction::glSelection2(unsigned mouse_x, unsigned mouse_y, unsigned& hittedTriangle)
{
	if (rgbBuffer == NULL) return false;

	unsigned pixelID = (mouse_x + (windowSize_y - mouse_y)*windowSize_x) * 3;
	if (pixelID >= windowSize_y*windowSize_x * 3) return false;

	int r = int(*(rgbBuffer + pixelID + 0));
	int g = int(*(rgbBuffer + pixelID + 1));
	int b = int(*(rgbBuffer + pixelID + 2));

	unsigned hitId = b + g * 256 + r * 256 * 256;

	if (hitId > glObjColors.size() || hitId == 0)
		return false;
	else{
		hittedTriangle = hitId;
		return true;
	}
}
inline bool MyInteraction::glSelection2_Rectangle(unsigned mouse_x, unsigned mouse_y, unsigned mouse_x2, unsigned mouse_y2, std::set<unsigned>& hittedTriangles)
{
	if (rgbBuffer == NULL) return false;

	if (mouse_x > mouse_x2) std::swap(mouse_x, mouse_x2);
	if (mouse_y > mouse_y2) std::swap(mouse_y, mouse_y2);
	unsigned width = mouse_x2 - mouse_x;
	unsigned height = mouse_y2 - mouse_y;

// 	std::cout << "width" << width << " height:" << height << std::endl;
	hittedTriangles.clear();
	for (unsigned i = 0; i < width; i++)
	{
		for (unsigned j = 0; j < height; j++)
		{
			//viewport pixel index is not the same as which in back-buffer;
			//viewport: left-below, back-buffer: left-above
			unsigned pixelID = (mouse_x + i + (windowSize_y - mouse_y - j)*windowSize_x) * 3;
			if (pixelID >= windowSize_y*windowSize_x * 3) continue;

			unsigned r = unsigned(*(rgbBuffer + pixelID + 0));
			unsigned g = unsigned(*(rgbBuffer + pixelID + 1));
			unsigned b = unsigned(*(rgbBuffer + pixelID + 2));

			unsigned hitId = b + g * 256 + r * 256 * 256;

			if (hitId > glObjColors.size())	continue;

			hittedTriangles.insert(hitId);
		}
	}
	if (!hittedTriangles.empty())
		return true;
	else
		return false;
}
inline unsigned MyInteraction::getClosestVertex(unsigned mouse_x, unsigned mouse_y, unsigned hittedTriangle)
{
	std::pair<Vec3, Vec3> rays;
	getRay(mouse_x, mouse_y, rays);
	Vec3&rayStr = rays.first;
	Vec3&rayEnd = rays.second;

	double minDistance = DBL_MAX;

	float p0p1LenSquared = (rayEnd - rayStr).dot(rayEnd - rayStr);

	auto f = mesh->getFIter()[hittedTriangle];
	unsigned vid = 0;
	for (int i = 0; i < 3; i++)
	{
		Vec3 &p = f->vertex_iter(i)->coordinate();
		double paramU = (
			((p[0] - rayStr[0])*(rayEnd[0] - rayStr[0])) +
			((p[1] - rayStr[1])*(rayEnd[1] - rayStr[1])) +
			((p[2] - rayStr[2])*(rayEnd[2] - rayStr[2]))
			) / p0p1LenSquared;
		Vec3 lineP = rayStr + (rayEnd - rayStr)*paramU;
		double diatance = (lineP - p).length();
		if (minDistance > diatance)
		{
			minDistance = diatance;
			vid = f->vertex_iter(i)->id();
		}
	}
	return vid;
}
inline double MyInteraction::getClosestVertex(unsigned mouse_x, unsigned mouse_y, const std::vector<std::vector<Vec3>>& curves,std::pair<unsigned,unsigned>& s)
{
	std::pair<Vec3, Vec3> rays;
	getRay(mouse_x, mouse_y, rays);
	Vec3&rayStr = rays.first;
	Vec3&rayEnd = rays.second;

	double minDistance = FLT_MAX;
	unsigned arcID, ind;
	float p0p1LenSquared = (rayEnd - rayStr).dot(rayEnd - rayStr);

	for (int i = 0; i < curves.size(); i++)
	{
		for (int j = 0; j < curves[i].size(); j++)
		{
			Vec3 p = curves[i][j];
			double paramU = (
				((p[0] - rayStr[0])*(rayEnd[0] - rayStr[0])) +
				((p[1] - rayStr[1])*(rayEnd[1] - rayStr[1])) +
				((p[2] - rayStr[2])*(rayEnd[2] - rayStr[2]))
				) / p0p1LenSquared;

			Vec3 lineP = rayStr + (rayEnd - rayStr)*paramU;
			double diatance = (lineP - p).length();
			if (minDistance > diatance)
			{
				minDistance = diatance;
				arcID = i; ind = j;
			}
		}
	}

	s.first = arcID; s.second = ind;

	return minDistance;
}
inline void MyInteraction::genStraightScreenLine(unsigned str_x, unsigned str_y, unsigned end_x, unsigned end_y, std::vector<std::pair<unsigned, unsigned>>& sequence)
{
	double x = (double)str_x;
	double y = (double)str_y;

	int steps = abs((int)end_x - (int)str_x)>abs((int)end_y - (int)str_y) ? abs((int)end_x - (int)str_x) : abs((int)end_y - (int)str_y);

	double cx = (double)((int)end_x - (int)str_x) / (double)steps;
	double cy = (double)((int)end_y - (int)str_y) / (double)steps;

	for (int i = 0; i < steps; i++)
	{
		x = x + cx;
		y = y + cy;

		sequence.push_back(std::pair<unsigned, unsigned>((unsigned)(x + 0.5), (unsigned)(y + 0.5)));
	}
}
inline double MyInteraction::pointToPlane(Vec3& p, Vec3& pl_n, Vec3& pl_p)
{
	p -= pl_p;
	return pl_n.dot(p);
}
inline void MyInteraction::connectTwoTriangles(Vec3& p, Vec3& norm, Vec3& dir, unsigned f1, unsigned f2, std::vector<Vec3>& pts, std::vector<unsigned>& fs)
{
	const auto& fstr = mesh->getFIter()[f1];
	const auto& fend = mesh->getFIter()[f2];
	//check if two face are connected; yes, then return; no, then upsampling a triangle strip to connect them.
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (fstr->edge_iter(i)->id() == fend->edge_iter(j)->id())
				return;
		}
	}

	std::vector<double> signDist(3, -10.0);
	for (int i = 0; i < 3; i++)
	{
		auto f_it = fstr->edge_iter(i)->opposite_face(fstr);
		int count = 0;
		for (int j = 0; j < 3; j++)
		{
			auto q = f_it->vertex_iter(j)->coordinate();
			if (pointToPlane(q, norm, p)>0) count++;
		}

		if (count == 1 || count == 2)
		{
			auto q = f_it->center_point() - p;
			//cout<<" c:"<<count<<" dis:"<<q.dot(dir);
			signDist[i] = q.dot(dir);
			if (q.dot(dir) > 0)
			{
				if (std::find(fs.begin(), fs.end(), f_it->id()) == fs.end())
				{
					pts.push_back(f_it->center_point());
					fs.push_back(f_it->id());
					connectTwoTriangles(p, norm, dir, f_it->id(), f2, pts, fs);
				}
				//cout<<" f:"<<f_it->id();
				return;
			}
		}
	}
	int ind = std::max_element(signDist.begin(), signDist.end()) - signDist.begin();
	if (signDist[ind] > -10.0)
	{
		auto f_it = fstr->edge_iter(ind)->opposite_face(fstr);
		if (std::find(fs.begin(), fs.end(), f_it->id()) != fs.end())
		{
			ind = (ind + 1) % 3;
			if (signDist[ind] == -10.0)
				ind = (ind + 1) % 3;
		}
		pts.push_back(f_it->center_point());
		fs.push_back(f_it->id());
		//cout<<" f:"<<f_it->id();
		connectTwoTriangles(p, norm, dir, f_it->id(), f2, pts, fs);
	}
}
inline bool MyInteraction::scribleOnMesh(unsigned mouseX, unsigned mouseY, Hit& lastHit, Hit& newhit, bool isSingleHit, std::vector<Hit>& hits)
{
	//show the vertex under mouse; save to selectTriangle, hittedVertex
	//if scribbling, then draw the whole path;  save to scribbling, scribbleTriangleStrip
	bool hitted = glSelection2(mouseX, mouseY, newhit.hittedTriangle); 

	if (hitted)
	{
		//store intersection point;
		std::pair<Vec3, Vec3> Ray;
		getRay(mouseX, mouseY, Ray);

		std::vector<Vec3> Triangle(3);
		auto f = mesh->getFIter()[newhit.hittedTriangle];
		for (unsigned i = 0; i < 3; i++)
			Triangle[i] = f->vertex_iter(i)->coordinate();

		if (triangleLineIntersection(Ray, Triangle, newhit.hittedVertex) == 1)
		{
			if (!isSingleHit) //draw a path;
			{
				if (lastHit.hittedTriangle != newhit.hittedTriangle) //hits a new vertex, generates a strip of triangles connecting new and last vertex.
				{
// 					auto axix = newhit.hittedNearVertex - lastHit.hittedVertex;
// 					auto axiy = newhit.hittedVertex - lastHit.hittedVertex; axiy.normalize();
// 					auto norm = axix.cross(axiy); norm.normalize();
// 
// 					std::vector<Vec3> pts;
// 					std::vector<unsigned> fs;
// 					connectTwoTriangles(lastHit.hittedVertex, norm, axiy, lastHit.hittedTriangle, newhit.hittedTriangle, pts, fs);
// 
// 					hits.clear(); hits.resize(pts.size());
// 	
// 					for (unsigned i = 0; i < pts.size(); i++)
// 					{
// 						hits[i].hittedVertex = pts[i];
// 						hits[i].hittedTriangle = fs[i];
// 					}
					hits.push_back(newhit);
				}
			}
			//newhit.hittedNearVertex = newhit.hittedVertex; 
			//display hitted vertex near screen(close Ray.first) otherwise, may be obstructed by mesh.
			newhit.hittedNearVertex = 0.9*newhit.hittedVertex + 0.1*Ray.first;
			return true;
		}
		else
		{
// 			newhit.hittedNearVertex = 0.001*newhit.hittedVertex + 0.999*Ray.first;
			return false;
		}
	}
	else
	{
		return false;
	}


}
inline bool MyInteraction::brushOnMesh(unsigned mouseX, unsigned mouseY, unsigned nbsize, std::set<unsigned>& addToTriangles, Hit& hit)
{
	//grow a region from mouse location;
// 	std::set<std::pair<unsigned, unsigned> > mouseRegion;
// 	mouseRegion.insert(std::pair<unsigned, unsigned>(mouseX, mouseY));
// 	std::set<std::pair<unsigned, unsigned> > frontRegion = mouseRegion;
// 	while (nbsize > 0)
// 	{
// 		nbsize--;
// 		while (!frontRegion.empty())
// 		{
// 			auto p = *frontRegion.begin(); frontRegion.erase(frontRegion.begin());
// 
// 			//four neighboures;
// 			std::vector<std::pair<unsigned, unsigned>> ng(4, p);
// 			ng[0].first--; //left
// 			ng[1].first++; //right
// 			ng[2].second--;//up
// 			ng[3].second++;//down
// 			for (auto i = ng.begin(); i != ng.end(); i++)
// 			{
// 				mouseRegion.insert(*i);
// 			}
// 		}
// 
// 		frontRegion = mouseRegion;
// 	}

	int isize = (int)nbsize;
	int x = (int)mouseX;
	int y = (int)mouseY;
	std::vector<std::pair<int, int> > mouseRegion;
	for (int i = -isize; i < isize; i++)
	{
		for (int j = - (int)nbsize; j < (int)nbsize; j++)
		{
			if (abs(j) <= abs(isize - abs(i)))
			{
				mouseRegion.push_back(std::pair<int, int>(x + i, y + j));
			}
		}
	}

	std::set<unsigned> hittedTriangles;
	for (auto i = mouseRegion.begin(); i != mouseRegion.end(); i++)
	{
		unsigned hittedTriangle;
		bool hitted = glSelection2(i->first, i->second, hittedTriangle);
		if (hitted) hittedTriangles.insert(hittedTriangle);
	}

	for (auto i = hittedTriangles.begin(); i != hittedTriangles.end(); i++)
	{
		addToTriangles.insert(*i);
	}

	return scribleOnMesh(mouseX, mouseY, Hit(), hit); //update intersection position;
}
inline void MyInteraction::pushToMultipleScribble()
{
	if (!scribbleCurve.empty())
	{
		scribbleCurves.push_back(scribbleCurve);
		scribbleCurve.clear();
	}
	if (!scribbleTriangleStrip.empty())
	{
		scribbleTriangleStrips.push_back(scribbleTriangleStrip);
		scribbleTriangleStrip.clear();
	}
}

inline void MyInteraction::processSketchForSeg(std::vector<unsigned>& sketchFaces, std::vector<Vec3>& sketchCurves)
{
	if (scribbleTriangleStrips.empty()) return;
	sketchFaces.clear();
	sketchCurves.clear();
	if (scribbleType == Interaction::ADD)
	{
		for (unsigned i = 0; i < scribbleTriangleStrips.size(); i++)
		{
			if (!scribbleTriangleStrips[i].empty())
			{
				sketchFaces.push_back(scribbleTriangleStrips[i][0]);
				sketchCurves.push_back(scribbleCurves[i][0]);
			}
			for (unsigned j = 1; j < scribbleTriangleStrips[i].size(); j++)
			{
				if (scribbleTriangleStrips[i][j - 1] == scribbleTriangleStrips[i][j])
					continue;
				sketchFaces.push_back(scribbleTriangleStrips[i][j]);
				sketchCurves.push_back(scribbleCurves[i][j]);
			}
		}
	}
	else
	{
		for (unsigned i = 0; i < scribbleTriangleStrips.size(); i++)
		{
			for (unsigned j = 0; j < scribbleTriangleStrips[i].size(); j++)
			{
				sketchFaces.push_back(scribbleTriangleStrips[i][j]);
			}
		}
	}
}
inline void MyInteraction::userSketchSamplingFromSketchFaces(std::vector<unsigned>& sketchFaces, std::vector<Vec3>& sketchCurves, std::vector<SketchFace>& sketchElements, std::vector<unsigned>& eis)
{
	if (sketchFaces.empty() || sketchCurves.empty()) return;

	//the sketches are sequence of triangles, with each triangle a hitted vertex.
	//an usable feature edge lies in a triangle, with two ends at two triangle edges.
	//sampling from sketches before use it.

	//add disconnectted sketch;
	{
		std::vector<unsigned> newSketchFaces(1, sketchFaces.front());
		std::vector<Vec3> newSketchCurves(1, sketchCurves.front());
		std::vector<bool> usedFaces(mesh->getFaces().size(), false);
		usedFaces[sketchFaces[0]] = true;
		unsigned currentId = 0;
		for (unsigned i = 1; i < sketchFaces.size(); i++)
		{
			if (usedFaces[sketchFaces[i]]) continue;
			usedFaces[sketchFaces[i]] = true;

			auto face_str = mesh->getFIter()[sketchFaces[currentId]];
			auto face_end = mesh->getFIter()[sketchFaces[i]];
			//find common vertex;
			MyMesh::VertexIter v;
			unsigned c = 0;
			for (unsigned j = 0; j < 3; j++)
			{
				for (unsigned k = 0; k < 3; k++)
				{
					if (face_str->vertex_iter(j)->id() == face_end->vertex_iter(k)->id())
					{
						v = face_str->vertex_iter(j); c++;
					}
				}
			}

			if (c == 0)
			{
				std::cout << "isolated sketch ";
			}
			else if (c == 1)
			{
				//unfolding from f1 to f2; {f1,e1,v} -> {fnext,enext,v}... until fnext==f2;
				//find e1, and it's representative
				auto v2 = face_str->next_vertex(v);
				auto v3 = face_str->next_vertex(v2);

				std::vector<unsigned> face_clockwise;
				std::vector<unsigned> edge_clockwise;
				auto edge_next = face_str->opposite_edge(v2);
				edge_clockwise.push_back(edge_next->id());
				auto face_next = face_str;
				while (true)
				{
					face_next = edge_next->opposite_face(face_next);
					auto v1 = edge_next->vertex_iter(0);
					if (edge_next->vertex_iter(0)->id() == v->id())
						v1 = edge_next->vertex_iter(1);
					edge_next = face_next->opposite_edge(v1);

					if (face_next->id() == face_end->id())
						break;

					face_clockwise.push_back(face_next->id());
					edge_clockwise.push_back(edge_next->id());
				}
				std::vector<Vec3> edge_clockwise_pos;
				double errAccu_clockwise = 0;
				for (unsigned j = 0; j < edge_clockwise.size(); j++)
				{
					auto e = mesh->getEIter()[edge_clockwise[j]];
					auto vl = e->vertex_iter(0)->coordinate();
					auto vr = e->vertex_iter(1)->coordinate();

					Vec3 Il, Ir;
					double dl, dr;
					vertexLineIntersection(std::pair<Vec3, Vec3>(vl, vr), sketchCurves[currentId], Il, dl);
					vertexLineIntersection(std::pair<Vec3, Vec3>(vl, vr), sketchCurves[i], Ir, dr);

					double d = dl + dr;
					dl /= d; dr /= d;

					Vec3 I = Il*dr + Ir*dl;
					edge_clockwise_pos.push_back(I);

					//check if the Intersection is inside edge or not;
					if ((I - vl).dot(vr - I) < 0)std::cout << "bad intersection" << std::endl;
					errAccu_clockwise += (I - sketchCurves[currentId]).length() + (sketchCurves[i] - I).length() - (sketchCurves[currentId] - sketchCurves[i]).length();
				}

				std::vector<unsigned> face_anticlockwise;
				std::vector<unsigned> edge_anticlockwise;
				edge_next = face_str->opposite_edge(v3);
				edge_anticlockwise.push_back(edge_next->id());
				face_next = face_str;
				while (true)
				{
					face_next = edge_next->opposite_face(face_next);
					auto v1 = edge_next->vertex_iter(0);
					if (edge_next->vertex_iter(0)->id() == v->id())
						v1 = edge_next->vertex_iter(1);
					edge_next = face_next->opposite_edge(v1);

					if (face_next->id() == face_end->id())
						break;

					face_anticlockwise.push_back(face_next->id());
					edge_anticlockwise.push_back(edge_next->id());
				}
				std::vector<Vec3> edge_anticlockwise_pos;
				double errAccu_anticlockwise = 0;
				for (unsigned j = 0; j < edge_anticlockwise.size(); j++)
				{
					auto e = mesh->getEIter()[edge_anticlockwise[j]];
					auto vl = e->vertex_iter(0)->coordinate();
					auto vr = e->vertex_iter(1)->coordinate();

					Vec3 Il, Ir;
					double dl, dr;
					vertexLineIntersection(std::pair<Vec3, Vec3>(vl, vr), sketchCurves[currentId], Il, dl);
					vertexLineIntersection(std::pair<Vec3, Vec3>(vl, vr), sketchCurves[i], Ir, dr);

					double d = dl + dr;
					dl /= d; dr /= d;

					Vec3 I = Il*dr + Ir*dl;
					edge_anticlockwise_pos.push_back(I);

					//check if the Intersection is inside edge or not;
					if ((I - vl).dot(vr - I) < 0)std::cout << "bad intersection" << std::endl;
					errAccu_anticlockwise += (I - sketchCurves[currentId]).length() + (sketchCurves[i] - I).length() - (sketchCurves[currentId] - sketchCurves[i]).length();
				}


				if (errAccu_clockwise > errAccu_anticlockwise)
				{
					for (unsigned j = 0; j < face_anticlockwise.size(); j++)
					{
						newSketchFaces.push_back(face_anticlockwise[j]);
						newSketchCurves.push_back((edge_anticlockwise_pos[j] + edge_anticlockwise_pos[j + 1]) / 2.0);
						usedFaces[face_anticlockwise[j]] = true;
					}
				}
				else
				{
					for (unsigned j = 0; j < face_clockwise.size(); j++)
					{
						newSketchFaces.push_back(face_clockwise[j]);
						newSketchCurves.push_back((edge_clockwise_pos[j] + edge_clockwise_pos[j + 1]) / 2.0);
						usedFaces[face_clockwise[j]] = true;
					}
				}
				newSketchFaces.push_back(sketchFaces[i]);
				newSketchCurves.push_back(sketchCurves[i]);
				usedFaces[sketchFaces[i]] = true;
			}
			else
			{
				newSketchFaces.push_back(sketchFaces[i]);
				newSketchCurves.push_back(sketchCurves[i]);
			}
			currentId = i;
		}
		newSketchFaces.swap(sketchFaces);
		newSketchCurves.swap(sketchCurves);
	}

	//processing
	std::map < unsigned, std::pair<Vec3, double> > fstat;
	for (unsigned i = 0; i < sketchFaces.size(); i++)
	{
		if (fstat.find(sketchFaces[i]) == fstat.end())
		{
			fstat[sketchFaces[i]] = std::pair<Vec3, double>(sketchCurves[i], 1);
		}
		else
		{
			fstat[sketchFaces[i]] = std::pair<Vec3, double>(sketchCurves[i] + fstat[sketchFaces[i]].first, 1 + fstat[sketchFaces[i]].second);
		}
	}

	struct sketchEdgeElement
	{
		sketchEdgeElement() :count(0){}
		sketchEdgeElement(const Vec3& v1) :l(v1), count(1){}
		sketchEdgeElement(const Vec3& v1, const Vec3&v2) :l(v1), r(v2), count(2){}
		Vec3 l;
		Vec3 r;
		unsigned count;
	};
	std::map<unsigned, sketchEdgeElement> edgeShared;
	for (auto i = fstat.begin(); i != fstat.end(); i++)
	{
		auto tpos = i->second.first / i->second.second;
		auto f = mesh->getFIter()[i->first];
		for (unsigned j = 0; j < 3; j++)
		{
			unsigned eid = f->edge_iter(j)->id();
			if (edgeShared.find(eid) == edgeShared.end())
			{
				edgeShared[eid] = sketchEdgeElement(tpos);
			}
			else
			{
				if (edgeShared[eid].count < 2)
					edgeShared[eid] = sketchEdgeElement(edgeShared[eid].l, tpos);
			}
		}
	}

	eis.clear();
	std::map<unsigned, Vec3> edgeRepresent;
	for (auto i = edgeShared.begin(); i != edgeShared.end(); i++)
	{
		if (i->second.count == 2)
		{
			eis.push_back(i->first);

			//get the intersection; this method works only when the two lines are nearly co-planar; otherwise, it depends on the normal plane.
			auto e = mesh->getEIter()[i->first];
			auto vl = e->vertex_iter(0)->coordinate();
			auto vr = e->vertex_iter(1)->coordinate();

			Vec3 Il, Ir;
			double dl, dr;
			vertexLineIntersection(std::pair<Vec3, Vec3>(vl, vr), i->second.l, Il, dl);
			vertexLineIntersection(std::pair<Vec3, Vec3>(vl, vr), i->second.r, Ir, dr);

			double d = dl + dr;
			dl /= d; dr /= d;

			Vec3 I = Il*dr + Ir*dl;

			//check if the Intersection is inside edge or not;
			if ((I - vl).dot(vr - I) < 0)std::cout << "bad intersection" << std::endl;
			edgeRepresent[i->first] = I;
		}
	}

	std::map<unsigned, unsigned> fmap;
	sketchElements.clear();
	for (auto i = fstat.begin(); i != fstat.end(); i++)
	{
		SketchFace sk;
		sk.fid = i->first;
		sk.fpos = i->second.first / i->second.second;

		auto f = mesh->getFIter()[i->first];
		unsigned tid = 0;
		for (unsigned j = 0; j < 3; j++)
		{
			unsigned eid = f->edge_iter(j)->id();
			if (std::find(eis.begin(), eis.end(), eid) != eis.end())
			{
				if (tid == 0)
				{
					sk.e1 = eid;
					sk.ep1 = edgeRepresent[eid];
				}
				else
				{
					sk.e2 = eid;
					sk.ep2 = edgeRepresent[eid];
				}
				tid++;
			}
		}

		if (tid == 2)
		{
			sk.fpos = (sk.ep1 + sk.ep2) / 2.0;
			sketchElements.push_back(sk);
			fmap[sk.fid] = sketchElements.size() - 1;
		}
	}
}

#endif