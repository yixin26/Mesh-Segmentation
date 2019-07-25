#pragma once

#include <QtWidgets/QOpenGLWidget>
#include <Qt3DInput/QMouseEvent>

#include <time.h>
#include "myMesh.h"


class MyInteraction;
class MeshSegment;

class MyGLCanvas : public QOpenGLWidget
{
	Q_OBJECT

public:
	MyGLCanvas(QWidget* parent = 0, const char* name = 0, bool fs = false);
	~MyGLCanvas();

	MyInteraction** getInteractionUtils() { return &m_interaction; }
	MeshSegment** getMeshSegment() { return &m_meshSegment; }

	void paintGL();
	void resizeGL(int width, int height);
	void initializeGL();

	void mousePressEvent(QMouseEvent *event);  
	void mouseReleaseEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *e);

	void keyPressEvent(QKeyEvent *e);
	void keyReleaseEvent(QKeyEvent *e);

	int pos_y{ 0 };

private:

	MyInteraction * m_interaction;
	MeshSegment* m_meshSegment;

	bool m_initialized;

	int m_width;
	int m_height;

	//interaction configuration
	bool startDraw;
	bool isChangeView;

	bool isLeftDown{ false };
	bool isRightDown{ false };
	bool isMiddleDown{ false };
	int m_lastx;
	int m_lasty;

	int m_selectx;
	int m_selecty;
	std::vector<std::pair<int, int> > vts;


	//opengl parameter
	double m_RotationMatrix[16];
	Vec3 m_Eye;
	double m_PerspectiveAngleDegrees;
	double m_NearPlane;
	double m_FarPlaneOffset;

	void RevealFeature();
	void ConcealFeature();
	void processSketch(QMouseEvent *event);
	void RenderSelectionBuffer();
	void ScreenShot();
};

