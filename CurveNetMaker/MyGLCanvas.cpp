#include "MyGLCanvas.h"
#include <Qt3DInput/QKeyEvent>
#include <Qt3DInput/QMouseEvent>
#include <QtOpenGL/QGLWidget>

#include <GL/glu.h>
#include "core/Segmentation.h"

MyGLCanvas::MyGLCanvas(QWidget* parent, const char* name, bool fs) 
	:QOpenGLWidget(parent)
{
	m_initialized = false;
	startDraw = false;
	isChangeView = true;
}
MyGLCanvas::~MyGLCanvas()
{
}

void MyGLCanvas::initializeGL()
{
	setFocus();

	m_width = this->frameGeometry().width();
	m_height = this->frameGeometry().height();

	m_initialized = true;

	m_Eye[0] = 0.0f; m_Eye[1] = 0.0f; m_Eye[2] = -2.0f; //Actual code
	m_PerspectiveAngleDegrees = 45.0f;
	m_NearPlane = 0.01f;
	m_FarPlaneOffset = 100.0f;

	glViewport(0, 0, (GLint)m_width, (GLint)m_height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(-2, 2, -2, 2, m_NearPlane, m_NearPlane + m_FarPlaneOffset);

	gluPerspective(m_PerspectiveAngleDegrees,
		(GLfloat)m_width / (GLfloat)m_height,
		m_NearPlane,
		m_NearPlane + m_FarPlaneOffset);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//m_RotationMatrix.setIdentity(); //NO_MATRIX
	memset(m_RotationMatrix, 0.0, sizeof(double) * 16);
	for (unsigned int i = 0; i< 4; ++i) {
		m_RotationMatrix[i + (4 * i)] = 1.0;
	}

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading

	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glClearDepth(1.0f);									// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
														//glDepthFunc(GL_LESS);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);			// Really Nice Point Smoothing
	glEnable(GL_LINE_SMOOTH); //anti-aliased lines
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines


// **** rest of function is code from initGL function of Repousse standalone ****

	GLfloat light_position0[] = { -2.0f, -1.0f, 5.0f, 0.0f };
	GLfloat diffuse_light0[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat ambient_light0[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat specular_light0[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat As[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, As);

	glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_light0);


	GLfloat light_position1[] = { 2.0f, 1.0f, 0.5f, 0.0f };
	GLfloat diffuse_light1[] = { 0.2f,0.2f,0.2f,1.0f };
	GLfloat ambient_light1[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	GLfloat specular_light1[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse_light1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular_light1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient_light1);


	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	// enable color tracking
	bool isTwoSided = true;

	glEnable(GL_COLOR_MATERIAL);
	if (isTwoSided) {
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	}
	else {
		glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	}


	// ****** Material properties ******
	GLfloat mat_specular[] = { 0.6f, 0.6f, 0.6f, 1.0f };
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);

	GLfloat mat_shininess[] = { 8 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

	//GLfloat mat_diffuseB[] = { 0.6f, 0.6f, 0.8f, 1.0f };
	GLfloat mat_diffuseB[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuseB);

	GLfloat mat_ambient_front[] = { 0.1f, 0.1f, .1f, 1.0f };
	//GLfloat mat_ambient_front[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient_front);

	GLfloat mat_ambient_back[] = { 0.1f, 0.1f, .1f, 1.0f };//{ 1.0f, 1.0f, 1.0f, 1.0f };
	glMaterialfv(GL_BACK, GL_AMBIENT, mat_ambient_back);

	/*
	GLfloat mat_emit[] = { 0.1f, 0.1f, .1f, 1.0f };
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION , mat_emit);
	*/


	if (isTwoSided) {
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	}
	else {
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	}


	//glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);
	glEnable(GL_LIGHTING);

	double bsize = m_interaction->brushSize;
	GLint	viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	m_interaction->scribleRadius = bsize*1.2 / double(viewport[3]);

}
void MyGLCanvas::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(m_Eye[0], m_Eye[1], m_Eye[2]);

	glMultMatrixd(m_RotationMatrix);

	//draw surface
	m_FarPlaneOffset = 100.0f;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(m_PerspectiveAngleDegrees,
		(GLfloat)m_width / (GLfloat)m_height,
		m_NearPlane,
		m_NearPlane + m_FarPlaneOffset);

	glMatrixMode(GL_MODELVIEW);
	
	m_meshSegment->DrawGraph();

	//draw lines
	m_FarPlaneOffset = 200.0f;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(m_PerspectiveAngleDegrees,
		(GLfloat)m_width / (GLfloat)m_height,
		m_NearPlane,
		m_NearPlane + m_FarPlaneOffset);

	glMatrixMode(GL_MODELVIEW);
	m_meshSegment->DrawCurve();

	update();
}
void MyGLCanvas::resizeGL(int width, int height)
{
	if (!m_initialized)
	{
		initializeGL();
		m_initialized = true;
	}

	m_width = this->frameGeometry().width();
	m_height = this->frameGeometry().height();

	glViewport(0, 0, (GLint)m_width, (GLint)m_height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(m_PerspectiveAngleDegrees,
		(GLfloat)(m_width) / (GLfloat)m_height,
		m_NearPlane,
		m_NearPlane + m_FarPlaneOffset);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glMultMatrixd(m_RotationMatrix);

}

void MyGLCanvas::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_Shift:
		m_interaction->isShiftPress = true;
		break;
	case Qt::Key_Control:
		m_interaction->isControlPress = true; 
		break;
	case Qt::Key_Alt:
		m_interaction->isAltPress = true;
		break;
	case Qt::Key_Escape:
		exit(0);
		break;
	}
}
void MyGLCanvas::keyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_Shift:
		m_interaction->isShiftPress = false;
		break;
	case Qt::Key_Control:
		m_interaction->isControlPress = false;
		break;
	case Qt::Key_Alt:
		m_interaction->isAltPress = false;
		break;
	case Qt::Key_S:
		ScreenShot();
		break;
	}

	if (!m_interaction->isAltPress && !m_interaction->isShiftPress && !m_interaction->isControlPress)
	{
		m_interaction->isScrolling = false;
		startDraw = false;

		m_interaction->pushToMultipleScribble();
		m_interaction->processSketchForSeg(m_meshSegment->sketchFaces, m_meshSegment->sketchCurves);
		m_meshSegment->modifySegmentBySketch(m_interaction->scribbleType);

		m_interaction->clearSketches();
		m_meshSegment->realTimeSelectedTriangles.clear();
		m_meshSegment->sketchFaces.clear();
	}

}

void MyGLCanvas::mousePressEvent(QMouseEvent *event)
{
	switch (event->button())
	{
	case Qt::LeftButton:isLeftDown = true; break;
	case Qt::RightButton:isRightDown = true; break;
	case Qt::MiddleButton:isMiddleDown = true; break;
	default:break;
	}	
	m_lastx = event->pos().x();
	m_lasty = event->pos().y();

	paintGL();

}
void MyGLCanvas::mouseReleaseEvent(QMouseEvent *event) 
{
	//isLeftDown = isRightDown = isMiddleDown = false;
	switch (event->button())
	{
	case Qt::LeftButton:isLeftDown = false; break;
	case Qt::RightButton:isRightDown = false; break;
	case Qt::MiddleButton:isMiddleDown = false; break;
	default:break;
	}

	if (m_interaction->isAltPress) //split
	{
		if (event->button() == Qt::LeftButton)//interval
		{
			if (m_interaction->sr == Interaction::PIXEL)
			{
				m_meshSegment->smoothScribble(m_interaction->scribbleTriangleStrip, m_interaction->scribbleCurve);
				m_interaction->pushToMultipleScribble();
			}
			startDraw = false;
		}
		else if (event->button() == Qt::RightButton)
		{
			m_interaction->clearSketches();
			m_meshSegment->realTimeSelectedTriangles.clear();
		}
	}
	else if (m_interaction->isControlPress)
	{
		if (event->button() == Qt::LeftButton)//interval
		{
			if (m_interaction->sr == Interaction::PIXEL)
			{
				m_interaction->pushToMultipleScribble();
			}
			else if (m_interaction->sr == Interaction::RECTANGLE)
			{

			}
			startDraw = false;
		}
		else if (event->button() == Qt::RightButton)
		{
			m_interaction->clearSketches();
			m_meshSegment->realTimeSelectedTriangles.clear();
		}
	}
	else if (m_interaction->isShiftPress)
	{
		if (event->button() == Qt::LeftButton)//interval
		{
			if (m_interaction->sr == Interaction::PIXEL)
			{
				m_interaction->pushToMultipleScribble();
			}
			startDraw = false;
		}
		else if (event->button() == Qt::RightButton)
		{
			m_interaction->isScrolling = false;
			m_interaction->clearSketches();
			m_meshSegment->realTimeSelectedTriangles.clear();
		}
	}

	if (!m_interaction->isAltPress && !m_interaction->isShiftPress && !m_interaction->isControlPress)
	{
		m_interaction->isScrolling = false;
		startDraw = false;

		m_interaction->pushToMultipleScribble();
		m_interaction->processSketchForSeg(m_meshSegment->sketchFaces, m_meshSegment->sketchCurves);
		m_meshSegment->modifySegmentBySketch(m_interaction->scribbleType);

		m_interaction->clearSketches();
		m_meshSegment->realTimeSelectedTriangles.clear();
		m_meshSegment->sketchFaces.clear();
	}

	paintGL();
}

void MyGLCanvas::mouseMoveEvent(QMouseEvent *event)
{
	int x = event->pos().x();
	int y = event->pos().y();
	//if (x < 0 || y < 0 || x>m_width || y>m_height) return;
	
	if (!m_interaction->isAltPress && !m_interaction->isShiftPress && !m_interaction->isControlPress)
	{
		if (isLeftDown)
		{
			//convert the mouse clicked locations to lie between [-1,1] for both X and Y
			double halfWidth = m_width / 2.0;
			double halfHeight = m_height / 2.0;
			double xNormalized = (x - halfWidth) / halfWidth;
			double yNormalized = (halfHeight - y) / halfHeight;
			double oldXNormalized = (m_lastx - halfWidth) / halfWidth;
			double oldYNormalized = (halfHeight - m_lasty) / halfHeight;

			// rotates screen
			float rot[3] = { 0 };
			//rot[0] -= (m_lasty - y) * 0.5;
			//rot[1] -= (m_lastx - x) * 0.5;
			rot[1] -= (m_lasty - y) * 0.5;
			rot[0] -= (m_lastx - x) * 0.5;
			//------------------------------------------------------------------------
			// If rotation angle is greater of 360 or lesser than -360,
			// resets it back to zero.
			//------------------------------------------------------------------------
			for (unsigned int i = 0; i < 3; i++)
				if (rot[i] > 360 || rot[i] < -360)
					rot[i] = 0;

			//transfer this rotation into the rotation matrix for the scene
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glLoadIdentity();
			//glRotated(-rot [0], 0,1.0,0);
			//glRotated(-rot [1], 1.0,0,0);
			glRotated(rot[0], 0, 1.0, 0);
			glRotated(rot[1], 1.0, 0, 0);
			glMultMatrixd(m_RotationMatrix); //accummulate the previous transformations
			glGetDoublev(GL_MODELVIEW_MATRIX, m_RotationMatrix); //update
			glPopMatrix();
		}
		else if (isRightDown)
		{
			// translate the screen, z axis gives the idea of zooming in and out
			m_Eye[2] -= (m_lasty - y) * 0.05; // here I multiply by a 0.05 factor to slow down the zoom
		}
		else if (isMiddleDown)
		{
			m_Eye[0] -= (m_lastx - x) * 0.01f;
			m_Eye[1] += (m_lasty - y) * 0.01f;
		}
	
		if (isLeftDown || isMiddleDown || isRightDown)
		{
			paintGL();
			isChangeView = true; //change of view pos will cause re-render the back buffer for opengl selection;
		}
	}

	processSketch(event);

	m_lastx = x;
	m_lasty = y;
}

void MyGLCanvas::wheelEvent(QWheelEvent *e)
{
	if (!m_interaction->isAltPress && !m_interaction->isShiftPress && !m_interaction->isControlPress)
	{
		if (e->delta() > 0)
		{
			m_Eye[2] -= 0.01; // here I multiply by a 0.05 factor to slow down the zoom
		}
		else
		{
			m_Eye[2] += 0.05; // here I multiply by a 0.05 factor to slow down the zoom
		}
	}
	else
	{
		if (m_interaction->isControlPress && m_interaction->isShiftPress)
		{
			if (e->delta() > 0)
			{
				m_interaction->brushSize++;
			}
			else
			{
				m_interaction->brushSize--;
				if (m_interaction->brushSize < 1)
					m_interaction->brushSize = 1;
			}

			double bsize = m_interaction->brushSize;
			GLint	viewport[4];
			glGetIntegerv(GL_VIEWPORT, viewport);
			m_interaction->scribleRadius = bsize*1.2 / double(viewport[3]);
		}
		else if (m_interaction->isShiftPress)//&& isRightDown
		{
			if (e->delta() > 0)
			{
				m_interaction->scribbleType = Interaction::REVEAL;
			}
			else
			{
				m_interaction->scribbleType = Interaction::CONCEAL;
			}
			m_interaction->isScrolling = true;

			m_interaction->pushToMultipleScribble();
			m_interaction->processSketchForSeg(m_meshSegment->sketchFaces, m_meshSegment->sketchCurves);

			if (m_meshSegment->sketchFaces.empty())
			{
				if (m_interaction->scribbleType == Interaction::REVEAL)
					RevealFeature();
				else if (m_interaction->scribbleType == Interaction::CONCEAL)
					ConcealFeature();
			}
			else
			{
				m_meshSegment->modifySegmentBySketch(m_interaction->scribbleType);
			}
		}
	}
	paintGL();
}

void MyGLCanvas::RevealFeature()
{
	m_meshSegment->updateGlabalAwardAlpha();
	m_meshSegment->mergePartition();
	paintGL();
}
void MyGLCanvas::ConcealFeature()
{
	m_meshSegment->updateGlabalAwardAlpha();
	m_meshSegment->mergePartition();
	paintGL();
}

void MyGLCanvas::processSketch(QMouseEvent *event)//mouse move
{
	int x = event->pos().x();
	int y = event->pos().y();

	if (m_interaction->isAltPress) //split
	{
		m_interaction->isAltPress = true;
		m_interaction->scribbleType = Interaction::ADD;

		if (isChangeView)
		{
			RenderSelectionBuffer();
			isChangeView = false;
		}

		//m_selectx,m_selecty is last position; x,y is current position; only process when mouse moves;
		if (isLeftDown && m_selectx != x&&m_selecty != y)//on-drawing
		{
			if (m_interaction->sr == Interaction::PIXEL)
			{
				std::vector<std::pair<unsigned, unsigned> > sequence;
				if (startDraw)
				{
					m_interaction->genStraightScreenLine(m_selectx, m_selecty, x, y, sequence);

					for (int i = 0; i < sequence.size(); i++)
					{
						MyInteraction::Hit lastHit;
						bool isSingleHit = true;
						std::vector<MyInteraction::Hit> hits;
						if (!m_interaction->scribbleCurve.empty())
						{
							lastHit.hittedVertex = m_interaction->scribbleCurve.back();
							lastHit.hittedTriangle = m_interaction->scribbleTriangleStrip.back();
							isSingleHit = false;
						}
						m_interaction->scribleOnMesh(sequence[i].first, sequence[i].second, lastHit, m_interaction->hit, isSingleHit, hits);
						for (auto i = hits.begin(); i < hits.end(); i++)
						{
							if (std::find(m_interaction->scribbleTriangleStrip.begin(), m_interaction->scribbleTriangleStrip.end(), i->hittedTriangle) != m_interaction->scribbleTriangleStrip.end())continue;

							m_interaction->scribbleCurve.push_back(i->hittedVertex);
							m_interaction->scribbleTriangleStrip.push_back(i->hittedTriangle);
						}
						if (isSingleHit)
						{
							m_interaction->scribbleCurve.push_back(m_interaction->hit.hittedVertex);
							m_interaction->scribbleTriangleStrip.push_back(m_interaction->hit.hittedTriangle);
						}
					}
				}
				else
				{
					startDraw = true;
				}

				m_selectx = x;
				m_selecty = y;
			}
		}
		else//navigate
		{
			bool hitted = m_interaction->scribleOnMesh(x, y, MyInteraction::Hit(), m_interaction->hit);//not save the path;
		}
	}
	else if (m_interaction->isControlPress) //merge
	{
		m_interaction->isControlPress = true;
		m_interaction->scribbleType = Interaction::ERASE;
		if (isChangeView)
		{
			RenderSelectionBuffer();
			isChangeView = false;
		}

		//m_selectx,m_selecty is last position; x,y is current position; only process when mouse moves;
		if (isLeftDown && m_selectx != x&&m_selecty != y)//on-drawing
		{
			if (m_interaction->sr == Interaction::PIXEL)
			{
				std::vector<std::pair<unsigned, unsigned> > sequence;
				if (startDraw)
				{
					m_interaction->genStraightScreenLine(m_selectx, m_selecty, x, y, sequence);

					for (int i = 0; i < sequence.size(); i++)
					{
						m_interaction->brushOnMesh(sequence[i].first, sequence[i].second, m_interaction->brushSize, m_meshSegment->realTimeSelectedTriangles, m_interaction->hit);
					}
				}
				else
				{
					startDraw = true;
				}

				m_selectx = x;
				m_selecty = y;
			}

			m_interaction->scribbleTriangleStrip.clear();
			m_interaction->scribbleTriangleStrip.reserve(m_meshSegment->realTimeSelectedTriangles.size());
			for (auto i = m_meshSegment->realTimeSelectedTriangles.begin(); i != m_meshSegment->realTimeSelectedTriangles.end(); i++)
			{
				m_interaction->scribbleTriangleStrip.push_back(*i);
			}
		}
		else//navigate
		{
			bool hitted = m_interaction->scribleOnMesh(x, y, MyInteraction::Hit(), m_interaction->hit);//not save the path;
		}
	}
	else if (m_interaction->isShiftPress) //boost/decrease
	{
		m_interaction->isShiftPress = true;
		if (m_interaction->scribbleType != Interaction::REVEAL && m_interaction->scribbleType != Interaction::CONCEAL)
			m_interaction->scribbleType = Interaction::REVEAL;
		if (isChangeView)
		{
			RenderSelectionBuffer();
			isChangeView = false;
		}

		//m_selectx,m_selecty is last position; x,y is current position; only process when mouse moves;
		if (isLeftDown && m_selectx != x&&m_selecty != y)//on-drawing
		{
			if (m_interaction->sr == Interaction::PIXEL)
			{
				std::vector<std::pair<unsigned, unsigned> > sequence;
				if (startDraw)
				{
					m_interaction->genStraightScreenLine(m_selectx, m_selecty, x, y, sequence);

					for (int i = 0; i < sequence.size(); i++)
					{
						m_interaction->brushOnMesh(sequence[i].first, sequence[i].second, m_interaction->brushSize, m_meshSegment->realTimeSelectedTriangles, m_interaction->hit);
					}
				}
				else
				{
					startDraw = true;
				}

				m_selectx = x;
				m_selecty = y;
			}

			m_interaction->scribbleTriangleStrip.clear();
			m_interaction->scribbleTriangleStrip.reserve(m_meshSegment->realTimeSelectedTriangles.size());
			for (auto i = m_meshSegment->realTimeSelectedTriangles.begin(); i != m_meshSegment->realTimeSelectedTriangles.end(); i++)
			{
				m_interaction->scribbleTriangleStrip.push_back(*i);
			}
		}
		else//navigate
		{
			bool hitted = m_interaction->scribleOnMesh(x, y, MyInteraction::Hit(), m_interaction->hit);//not save the path;
		}
	}

	if (!m_interaction->isAltPress && !m_interaction->isShiftPress && !m_interaction->isControlPress)
	{
		m_interaction->isScrolling = false;
		startDraw = false;

		m_interaction->pushToMultipleScribble();
		m_interaction->processSketchForSeg(m_meshSegment->sketchFaces, m_meshSegment->sketchCurves);
		m_meshSegment->modifySegmentBySketch(m_interaction->scribbleType);

		m_interaction->clearSketches();
		m_meshSegment->realTimeSelectedTriangles.clear();
		m_meshSegment->sketchFaces.clear();
	}

	paintGL();
}
void MyGLCanvas::RenderSelectionBuffer()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_LIGHTING);
	//glDisable(GL_TEXTURE_2D);
	glDisable(GL_COLOR_MATERIAL);
	glShadeModel(GL_FLAT);
	glDisable(GL_POINT_SMOOTH);
	glDisable(GL_LINE_SMOOTH); //anti-aliased lines
	glDisable(GL_BLEND);
	
	//draw back buffer for opengl selection;
	glViewport(0, 0, (GLint)m_width, (GLint)m_height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(m_Eye[0], m_Eye[1], m_Eye[2]);

	glMultMatrixd(m_RotationMatrix);

	//draw surface
	m_FarPlaneOffset = 100.0f;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(m_PerspectiveAngleDegrees,
		(GLfloat)m_width / (GLfloat)m_height,
		m_NearPlane,
		m_NearPlane + m_FarPlaneOffset);

	glMatrixMode(GL_MODELVIEW);

	m_interaction->updateBackBuffer(0,0, this->frameGeometry().width(),this->frameGeometry().height());

	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	//glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH); //anti-aliased lines
	glEnable(GL_BLEND);
}

void MyGLCanvas::ScreenShot()
{
	setFocus();

	//GLint	viewport[4];
	//glGetIntegerv(GL_VIEWPORT, viewport);
	//windowWidth = viewport[2]; windowHeight = viewport[3];
	//windowWidth -= 10;
	int windowWidth, windowHeight;
	windowWidth = this->frameGeometry().width();
	windowHeight = this->frameGeometry().height();

	byte*  bmpBuffer = (byte*)malloc(windowWidth*windowHeight * 3);
	char fileName[400];
	if (bmpBuffer != NULL) {
		glReadBuffer(GL_BACK);
		glReadPixels((GLint)0, (GLint)pos_y, (GLint)windowWidth, (GLint)windowHeight, GL_RGB, GL_UNSIGNED_BYTE, bmpBuffer);
		for (int i = 0; i < windowWidth*windowHeight * 3; i += 3) 
		{
			byte tmp = *(bmpBuffer + i);
			*(bmpBuffer + i) = *(bmpBuffer + i + 2);
			*(bmpBuffer + i + 2) = tmp;
		}

		QString imageFile = m_meshSegment->m_filename;
		imageFile.replace(imageFile.lastIndexOf(".") + 1, 3, "bmp");
		FILE *filePtr = fopen(imageFile.toStdString().data(), "wb");
		cout << "save screen to " << imageFile.toStdString().data() << endl;

		if (filePtr != NULL) {
			BITMAPFILEHEADER  bitmapFileHeader;
			BITMAPINFOHEADER  bitmapInfoHeader;
			bitmapFileHeader.bfType = ((WORD)('M' << 8) | 'B'); ;  //"BM"
			bitmapFileHeader.bfSize = windowWidth*windowHeight * 3 + sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
			bitmapFileHeader.bfReserved1 = 0;
			bitmapFileHeader.bfReserved2 = 0;
			bitmapFileHeader.bfOffBits =
				sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
			bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
			bitmapInfoHeader.biWidth = windowWidth - 1;
			bitmapInfoHeader.biHeight = windowHeight - 1;
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
			for (int i = 0; i<windowWidth*windowHeight * 3; i++)
				fputc(bmpBuffer[i], filePtr);
			//	fwrite(bmpBuffer, windowWidth*windowHeight*3, 1,  filePtr);
			fflush(filePtr);
			fclose(filePtr);
		}
		free(bmpBuffer);
	}
}
