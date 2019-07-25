/***************************************************************
* Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
**************************************************************/

#include "core/Segmentation.h"
#include <algorithm>

void FeatureLine::drawCrestLine()
{
	if (unusedCrestline.empty())
		unusedCrestline.resize(crestEdges.size(), false);

	glLineWidth(featureWidth);
	glBegin(GL_LINES);
	for (unsigned i = 0; i < crestLineoffset; i++)
	{
		if (crestEdgesVisible[i] == false)continue;

		if (unusedCrestline[i] == true)
		{
// 			glColor3d(0.3, 0.3, 0.3);
			glColor3d(0.8, 0.8, 0.0);
		}
		else
		{
// 			glColor3d(0.3, 0.3, 0.3);
			glColor3d(1, 0, 0);
		}
		glVertex3d(crestPoints[crestEdges[i][0]].x, crestPoints[crestEdges[i][0]].y, crestPoints[crestEdges[i][0]].z);
		glVertex3d(crestPoints[crestEdges[i][1]].x, crestPoints[crestEdges[i][1]].y, crestPoints[crestEdges[i][1]].z);
	}
	for (unsigned i = crestLineoffset; i < crestEdges.size(); i++)
	{
		if (crestEdgesVisible[i] == false)continue;

		if (unusedCrestline[i] == true)
		{
// 			glColor3d(0.3, 0.3, 0.3);
			glColor3d(0.8, 0.8, 0.0);
		}
		else
		{
// 			glColor3d(0.3, 0.3, 0.3);
			glColor3d(0, 0, 1);
		}
		glVertex3d(crestPoints[crestEdges[i][0]].x, crestPoints[crestEdges[i][0]].y, crestPoints[crestEdges[i][0]].z);
		glVertex3d(crestPoints[crestEdges[i][1]].x, crestPoints[crestEdges[i][1]].y, crestPoints[crestEdges[i][1]].z);
	}
	glEnd();

	glLineWidth(1.);
}
void MeshSegment::drawUserSketch()
{
	glLineWidth(featureWidth);
	glBegin(GL_LINES);
	glColor3d(0.15, 0.7, 0.3);
	/*
	for(int i=0;i<sketchEdges.size();i++){
	glVertex3d(sketchPoints[sketchEdges[i][0]].x,sketchPoints[sketchEdges[i][0]].y,sketchPoints[sketchEdges[i][0]].z);
	glVertex3d(sketchPoints[sketchEdges[i][1]].x,sketchPoints[sketchEdges[i][1]].y,sketchPoints[sketchEdges[i][1]].z);
	}
	*/
	for (int i = 0; i<userSketches.size(); i++)
	{
		glVertex3d(userSketches[i].ep1.x, userSketches[i].ep1.y, userSketches[i].ep1.z);
		glVertex3d(userSketches[i].ep2.x, userSketches[i].ep2.y, userSketches[i].ep2.z);
	}

	glEnd();
}

void MeshSegment::drawGroupFeatures()
{
	if (vertexLabel.empty()) return;
	glLineWidth(1.0);
	glColor3d(0.5, 0.5, 0.5);
	glBegin(GL_LINES);
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		if (graphFeature.ef[e_it->id()].lab == -1)continue;
		if (graphFeature.ef[e_it->id()].isCut)continue;

		if (e_it->manifold())
		{
			unsigned fs[] = { e_it->face_iter(0)->id(), e_it->face_iter(1)->id() };
			glVertex3d(graphFeature.cf[fs[0]].rep.x, graphFeature.cf[fs[0]].rep.y, graphFeature.cf[fs[0]].rep.z);
			glVertex3d(graphFeature.ef[e_it->id()].rep.x, graphFeature.ef[e_it->id()].rep.y, graphFeature.ef[e_it->id()].rep.z);
			glVertex3d(graphFeature.ef[e_it->id()].rep.x, graphFeature.ef[e_it->id()].rep.y, graphFeature.ef[e_it->id()].rep.z);
			glVertex3d(graphFeature.cf[fs[1]].rep.x, graphFeature.cf[fs[1]].rep.y, graphFeature.cf[fs[1]].rep.z);
		}
		else
		{
			unsigned fs[] = { e_it->face_iter(0)->id() };
			glVertex3d(graphFeature.cf[fs[0]].rep.x, graphFeature.cf[fs[0]].rep.y, graphFeature.cf[fs[0]].rep.z);
			glVertex3d(graphFeature.ef[e_it->id()].rep.x, graphFeature.ef[e_it->id()].rep.y, graphFeature.ef[e_it->id()].rep.z);
		}
	}
	glEnd();
	glLineWidth(1.0);
}

void MeshSegment::genSegPatchColor(bool val)
{
	genPatchColor(patchColors, segNumber, val);
}
void MeshSegment::patchColorMinChanged(std::vector<unsigned>& tlabels, int patNum)
{
	std::vector<std::set < int > > multi_labelMap(segNumber);
	for (int i = 0; i < vertexLabel.size(); i++)
	{
		multi_labelMap[vertexLabel[i]].insert(tlabels[i]);
	}
	std::map<int, int> labelMap;
	std::vector<bool> labelMapped(segNumber, false);
	std::vector<bool> labelUsed(max(patNum, segNumber), false);
	for (int i = 0; i < multi_labelMap.size(); i++)
	{
		int tid = *multi_labelMap[i].begin();
		if (labelUsed[tid] == false)
		{
			labelMap[i] = tid;
			labelUsed[tid] = true;
			labelMapped[i] = true;
		}
	}
	for (int i = 0; i < labelMapped.size(); i++)
	{
		if (!labelMapped[i])
		{
			int newLabel = 0;
			while (labelUsed[newLabel])
			{
				newLabel++;
			}
			labelUsed[newLabel] = true;
			labelMapped[i] = true;
			labelMap[i] = newLabel;
		}
	}
	if (segNumber>patchColors.size())
	{
		genPatchColor(patchColors, max(segNumber, patNum), false);
	}
	std::vector<ColorEngine::color> tColors(segNumber);
	for (int i = 0; i < labelMap.size(); i++)
	{
		tColors[i] = patchColors[labelMap[i]];
	}
	patchColors = tColors;
}
void MeshSegment::genPatchColor(std::vector<ColorEngine::color>& patchColors, int segNumber, bool reGenerate)
{
	if (reGenerate)
	{
		patchColors.clear();
		unsigned renum = segNumber;
		while (renum > 0)
		{
			renum--;

			ColorEngine::color newcolor;
			newcolor.r = 1.0 - double(rand() % 1000) / 2000.0; //range from 0.5-1.0
			newcolor.g = 1.0 - double(rand() % 1000) / 2000.0;
			newcolor.b = 1.0 - double(rand() % 1000) / 2000.0;

			patchColors.push_back(newcolor);
		}
	}
	else
	{
		if (patchColors.size() > segNumber)
			return;
		else
		{
			unsigned renum = segNumber - patchColors.size();
			while (renum>0)
			{
				renum--;

				ColorEngine::color newcolor;
				newcolor.r = 1.0 - double(rand() % 1000) / 2000.0; //range from 0.5-1.0
				newcolor.g = 1.0 - double(rand() % 1000) / 2000.0;
				newcolor.b = 1.0 - double(rand() % 1000) / 2000.0;

				patchColors.push_back(newcolor);
			}
		}
	}
}
void MeshSegment::drawSegmentation()
{
	bool isSegReady = false;
	if (!vertexLabel.empty())
	{
		if (patchColors.size() != segNumber)
		{
			genSegPatchColor(false);
		}
		isSegReady = true;
	}

	GLfloat specular_light0[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light0);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular_light0);
	GLfloat As[4] = { 0.05f, 0.05f, 0.05f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, As);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);

	glColor3d(.5, 0.539, 0.527);

	for (auto f_it = myMesh->getFaces().begin(); f_it != myMesh->getFaces().end(); f_it++)
	{
		MyMesh::VertexIter v[] = { f_it->vertex_iter(0), f_it->vertex_iter(1), f_it->vertex_iter(2) };
		for (unsigned j = 0; j < 3; j++)
		{
			glNormal3f(v[j]->normal().x, v[j]->normal().y, v[j]->normal().z);
			if (isSegReady)
			{
				ColorEngine::color &tcolor = patchColors[vertexLabel[v[j]->id()]];
				glColor3d(tcolor.r, tcolor.b, tcolor.g);
			}
			glVertex3f(v[j]->coordinate().x, v[j]->coordinate().y, v[j]->coordinate().z);
		}
	}

	glEnd();

	GLfloat specular_light[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular_light);
	GLfloat As2[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, As2);
}

void MeshSegment::drawSegBoundary()
{
	if (vertexLabel.empty()) return;
	glLineWidth(boundaryWidth);
	if (boundaryCurves.empty())
	{
		std::vector<MyMesh::EdgeIter> selId;
		glBegin(GL_LINES);
		for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
		{
			if (graphFeature.ef[e_it->id()].isCut == false)continue;
			glColor3d(0, 0, 0);

			if (e_it->manifold())
			{
				unsigned fs[] = { e_it->face_iter(0)->id(), e_it->face_iter(1)->id() };
				glVertex3d(graphFeature.cf[fs[0]].rep.x, graphFeature.cf[fs[0]].rep.y, graphFeature.cf[fs[0]].rep.z);
				glVertex3d(graphFeature.ef[e_it->id()].rep.x, graphFeature.ef[e_it->id()].rep.y, graphFeature.ef[e_it->id()].rep.z);
				glVertex3d(graphFeature.ef[e_it->id()].rep.x, graphFeature.ef[e_it->id()].rep.y, graphFeature.ef[e_it->id()].rep.z);
				glVertex3d(graphFeature.cf[fs[1]].rep.x, graphFeature.cf[fs[1]].rep.y, graphFeature.cf[fs[1]].rep.z);
			}
			else
			{
				unsigned fs[] = { e_it->face_iter(0)->id() };
				glVertex3d(graphFeature.cf[fs[0]].rep.x, graphFeature.cf[fs[0]].rep.y, graphFeature.cf[fs[0]].rep.z);
				glVertex3d(graphFeature.ef[e_it->id()].rep.x, graphFeature.ef[e_it->id()].rep.y, graphFeature.ef[e_it->id()].rep.z);
			}
		}
		glEnd();
	}
	else
	{
		glColor3d(0, 0, 0);

		for (int i = 0; i < boundaryCurves.size(); i++)
		{
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j < boundaryCurves[i].size(); j++)
				glVertex3d(boundaryCurves[i][j].x, boundaryCurves[i][j].y, boundaryCurves[i][j].z);
			glEnd();
		}
	}
	glLineWidth(1.0);

}

//interaction visualization
void MeshSegment::drawInteraction()
{
	if (m_interaction->scribbleType == Interaction::ERASE || m_interaction->scribbleType == Interaction::REVEAL || m_interaction->scribbleType == Interaction::CONCEAL)
	{
		if (m_interaction->scribbleType == Interaction::REVEAL || m_interaction->scribbleType == Interaction::CONCEAL)
			glColor4d(0.4, 0.4, 0.4, 0.4);
		else
			glColor4d(1, 1, 0, 0.4);

		if (m_interaction->isScrolling == true)
		{
			if (m_interaction->scribbleType == Interaction::REVEAL)
				glColor4d(1, 0, 0, 0.4);
			else if (m_interaction->scribbleType == Interaction::CONCEAL)
				glColor4d(0, 0, 1, 0.4);
		}

		//draw scribbled curves;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBegin(GL_TRIANGLES);
		for (auto f = m_interaction->scribbleTriangleStrip.begin(); f != m_interaction->scribbleTriangleStrip.end(); f++)
		{
			auto f_it = myMesh->getFIter()[*f];
			MyMesh::VertexIter v[] = { f_it->vertex_iter(0), f_it->vertex_iter(1), f_it->vertex_iter(2) };
			for (unsigned j = 0; j < 3; j++)
			{
				glNormal3f(v[j]->normal().x, v[j]->normal().y, v[j]->normal().z);
				glVertex3f(v[j]->coordinate().x, v[j]->coordinate().y, v[j]->coordinate().z);
			}
		}
		for (auto s = m_interaction->scribbleTriangleStrips.begin(); s != m_interaction->scribbleTriangleStrips.end(); s++)
		{
			for (auto f = s->begin(); f != s->end(); f++)
			{
				auto f_it = myMesh->getFIter()[*f];
				MyMesh::VertexIter v[] = { f_it->vertex_iter(0), f_it->vertex_iter(1), f_it->vertex_iter(2) };
				for (unsigned j = 0; j < 3; j++)
				{
					glNormal3f(v[j]->normal().x, v[j]->normal().y, v[j]->normal().z);
					glVertex3f(v[j]->coordinate().x, v[j]->coordinate().y, v[j]->coordinate().z);
				}
			}
		}
		glEnd();

		//draw mouse location
		if (m_interaction->scribleRadius == 0)
		{
			double bsize = m_interaction->scribbleType == Interaction::ERASE ? 5 : m_interaction->brushSize;
			GLint	viewport[4];
			glGetIntegerv(GL_VIEWPORT, viewport);
			m_interaction->scribleRadius = bsize*1.2 / double(viewport[3]);
		}

		GLUquadric *quad = gluNewQuadric();

		glPushMatrix();
		glTranslated(m_interaction->hit.hittedNearVertex.x, m_interaction->hit.hittedNearVertex.y, m_interaction->hit.hittedNearVertex.z);

		glMatrixMode(GL_MODELVIEW);
		float modelview_matrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix);
		for (int i = 0; i < 12; i++)
		{
			if (i % 5 == 0) modelview_matrix[i] = 1;
			else modelview_matrix[i] = 0;
		}
		modelview_matrix[15] = 1.0;
		glLoadIdentity();
		glMultMatrixf(modelview_matrix);

		double eyeChange = modelview_matrix[14] / -2.0;
		glScalef(eyeChange, eyeChange, eyeChange);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		gluDisk(quad, 0., m_interaction->scribleRadius, 18, 1);
		glPopMatrix();
	}
	else if (m_interaction->scribbleType == Interaction::ADD)
	{
		//draw sketch curves on mesh
		glColor3d(0, 0, 0);
		glLineWidth(1.0);

		glBegin(GL_LINE_STRIP);
		for (unsigned i = 0; i<m_interaction->scribbleCurve.size(); i++)
		{
			glVertex3d(m_interaction->scribbleCurve[i].x, m_interaction->scribbleCurve[i].y, m_interaction->scribbleCurve[i].z);
		}
		glEnd();

		for (unsigned c = 0; c<m_interaction->scribbleCurves.size(); c++)
		{
			glBegin(GL_LINE_STRIP);
			for (unsigned i = 0; i<m_interaction->scribbleCurves[c].size(); i++)
			{
				glVertex3d(m_interaction->scribbleCurves[c][i].x, m_interaction->scribbleCurves[c][i].y, m_interaction->scribbleCurves[c][i].z);
			}
			glEnd();
		}


		//draw cross sign, indicating the mouse location
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();

		glTranslated(m_interaction->hit.hittedNearVertex.x, m_interaction->hit.hittedNearVertex.y, m_interaction->hit.hittedNearVertex.z);
		float modelview_matrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix);
		for (int i = 0; i<12; i++)
		{
			if (i % 5 == 0) modelview_matrix[i] = 1;
			else modelview_matrix[i] = 0;
		}
		modelview_matrix[15] = 1;
		glLoadIdentity();
		glMultMatrixf(modelview_matrix);
		double eyeChange = modelview_matrix[14] / -2.0;
		glScalef(eyeChange, eyeChange, eyeChange);

		glBegin(GL_LINES);
		glVertex2f(-0.015, 0);
		glVertex2f(+0.015, 0);
		glVertex2f(0, -0.015);
		glVertex2f(0, +0.015);
		glEnd();

		glPopMatrix();
	}
}

void MeshSegment::DrawCurve()
{
	if (myMesh == NULL) return;

	//the curves feed to our algorithm contains crestlines and user inputs.
	if (showCrestLine)
	{
		drawCrestLine(); drawUserSketch();
	} 	

	// draw all feature lines use a light grey color helping user to add new sketched 
	if (m_interaction->isShiftPress) 
		drawGroupFeatures();

	// the result curve network -- the boundary between patches.
	if (showSegmentationBoundary)
		drawSegBoundary();

	if (m_interaction->isShiftPress || m_interaction->isControlPress || m_interaction->isAltPress)
		drawInteraction();

}
void MeshSegment::DrawGraph()
{
	if (myMesh == NULL) return;

	if (showMesh)
	{
		if (showSegmentation)
			drawSegmentation();
	}
}