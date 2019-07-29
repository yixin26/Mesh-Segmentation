#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_CurveNetMaker.h"


class MyInteraction;
class MyGLCanvas;
class MeshSegment;

class CurveNetMaker : public QMainWindow
{
	Q_OBJECT

public:
	CurveNetMaker(QWidget *parent = Q_NULLPTR);

private:
	Ui::CurveNetMakerClass ui;


public slots:

	void OpenMeshFile();
	void SaveCurveFile();
	void AlphaChanged();
	void ShowFeatureLines();
	void SmoothBoundary();

public:

	MyInteraction* m_interaction;
	MeshSegment* m_meshSegment;

private:

	void ResetAll();
	void ConstructMesh(const char* fileName, const char* fileType);
	void ComputeCrestLine();
	void OverSegmentation();


	void keyPressEvent(QKeyEvent *e);
	void keyReleaseEvent(QKeyEvent *e);

	void wheelEvent(QWheelEvent *e);

	void resizeEvent(QResizeEvent* event);

};
