#include "CurveNetMaker.h"
#include <QtGui/QtGui>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <iostream>
#include <QtWidgets/QInputDialog>

#include "core/Interaction.h"
#include "core/Segmentation.h"
#include "MyGLCanvas.h"


CurveNetMaker::CurveNetMaker(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);


	//before everything get started, initialize the mesh, interaction and render tool.
	m_meshSegment = new MeshSegment();
	*(ui.openGLWidget->getMeshSegment()) = m_meshSegment;

	m_interaction = new MyInteraction();
	*(ui.openGLWidget->getInteractionUtils()) = m_interaction;
	*(m_meshSegment->getInteractionUtils()) = m_interaction;

	//ui.openGLWidget->initializeGL();
	ui.openGLWidget->pos_y = this->frameGeometry().height() - ui.openGLWidget->frameGeometry().height() - ui.menuBar->frameGeometry().height();

	ui.AlphaSpin->setVisible(false);
	ui.AlphaTxt->setVisible(false);
	ui.FeatureCheckBox->setVisible(false);
	ui.SmoothBoundaryCheckBox->setVisible(false);

}

void CurveNetMaker::ResetAll()
{
	delete m_meshSegment;
	delete m_interaction;

	m_meshSegment = new MeshSegment();
	*(ui.openGLWidget->getMeshSegment()) = m_meshSegment;

	m_interaction = new MyInteraction();
	*(m_meshSegment->getInteractionUtils()) = m_interaction;
	*(ui.openGLWidget->getInteractionUtils()) = m_interaction;

	ui.openGLWidget->initializeGL();
	ui.openGLWidget->setMouseTracking(true);
}

//first step: load a mesh, with extension of off,obj,ply...
// store to myMesh, m_meshSegment
void CurveNetMaker::OpenMeshFile()
{
	QString file_name = QFileDialog::getOpenFileName(this,
		tr("Open File"),
		"",
		"*.off;*.ply;*.obj;;*.off;;*.ply;;*.obj",
		0);

	if (!file_name.isNull())
	{
		char tFileName[400];
		strcpy(tFileName, file_name.toStdString().data());
		char FileType[400];
		for (int i = strlen(tFileName) - 1; i >= 0; i--) {
			if (tFileName[i] == '.') {
				i++;
				unsigned tsize = strlen(tFileName) - i;
				for (unsigned j = 0; j < tsize; j++) {
					FileType[j] = tFileName[i + j];
				}
				FileType[tsize] = '\0';
				break;
			}
		}

		ConstructMesh(file_name.toStdString().data(), FileType);
		ui.openGLWidget->paintGL();

		ComputeCrestLine();
		ui.openGLWidget->paintGL();

		OverSegmentation();
		ui.openGLWidget->paintGL();

		ui.AlphaSpin->setVisible(true);
		ui.AlphaTxt->setVisible(true);
		ui.FeatureCheckBox->setVisible(true);
		ui.SmoothBoundaryCheckBox->setVisible(true);
	}
}
void CurveNetMaker::ConstructMesh(const char* fileName, const char* fileType)
{
	ResetAll();

	m_meshSegment->m_filename = fileName;

	if (strcmp(fileType, "ply") == 0 || strcmp(fileType, "off") == 0 || strcmp(fileType, "obj") == 0) 
	{
		m_meshSegment->readTrimesh(fileName);
	}
	else
		return;

	*(m_interaction->getMyMesh()) = *(m_meshSegment->getMyMesh());
	m_interaction->initSelectionColor();

	srand(time(NULL));
}

void CurveNetMaker::ComputeCrestLine()
{
	m_meshSegment->computeCrestLine();
	m_meshSegment->computeMultiScaleFeatures();
}

void CurveNetMaker::OverSegmentation()
{
	//mitani's
	m_meshSegment->mitani_Watershed_Ringsize = 3;
	m_meshSegment->isAnisGeodesics_Watershed = true;
	m_meshSegment->useAllFeatureLine = true; //false if weak features are pruned
	m_meshSegment->gcIsMerge = true;
	m_meshSegment->overSegmentation();
	m_meshSegment->mergePartition();
}

void  CurveNetMaker::AlphaChanged()
{
	if (ui.AlphaSpin->value() == m_meshSegment->graphFeature.alpha) return;

	m_meshSegment->graphFeature.alpha = ui.AlphaSpin->value();
	m_meshSegment->edgeWeightParameter.clear();
	//OverSegmentation();
	m_meshSegment->updateGlabalAwardAlpha();
	m_meshSegment->mergePartition();
	ui.openGLWidget->paintGL();
}
void CurveNetMaker::ShowFeatureLines()
{
	m_meshSegment->showCrestLine = ui.FeatureCheckBox->isChecked();
	ui.openGLWidget->paintGL();
}
void CurveNetMaker::SmoothBoundary()
{
	m_meshSegment->autoBoundarySmooth = ui.SmoothBoundaryCheckBox->isChecked();
	m_meshSegment->boundarySmooth(m_meshSegment->autoBoundarySmooth);
	ui.openGLWidget->paintGL();
}
void CurveNetMaker::SaveCurveFile()
{
	m_meshSegment->saveSegmentation();

	QMessageBox message(QMessageBox::NoIcon, "Save File", QString("The files are save to the same directory. There are three files, including a per vertex labeling file(.label), two segmentation boundary files(.curve), with or without smoothing."));
	message.exec();
}

void CurveNetMaker::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_Q:
		m_meshSegment->crestline_scale = 1;
		break;
	case Qt::Key_W:
		m_meshSegment->crestline_scale = 2;
		break;
	case Qt::Key_E:
		m_meshSegment->crestline_scale = 3;
		break;
	}
	//if (e->key() < 4 && e->key() > 0)
	if (e->key() == Qt::Key_Q || e->key() == Qt::Key_W || e->key() == Qt::Key_E)
	{
		ComputeCrestLine();
		OverSegmentation();
		ui.openGLWidget->paintGL();
	}

	ui.openGLWidget->keyPressEvent(e);
}
void CurveNetMaker::keyReleaseEvent(QKeyEvent *e)
{
	ui.openGLWidget->keyReleaseEvent(e);
}

void CurveNetMaker::wheelEvent(QWheelEvent *e)
{
	//ui.openGLWidget->wheelEvent(e);

	if (ui.AlphaSpin->value() != m_meshSegment->graphFeature.alpha)
	{
		ui.AlphaSpin->setDisabled(true);
		//ui.AlphaSpin->setValue(m_meshSegment->graphFeature.alpha);
	}
}

void CurveNetMaker::resizeEvent(QResizeEvent* event)
{
	QMainWindow::resizeEvent(event);
	// Your code here
	//int width = this->frameGeometry().width();
	//int height = this->frameGeometry().height();
	int width = ui.centralWidget->frameGeometry().width();
	int height = ui.centralWidget->frameGeometry().height();
	ui.openGLWidget->setFixedSize(QSize(width, height));
	
	ui.openGLWidget->resizeGL(width, height);
	ui.openGLWidget->paintGL();
}