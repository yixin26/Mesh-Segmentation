/********************************************************************************
** Form generated from reading UI file 'CurveNetMaker.ui'
**
** Created by: Qt User Interface Compiler version 5.12.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CURVENETMAKER_H
#define UI_CURVENETMAKER_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QWidget>
#include "myglcanvas.h"

QT_BEGIN_NAMESPACE

class Ui_CurveNetMakerClass
{
public:
    QAction *OpenMesh;
    QAction *SaveSegmentation;
    QAction *OpenMeshFile;
    QWidget *centralWidget;
    MyGLCanvas *openGLWidget;
    QCheckBox *FeatureCheckBox;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *AlphaTxt;
    QDoubleSpinBox *AlphaSpin;
    QCheckBox *SmoothBoundaryCheckBox;
    QMenuBar *menuBar;
    QMenu *menuOpen;
    QMenu *menuSave;

    void setupUi(QMainWindow *CurveNetMakerClass)
    {
        if (CurveNetMakerClass->objectName().isEmpty())
            CurveNetMakerClass->setObjectName(QString::fromUtf8("CurveNetMakerClass"));
        CurveNetMakerClass->resize(1600, 1000);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(CurveNetMakerClass->sizePolicy().hasHeightForWidth());
        CurveNetMakerClass->setSizePolicy(sizePolicy);
        CurveNetMakerClass->setMinimumSize(QSize(1600, 1000));
        CurveNetMakerClass->setMaximumSize(QSize(1600, 1000));
        CurveNetMakerClass->setStyleSheet(QString::fromUtf8(""));
        OpenMesh = new QAction(CurveNetMakerClass);
        OpenMesh->setObjectName(QString::fromUtf8("OpenMesh"));
        SaveSegmentation = new QAction(CurveNetMakerClass);
        SaveSegmentation->setObjectName(QString::fromUtf8("SaveSegmentation"));
        OpenMeshFile = new QAction(CurveNetMakerClass);
        OpenMeshFile->setObjectName(QString::fromUtf8("OpenMeshFile"));
        centralWidget = new QWidget(CurveNetMakerClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        sizePolicy.setHeightForWidth(centralWidget->sizePolicy().hasHeightForWidth());
        centralWidget->setSizePolicy(sizePolicy);
        openGLWidget = new MyGLCanvas(centralWidget);
        openGLWidget->setObjectName(QString::fromUtf8("openGLWidget"));
        openGLWidget->setGeometry(QRect(0, 0, 1600, 930));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(openGLWidget->sizePolicy().hasHeightForWidth());
        openGLWidget->setSizePolicy(sizePolicy1);
        FeatureCheckBox = new QCheckBox(centralWidget);
        FeatureCheckBox->setObjectName(QString::fromUtf8("FeatureCheckBox"));
        FeatureCheckBox->setGeometry(QRect(570, 942, 131, 21));
        QSizePolicy sizePolicy2(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(FeatureCheckBox->sizePolicy().hasHeightForWidth());
        FeatureCheckBox->setSizePolicy(sizePolicy2);
        FeatureCheckBox->setLayoutDirection(Qt::LeftToRight);
        FeatureCheckBox->setChecked(false);
        layoutWidget = new QWidget(centralWidget);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(740, 943, 111, 40));
        horizontalLayout = new QHBoxLayout(layoutWidget);
        horizontalLayout->setSpacing(0);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetFixedSize);
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        AlphaTxt = new QLabel(layoutWidget);
        AlphaTxt->setObjectName(QString::fromUtf8("AlphaTxt"));
        QSizePolicy sizePolicy3(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(AlphaTxt->sizePolicy().hasHeightForWidth());
        AlphaTxt->setSizePolicy(sizePolicy3);
        AlphaTxt->setMinimumSize(QSize(50, 0));
        AlphaTxt->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(AlphaTxt);

        AlphaSpin = new QDoubleSpinBox(layoutWidget);
        AlphaSpin->setObjectName(QString::fromUtf8("AlphaSpin"));
        sizePolicy.setHeightForWidth(AlphaSpin->sizePolicy().hasHeightForWidth());
        AlphaSpin->setSizePolicy(sizePolicy);
        AlphaSpin->setMaximum(10.000000000000000);
        AlphaSpin->setSingleStep(0.050000000000000);
        AlphaSpin->setValue(0.500000000000000);

        horizontalLayout->addWidget(AlphaSpin);

        SmoothBoundaryCheckBox = new QCheckBox(centralWidget);
        SmoothBoundaryCheckBox->setObjectName(QString::fromUtf8("SmoothBoundaryCheckBox"));
        SmoothBoundaryCheckBox->setGeometry(QRect(930, 942, 111, 21));
        SmoothBoundaryCheckBox->setChecked(true);
        CurveNetMakerClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(CurveNetMakerClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1600, 23));
        sizePolicy.setHeightForWidth(menuBar->sizePolicy().hasHeightForWidth());
        menuBar->setSizePolicy(sizePolicy);
        menuOpen = new QMenu(menuBar);
        menuOpen->setObjectName(QString::fromUtf8("menuOpen"));
        menuSave = new QMenu(menuBar);
        menuSave->setObjectName(QString::fromUtf8("menuSave"));
        CurveNetMakerClass->setMenuBar(menuBar);

        menuBar->addAction(menuOpen->menuAction());
        menuBar->addAction(menuSave->menuAction());
        menuOpen->addAction(OpenMeshFile);
        menuSave->addAction(SaveSegmentation);

        retranslateUi(CurveNetMakerClass);
        QObject::connect(OpenMeshFile, SIGNAL(triggered()), CurveNetMakerClass, SLOT(OpenMeshFile()));
        QObject::connect(SaveSegmentation, SIGNAL(triggered()), CurveNetMakerClass, SLOT(SaveCurveFile()));
        QObject::connect(AlphaSpin, SIGNAL(valueChanged(double)), CurveNetMakerClass, SLOT(AlphaChanged()));
        QObject::connect(FeatureCheckBox, SIGNAL(clicked()), CurveNetMakerClass, SLOT(ShowFeatureLines()));
        QObject::connect(SmoothBoundaryCheckBox, SIGNAL(clicked()), CurveNetMakerClass, SLOT(SmoothBoundary()));

        QMetaObject::connectSlotsByName(CurveNetMakerClass);
    } // setupUi

    void retranslateUi(QMainWindow *CurveNetMakerClass)
    {
        CurveNetMakerClass->setWindowTitle(QApplication::translate("CurveNetMakerClass", "CurveNetMaker", nullptr));
        OpenMesh->setText(QApplication::translate("CurveNetMakerClass", "Mesh File", nullptr));
        SaveSegmentation->setText(QApplication::translate("CurveNetMakerClass", "Segmentation", nullptr));
        OpenMeshFile->setText(QApplication::translate("CurveNetMakerClass", "Mesh", nullptr));
        FeatureCheckBox->setText(QApplication::translate("CurveNetMakerClass", "Show Feature Line", nullptr));
        AlphaTxt->setText(QApplication::translate("CurveNetMakerClass", "Alpha ", nullptr));
        SmoothBoundaryCheckBox->setText(QApplication::translate("CurveNetMakerClass", "Smooth Boundary", nullptr));
        menuOpen->setTitle(QApplication::translate("CurveNetMakerClass", "Open", nullptr));
        menuSave->setTitle(QApplication::translate("CurveNetMakerClass", "Save", nullptr));
    } // retranslateUi

};

namespace Ui {
    class CurveNetMakerClass: public Ui_CurveNetMakerClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CURVENETMAKER_H
