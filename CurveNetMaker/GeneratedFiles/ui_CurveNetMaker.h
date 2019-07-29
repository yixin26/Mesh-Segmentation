/********************************************************************************
** Form generated from reading UI file 'CurveNetMaker.ui'
**
** Created by: Qt User Interface Compiler version 5.12.4
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
#include <QtWidgets/QVBoxLayout>
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
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QCheckBox *FeatureCheckBox;
    QHBoxLayout *horizontalLayout;
    QLabel *AlphaTxt;
    QDoubleSpinBox *AlphaSpin;
    QLabel *label;
    QCheckBox *SmoothBoundaryCheckBox;
    QMenuBar *menuBar;
    QMenu *menuOpen;
    QMenu *menuSave;

    void setupUi(QMainWindow *CurveNetMakerClass)
    {
        if (CurveNetMakerClass->objectName().isEmpty())
            CurveNetMakerClass->setObjectName(QString::fromUtf8("CurveNetMakerClass"));
        CurveNetMakerClass->resize(1600, 1000);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(CurveNetMakerClass->sizePolicy().hasHeightForWidth());
        CurveNetMakerClass->setSizePolicy(sizePolicy);
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
        sizePolicy.setHeightForWidth(openGLWidget->sizePolicy().hasHeightForWidth());
        openGLWidget->setSizePolicy(sizePolicy);
        openGLWidget->setAutoFillBackground(false);
        verticalLayoutWidget = new QWidget(centralWidget);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 10, 202, 80));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        FeatureCheckBox = new QCheckBox(verticalLayoutWidget);
        FeatureCheckBox->setObjectName(QString::fromUtf8("FeatureCheckBox"));
        sizePolicy.setHeightForWidth(FeatureCheckBox->sizePolicy().hasHeightForWidth());
        FeatureCheckBox->setSizePolicy(sizePolicy);
        FeatureCheckBox->setMinimumSize(QSize(200, 0));
        FeatureCheckBox->setLayoutDirection(Qt::LeftToRight);
        FeatureCheckBox->setChecked(false);

        verticalLayout->addWidget(FeatureCheckBox);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(10);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetDefaultConstraint);
        AlphaTxt = new QLabel(verticalLayoutWidget);
        AlphaTxt->setObjectName(QString::fromUtf8("AlphaTxt"));
        sizePolicy.setHeightForWidth(AlphaTxt->sizePolicy().hasHeightForWidth());
        AlphaTxt->setSizePolicy(sizePolicy);
        AlphaTxt->setMinimumSize(QSize(0, 0));
        AlphaTxt->setMaximumSize(QSize(60, 16777215));
        AlphaTxt->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(AlphaTxt);

        AlphaSpin = new QDoubleSpinBox(verticalLayoutWidget);
        AlphaSpin->setObjectName(QString::fromUtf8("AlphaSpin"));
        sizePolicy.setHeightForWidth(AlphaSpin->sizePolicy().hasHeightForWidth());
        AlphaSpin->setSizePolicy(sizePolicy);
        AlphaSpin->setMaximumSize(QSize(60, 16777215));
        AlphaSpin->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);
        AlphaSpin->setMaximum(10.000000000000000);
        AlphaSpin->setSingleStep(0.050000000000000);
        AlphaSpin->setValue(0.500000000000000);

        horizontalLayout->addWidget(AlphaSpin);

        label = new QLabel(verticalLayoutWidget);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);


        verticalLayout->addLayout(horizontalLayout);

        SmoothBoundaryCheckBox = new QCheckBox(verticalLayoutWidget);
        SmoothBoundaryCheckBox->setObjectName(QString::fromUtf8("SmoothBoundaryCheckBox"));
        sizePolicy.setHeightForWidth(SmoothBoundaryCheckBox->sizePolicy().hasHeightForWidth());
        SmoothBoundaryCheckBox->setSizePolicy(sizePolicy);
        SmoothBoundaryCheckBox->setMinimumSize(QSize(200, 0));
        SmoothBoundaryCheckBox->setChecked(true);

        verticalLayout->addWidget(SmoothBoundaryCheckBox);

        CurveNetMakerClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(CurveNetMakerClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1600, 18));
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
        label->setText(QString());
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
