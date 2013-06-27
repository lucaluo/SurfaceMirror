/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

//	MRML includes
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLModelStorageNode.h>
#include <vtkMRMLNode.h>
#include <vtkMRMLScene.h>

//	Qt includes
#include <QDebug>
#include <QList>
#include <QMessageBox>
#include "QPlane.h"
#include <QString>

//	SlicerQt includes
#include <qMRMLThreeDView.h>
#include <qMRMLThreeDWidget.h>
#include <qSlicerApplication.h>
#include <qSlicerIOManager.h>
#include <qSlicerLayoutManager.h>
#include "qSlicerSurfaceMirrorModuleWidget.h"
#include "ui_qSlicerSurfaceMirrorModule.h"

//	vtk includes
#include <vtkActor.h>
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkMath.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include "vtkPlaneExtend.h"
#include <vtkPlaneSource.h>
#include <vtkPlaneWidget.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <vtkGeneralTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkReverseSense.h>

#include <math.h>

///	\ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSurfaceMirrorModuleWidgetPrivate: public Ui_qSlicerSurfaceMirrorModule
{
	Q_DECLARE_PUBLIC(qSlicerSurfaceMirrorModuleWidget);

public:
	qSlicerSurfaceMirrorModuleWidgetPrivate(qSlicerSurfaceMirrorModuleWidget& object);

protected:
	qSlicerSurfaceMirrorModuleWidget* const q_ptr;
};

//	qSlicerSurfaceMirrorModuleWidgetPrivate methods

qSlicerSurfaceMirrorModuleWidgetPrivate::
	qSlicerSurfaceMirrorModuleWidgetPrivate(qSlicerSurfaceMirrorModuleWidget& object)
	 : q_ptr(&object)
{
}

//	qSlicerSurfaceMirrorModuleWidget methods

qSlicerSurfaceMirrorModuleWidget::qSlicerSurfaceMirrorModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSurfaceMirrorModuleWidgetPrivate(*this) )
{
  //  return the reference to the renderWindow of the Slicer application singleton 
  //  and those to its renderer and renderWindowInteractor
	qSlicerApplication *app = qSlicerApplication::application();
	qSlicerLayoutManager *layoutManager = app->layoutManager();
	qMRMLThreeDWidget *threeDWidget = layoutManager->threeDWidget(0);
	qMRMLThreeDView *threeDView = threeDWidget->threeDView();
	this->renderWindow = threeDView->renderWindow();
	this->renderer = this->renderWindow->GetRenderers()->GetFirstRenderer();
	this->iren = this->renderWindow->GetInteractor();
	this->planeWidget = vtkSmartPointer<vtkPlaneWidget>::New();
	this->planeWidget->SetInteractor(this->iren);
	this->planeWidget->Off();
	this->planeWidget->SetOrigin(0, 0, 0);
	this->planeWidget->SetPoint1(0, 0, 200);
	this->planeWidget->SetPoint2(50, 0, 0);
	this->planeWidget->SetHandleSize(0.01);
	this->planeWidget->SetRepresentationToSurface();
	this->planeWidget->PlaceWidget();

	this->timeOfMirror = 0;
  	this->hasPlane = false;
}

qSlicerSurfaceMirrorModuleWidget::~qSlicerSurfaceMirrorModuleWidget()
{
}

void qSlicerSurfaceMirrorModuleWidget::createPlane()
{
	Q_D(qSlicerSurfaceMirrorModuleWidget);
	this->hasPlane = true;
	this->iren->Initialize();
	this->planeWidget->On();
	this->renderWindow->Render();
	this->iren->Start();

	d->createButton->setEnabled(false);
	d->hideButton->setEnabled(true);
  	d->mirrorButton->setEnabled(this->hasPlane);
}


void qSlicerSurfaceMirrorModuleWidget::hidePlane()
{
	Q_D(qSlicerSurfaceMirrorModuleWidget);

	if (this->hasPlane == true) 
	{
		this->hasPlane = false;
		d->hideButton->setText(tr("Show Plane"));
		this->planeWidget->Off();
	}
	else {
		this->hasPlane = true;
		d->hideButton->setText(tr("Hide Plane"));
		this->planeWidget->On();
	}
	d->createButton->setEnabled(false);
  	d->hideButton->setEnabled(true);
  	d->mirrorButton->setEnabled(this->hasPlane);
  	this->renderWindow->Render();
}


void qSlicerSurfaceMirrorModuleWidget::mirror()
{
	Q_D(qSlicerSurfaceMirrorModuleWidget);

	this->planeWidget->On();
	this->renderWindow->Render();
	this->timeOfMirror++;
  	d->createButton->setEnabled(false);
	d->hideButton->setEnabled(true);
  	d->mirrorButton->setEnabled(this->hasPlane);

	//	obtain the source model for clipping
  	qSlicerApplication *app = qSlicerApplication::application();
	vtkMRMLScene *mrmlScene = app->mrmlScene();
	QString nodeID;
	nodeID = d->mirrorNodeComboBox->currentNodeID();

	vtkMRMLNode *sourceNode = d->mirrorNodeComboBox->currentNode();
	vtkSmartPointer<vtkMRMLModelNode> sourceModel =
		vtkSmartPointer<vtkMRMLModelNode>::New();
	sourceModel->Copy(sourceNode);
	
	vtkSmartPointer<vtkPolyData> sourcePD = vtkSmartPointer<vtkPolyData>::New();
	sourcePD->DeepCopy(sourceModel->GetPolyData());
	vtkSmartPointer<vtkTransformPolyDataFilter> mirrorer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkGeneralTransform> transformer = vtkSmartPointer<vtkGeneralTransform>::New();
	vtkSmartPointer<vtkReverseSense> reverser = vtkSmartPointer<vtkReverseSense>::New();

	mirrorer->SetInput(sourcePD);

	double *origin = this->planeWidget->GetOrigin();
	double *point1 = this->planeWidget->GetPoint1();
	double *point2 = this->planeWidget->GetPoint2();
	double *center = this->planeWidget->GetCenter();
	double *normal = this->planeWidget->GetNormal();
	double *rotAxis = new double[3];
	double *xAxis = new double[3];
	double rotAngle;

	// Make transfromation
	*xAxis = 1;
	*(xAxis+1) = 0;
	*(xAxis+2) = 0;

	vtkMath::Cross(normal, xAxis, rotAxis);
	rotAngle = std::atan( std::sqrt( *(normal+2) * *(normal+2) + *(normal+1) * *(normal+1) ) / *(normal) )  
				/ vtkMath::Pi() * 180;

	transformer->Translate(*(center), *(center+1), *(center+2));
	transformer->RotateWXYZ( -rotAngle  ,rotAxis);	
	transformer->Scale(-1, 1, 1);
	transformer->RotateWXYZ( rotAngle  ,rotAxis);	
	transformer->Translate(-*(center), -*(center+1), -*(center+2));

	mirrorer->SetTransform(transformer);

	mirrorer->Update();

	// Reverse the normal vector
	reverser->SetInput(mirrorer->GetOutput());
	reverser->ReverseNormalsOn();
	reverser->Update();


	//	get the output and the clipped output
  	vtkSmartPointer<vtkPolyData> resultPD = vtkSmartPointer<vtkPolyData>::New();
	resultPD->DeepCopy(reverser->GetOutput());
	this->resultList.append(resultPD);

	//	display and store the result model
	vtkSmartPointer<vtkMRMLModelNode> resultModel =
		vtkSmartPointer<vtkMRMLModelNode>::New();
	QString resultName = tr("Result Part_");
	resultName.append(QString::number(this->timeOfMirror));
	resultModel->SetName(resultName.toLatin1().data());
	resultModel->SetAndObservePolyData(this->resultList.at(this->timeOfMirror - 1));
	mrmlScene->SaveStateForUndo();
	resultModel->SetScene(mrmlScene);

	vtkSmartPointer<vtkMRMLModelDisplayNode> resultDisplay =
		vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
	vtkSmartPointer<vtkMRMLModelStorageNode> resultStorage =
		vtkSmartPointer<vtkMRMLModelStorageNode>::New();
	resultDisplay->SetScene(mrmlScene);
	resultStorage->SetScene(mrmlScene);
	resultDisplay->SetInputPolyData(resultModel->GetPolyData());
	resultDisplay->SetColor(1.0, 0.0, 0.0);
	resultStorage->SetFileName(resultName.toLatin1().data());
	mrmlScene->AddNode(resultDisplay);
	mrmlScene->AddNode(resultStorage);
	resultModel->SetAndObserveDisplayNodeID(resultDisplay->GetID());
	resultModel->SetAndObserveStorageNodeID(resultStorage->GetID());

	mrmlScene->AddNode(resultModel);
}

void qSlicerSurfaceMirrorModuleWidget::setup()
{
	Q_D(qSlicerSurfaceMirrorModuleWidget);
	d->setupUi(this);
	this->Superclass::setup();

	QObject::connect(d->createButton, SIGNAL(clicked()), this, SLOT(createPlane()));
	QObject::connect(d->hideButton, SIGNAL(clicked()), this, SLOT(hidePlane()));
	QObject::connect(d->mirrorButton, SIGNAL(clicked()), this, SLOT(mirror()));
}
