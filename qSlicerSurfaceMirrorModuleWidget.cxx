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

#include <vtkMRMLHierarchyNode.h>
#include <vtkMRMLAnnotationHierarchyNode.h>
#include <vtkMRMLAnnotationFiducialNode.h>


//	Qt includes
#include <QDebug>
#include <QList>
#include <QMessageBox>
#include "QPlane.h"
#include <QString>
#include <QColor>  
#include <QPalette>


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

#include <vtkCollection.h>
#include <vtkGeneralTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkReverseSense.h>
#include <vtkObject.h>
#include <vtkCollectionIterator.h>

#include <math.h>
#include <float.h>

#include <iostream> // for debugging

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

	QPalette pal;
	pal.setColor(QPalette::Text, QColor( QString("red") ));
	d->hintText->setPalette(pal);
	d->hintText->setText("");

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

	QPalette pal;
	pal.setColor(QPalette::Text, QColor( QString("red") ));
	d->hintText->setPalette(pal);
	d->hintText->setText("");

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

	QPalette pal;
	pal.setColor(QPalette::Text, QColor( QString("red") ));
	d->hintText->setPalette(pal);

	//	obtain the source model for clipping
  	qSlicerApplication *app = qSlicerApplication::application();
	vtkMRMLScene *mrmlScene = app->mrmlScene();

	vtkMRMLNode *sourceNode = d->mirrorNodeComboBox->currentNode();
	if (sourceNode) {
		this->planeWidget->On();
		this->renderWindow->Render();
		this->timeOfMirror++;
	  	d->createButton->setEnabled(false);
		d->hideButton->setEnabled(true);
	  	d->mirrorButton->setEnabled(this->hasPlane);

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
		QString resultName = tr("Mirrored Part_");
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
		d->hintText->setText("");

	} else d->hintText->setText("Please select a model to mirror!");
}

void qSlicerSurfaceMirrorModuleWidget::createPlaneWithAnnotation()
{
	Q_D(qSlicerSurfaceMirrorModuleWidget);

	QPalette pal;
	pal.setColor(QPalette::Text, QColor( QString("red") ));
	d->hintText->setPalette(pal);
	d->hintText->setText("");

	vtkMRMLAnnotationHierarchyNode *annotationHierarchy =
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(d->annotationNodeComboBox->currentNode());

	if (annotationHierarchy) {
		this->hasPlane = true;
		this->iren->Initialize();
		this->renderWindow->Render();
		this->iren->Start();

		d->createButton->setEnabled(false);
		d->hideButton->setEnabled(true);
	  	d->mirrorButton->setEnabled(this->hasPlane);

		vtkCollection *annotationChildrenCollection = vtkCollection::New();
		annotationHierarchy->GetAllChildren(annotationChildrenCollection);

		int maxNumCollection = annotationChildrenCollection->GetNumberOfItems();
		int numCollection = 0;
		double **coord = new double*[maxNumCollection];
		for (int i = 0; i < maxNumCollection; i++)
	    {
	        coord[i]  = new double[3];
	    }
		double minXYZ[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
		double maxXYZ[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

		vtkCollectionIterator *iterator = annotationChildrenCollection->NewIterator();
		iterator->InitTraversal();
		vtkObject *object = NULL;
		vtkMRMLAnnotationFiducialNode *fiducial = NULL;

		while ( !(iterator->IsDoneWithTraversal()) ) {
			std::clog << "Start Traversal" << endl;
			object = iterator->GetCurrentObject();
			iterator->GoToNextItem();

			fiducial = vtkMRMLAnnotationFiducialNode::SafeDownCast(object);

			if (fiducial) {
				std::clog << "Fiducial Success" << endl;
				if (fiducial->GetFiducialCoordinates(coord[numCollection]) ){
					if (coord[numCollection][0] < minXYZ[0]) minXYZ[0] = coord[numCollection][0];
					if (coord[numCollection][0] > maxXYZ[0]) maxXYZ[0] = coord[numCollection][0];
					if (coord[numCollection][1] < minXYZ[1]) minXYZ[1] = coord[numCollection][1];
					if (coord[numCollection][1] > maxXYZ[1]) maxXYZ[1] = coord[numCollection][1];
					if (coord[numCollection][2] < minXYZ[2]) minXYZ[2] = coord[numCollection][2];
					if (coord[numCollection][2] > maxXYZ[2]) maxXYZ[2] = coord[numCollection][2];

					numCollection ++;	
	 			} else std::clog << "No Coord" << endl;

			} else std::clog << "Not Fiducial" << endl;
		}

		if (numCollection >= 3){
			double *parameter = new double[3];
			if (this->getRegressivePlanePara(numCollection, coord, parameter)) {
				std::clog << parameter[0] << ", " << parameter[1] << ", " << parameter[2] << endl;

				double newOrigin[3], newPoint1[3], newPoint2[3];

				// newPoint1 denotes the minimal point
				if (parameter[0] * minXYZ[0] + parameter[1] * minXYZ[1] + parameter[2] < minXYZ[2]) {
					newPoint1[0] = minXYZ[0];
					newPoint1[1] = minXYZ[1];
					newPoint1[2] = parameter[0] * minXYZ[0] + parameter[1] * minXYZ[1] + parameter[2];
				} else if ( (minXYZ[2] - parameter[0] * minXYZ[0] - parameter[2]) / parameter[1] < minXYZ[1]) {
					newPoint1[0] = minXYZ[0];
					newPoint1[2] = minXYZ[2];
					newPoint1[1] = (minXYZ[2] - parameter[0] * minXYZ[0] - parameter[2]) / parameter[1];
				} else {
					newPoint1[1] = minXYZ[1];
					newPoint1[2] = minXYZ[2];
					newPoint1[0] = (minXYZ[2] - parameter[1] * minXYZ[1] - parameter[2]) / parameter[0];
				}


				// newPoint2 denotes the maximal point
				if (parameter[0] * maxXYZ[0] + parameter[1] * maxXYZ[1] + parameter[2] > maxXYZ[2]) {
					newPoint2[0] = maxXYZ[0];
					newPoint2[1] = maxXYZ[1];
					newPoint2[2] = parameter[0] * maxXYZ[0] + parameter[1] * maxXYZ[1] + parameter[2];
				} else if ( (maxXYZ[2] - parameter[0] * maxXYZ[0] - parameter[2]) / parameter[1] > maxXYZ[1]) {
					newPoint2[0] = maxXYZ[0];
					newPoint2[2] = maxXYZ[2];
					newPoint2[1] = (maxXYZ[2] - parameter[0] * maxXYZ[0] - parameter[2]) / parameter[1];
				} else {
					newPoint2[1] = maxXYZ[1];
					newPoint2[2] = maxXYZ[2];
					newPoint2[0] = (maxXYZ[2] - parameter[1] * maxXYZ[1] - parameter[2]) / parameter[0];
				}

				double V[3];
				vtkMath::Subtract(newPoint2, newPoint1, V);
				double normal[3] = {parameter[0], parameter[1], -1.0};
				double double_d[3];
				vtkMath::Cross(normal, V, double_d);
				double_d[0] = double_d[0] / vtkMath::Norm(double_d) * vtkMath::Norm(V) / 2;
				double_d[1] = double_d[1] / vtkMath::Norm(double_d) * vtkMath::Norm(V) / 2;
				double_d[2] = double_d[2] / vtkMath::Norm(double_d) * vtkMath::Norm(V) / 2;
				double M[3];
				vtkMath::Add(newPoint2, newPoint1, M);
				M[0] /= 2;
				M[1] /= 2;
				M[2] /= 2;
				vtkMath::Add(M, double_d, newOrigin);

				this->planeWidget->SetOrigin(newOrigin);
				this->planeWidget->SetPoint1(newPoint1);
				this->planeWidget->SetPoint2(newPoint2);

				this->hasPlane = true;
				d->hideButton->setText(tr("Hide Plane"));
				this->planeWidget->On();
			}
		} else d->hintText->setText("At least three Fiducials needed!");

	} else {
		d->hintText->setText("Please select a list of Fiducials to create plane!");
	}
}

bool qSlicerSurfaceMirrorModuleWidget::getRegressivePlanePara(int numCollection, double *coord[3], double* parameter) {
	Q_D(qSlicerSurfaceMirrorModuleWidget);

	QPalette pal;
	pal.setColor(QPalette::Text, QColor( QString("red") ));
	d->hintText->setPalette(pal);

	for (int i = 0; i < numCollection; i++) {
		std::clog << coord[i][0] << ", " << coord[i][1] << ", " << coord[i][2] << endl;
	}

	double Y[3] = {0.0, 0.0, 0.0};
	parameter[0] = 0.0, parameter[1] = 0.0, parameter[2] = 0.0;

	double *Matrix[3], *IMatrix[3];
	for (int i = 0; i < 3; i++)
    {
        Matrix[i]  = new double[3];
        IMatrix[i] = new double[3];
    }
    for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++)
        {
            *(Matrix[i] + j) = 0.0;
        }
    }

    for (int i = 0; i < numCollection; i++)
    {
    	Matrix[0][0] += coord[i][0] * coord[i][0];
    	Matrix[1][1] += coord[i][1] * coord[i][1];
    	Matrix[1][0] += coord[i][1] * coord[i][0];
    	Matrix[0][1] += coord[i][1] * coord[i][0];
    	Matrix[2][0] += coord[i][0];
    	Matrix[0][2] += coord[i][0];
    	Matrix[2][1] += coord[i][1];
    	Matrix[1][2] += coord[i][1];
    	Matrix[2][2] += 1;
    	Y[0] += coord[i][0] * coord[i][2];
    	Y[1] += coord[i][1] * coord[i][2];
    	Y[2] += coord[i][2];
    }

    double double_d = this->Determinant(Matrix, 3);
    if (std::abs(double_d) < 0.0001)
    {
    	d->hintText->setText("Plane creating failed: Singular Matrix.");
        return false;
    }

    Inverse(Matrix, IMatrix, 3, double_d);
    for (int i = 0; i < 3; i++)
    {
        parameter[0] += *(IMatrix[0] + i) * Y[i];
        parameter[1] += *(IMatrix[1] + i) * Y[i];
        parameter[2] += *(IMatrix[2] + i) * Y[i];
    }

    for (int i = 0; i < 3; i++)
    {
        delete[] Matrix[i];
        delete[] IMatrix[i];
    }

	return true;
}

double qSlicerSurfaceMirrorModuleWidget::Determinant(double* matrix[], int n)  
{  
    double result = 0,temp;  
    int i;  
    if (n == 1)  
        result = (*matrix[0]);  
    else  
    {  
        for (i = 0; i < n; i++)  
        {  
            temp = AlCo(matrix, n, n-1, i);  
            result += ( *(matrix[n-1] + i) ) * temp;  
        }  
    }  
    return result;  
}  

double qSlicerSurfaceMirrorModuleWidget::AlCo(double* matrix[], int jie, int row, int column)  
{  
    double result; 
    if ((row+column) % 2 == 0) 
        result = Cofactor(matrix, jie, row, column);  
    else result = (-1) * Cofactor(matrix, jie, row, column); 
    return result;  
}  

double qSlicerSurfaceMirrorModuleWidget::Cofactor(double* matrix[], int jie, int row, int column)  
{  
    double result;  
    int i, j;  
    double* smallmatr[9];  
    for(i = 0; i < jie-1; i++)  
        smallmatr[i] = new double[jie - 1];
    for(i = 0; i < row; i++)  
        for(j = 0; j < column; j++)  
            *(smallmatr[i] + j) = *(matrix[i] + j);  
    for(i = row; i < jie-1; i++)  
        for(j = 0; j < column; j++)  
            *(smallmatr[i] + j) = *(matrix[i+1] + j);  
    for(i = 0; i < row; i++)  
        for(j = column; j < jie-1; j++)  
            *(smallmatr[i] + j) = *(matrix[i] + j + 1);  
    for(i = row; i < jie - 1; i++)  
        for(j = column; j < jie - 1; j++)  
            *(smallmatr[i] + j) = *(matrix[i + 1] + j + 1);  
    result = Determinant(smallmatr, jie - 1); 
    for(i = 0; i < jie - 1; i++)
        delete[] smallmatr[i];
    return result;   
}

void qSlicerSurfaceMirrorModuleWidget::Inverse(double *matrix1[], double *matrix2[], int n, double double_d) 
{ 
    int i, j; 
    for(i = 0; i < n; i++) 
        matrix2[i] = (double *)malloc( n * sizeof(double) ); 
    for(i = 0; i < n; i++) 
        for(j = 0; j < n; j++) 
            *(matrix2[j] + i) = (AlCo(matrix1, n, i, j) / double_d); 
} 



void qSlicerSurfaceMirrorModuleWidget::setup()
{
	Q_D(qSlicerSurfaceMirrorModuleWidget);
	d->setupUi(this);
	this->Superclass::setup();

	QObject::connect(d->createButton, SIGNAL(clicked()), this, SLOT(createPlane()));
	QObject::connect(d->hideButton, SIGNAL(clicked()), this, SLOT(hidePlane()));
	QObject::connect(d->mirrorButton, SIGNAL(clicked()), this, SLOT(mirror()));
	QObject::connect(d->annotationCreateButton, SIGNAL(clicked()), this, SLOT(createPlaneWithAnnotation()));
}
