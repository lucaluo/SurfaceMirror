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

#ifndef __qSlicerSurfaceMirrorModuleWidget_h
#define __qSlicerSurfaceMirrorModuleWidget_h

//	SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"
#include <qSlicerApplication.h>
#include "qSlicerSurfaceMirrorModuleExport.h"

//	Qt includes
#include "QPlane.h"

//	vtk includes
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include "vtkPlaneExtend.h"
#include <vtkPlaneWidget.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>

#include <vtkPolyLine.h>

#include <vtkGeneralTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkReverseSense.h>

#include <math.h>

class qSlicerSurfaceMirrorModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SURFACEMIRROR_EXPORT qSlicerSurfaceMirrorModuleWidget : 
	public qSlicerAbstractModuleWidget
{
	Q_OBJECT

public:
	typedef qSlicerAbstractModuleWidget Superclass;
	qSlicerSurfaceMirrorModuleWidget(QWidget *parent=0);
	virtual ~qSlicerSurfaceMirrorModuleWidget();

public slots:

protected:
	vtkRenderWindowInteractor* iren;
	vtkRenderWindow* renderWindow;
	vtkRenderer* renderer;
	vtkSmartPointer<vtkPlaneWidget> planeWidget;

	vtkSmartPointer<vtkPolyData> sourcePD;

	QList<vtkSmartPointer<vtkPolyData> > resultList;
	int timeOfMirror;
  	bool hasPlane;

	QScopedPointer<qSlicerSurfaceMirrorModuleWidgetPrivate> d_ptr;
	virtual void setup();

protected slots:
	void createPlane();
	void hidePlane();
	void mirror();
	void createPlaneWithAnnotation();

	bool getRegressivePlanePara(int numCollection, double *coord[3], double* parameter);
	double Determinant(double* matrix[], int n);
	double AlCo(double* matrix[], int jie, int row, int column);
	double Cofactor(double* matrix[], int jie, int row, int column);
	void Inverse(double *matrix1[], double *matrix2[], int n, double d);

private:
	Q_DECLARE_PRIVATE(qSlicerSurfaceMirrorModuleWidget);
	Q_DISABLE_COPY(qSlicerSurfaceMirrorModuleWidget);
};

#endif
