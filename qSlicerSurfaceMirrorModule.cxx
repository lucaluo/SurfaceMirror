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

// Qt includes
#include <QtPlugin>

// SurfaceMirror Logic includes
#include <vtkSlicerSurfaceMirrorLogic.h>

// SurfaceMirror includes
#include "qSlicerSurfaceMirrorModule.h"
#include "qSlicerSurfaceMirrorModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerSurfaceMirrorModule, qSlicerSurfaceMirrorModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSurfaceMirrorModulePrivate
{
public:
  qSlicerSurfaceMirrorModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSurfaceMirrorModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSurfaceMirrorModulePrivate
::qSlicerSurfaceMirrorModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSurfaceMirrorModule methods

//-----------------------------------------------------------------------------
qSlicerSurfaceMirrorModule
::qSlicerSurfaceMirrorModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSurfaceMirrorModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSurfaceMirrorModule::~qSlicerSurfaceMirrorModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSurfaceMirrorModule::helpText()const
{
  return "This is a loadable module designated to mirror a model with a surface. You can find contact information in the Acknowledgement.";
}

//-----------------------------------------------------------------------------
QString qSlicerSurfaceMirrorModule::acknowledgementText()const
{
  return "This module is done at Shanghai Jiao Tong University by Jiaxi Luo (JiaxiLuo@umich.edu) & Ruqing Ye (ruqing@umich.edu) directed by Xiaojun Chen (xiaojunchen@sjtu.edu.cn).";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSurfaceMirrorModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Jiaxi LUO, Ruqing YE, Xiaojun CHEN (SJTU)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSurfaceMirrorModule::icon()const
{
  return QIcon(":/Icons/SurfaceMirror.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSurfaceMirrorModule::categories() const
{
  return QStringList() << "Surface Models";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSurfaceMirrorModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSurfaceMirrorModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerSurfaceMirrorModule
::createWidgetRepresentation()
{
  return new qSlicerSurfaceMirrorModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSurfaceMirrorModule::createLogic()
{
  return vtkSlicerSurfaceMirrorLogic::New();
}
