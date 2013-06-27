/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerSurfaceMirrorFooBarWidget.h"
#include "ui_qSlicerSurfaceMirrorFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SurfaceMirror
class qSlicerSurfaceMirrorFooBarWidgetPrivate
  : public Ui_qSlicerSurfaceMirrorFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSurfaceMirrorFooBarWidget);
protected:
  qSlicerSurfaceMirrorFooBarWidget* const q_ptr;

public:
  qSlicerSurfaceMirrorFooBarWidgetPrivate(
    qSlicerSurfaceMirrorFooBarWidget& object);
  virtual void setupUi(qSlicerSurfaceMirrorFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSurfaceMirrorFooBarWidgetPrivate
::qSlicerSurfaceMirrorFooBarWidgetPrivate(
  qSlicerSurfaceMirrorFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSurfaceMirrorFooBarWidgetPrivate
::setupUi(qSlicerSurfaceMirrorFooBarWidget* widget)
{
  this->Ui_qSlicerSurfaceMirrorFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSurfaceMirrorFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSurfaceMirrorFooBarWidget
::qSlicerSurfaceMirrorFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSurfaceMirrorFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSurfaceMirrorFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSurfaceMirrorFooBarWidget
::~qSlicerSurfaceMirrorFooBarWidget()
{
}
