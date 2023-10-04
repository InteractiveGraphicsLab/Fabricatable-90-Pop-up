/* Fabricatable90Popup.cpp */
/* Applications' main program */
/*
* Copyright(C) 2023 Junpei Fujikawa (ma22121@shibaura-it.ac.jp)
*
* This program is free software; you can redistribute it and /or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.If not, see < https://www.gnu.org/licenses/>.
*/

#include "Fabricatable90Popup.h"


void Fabricatable90Popup::OpenLoadFileDialog()
{
  QString path = "./data/";
  QUrl url = QFileDialog::getOpenFileName(this, tr("Load File"), path, tr("*.txt"));

  if (url.isEmpty())
    return;

  std::string file_path = url.toString().toStdString();
  ui_.jqogl_widget_2d_->OrderLoadData(file_path);
}


void Fabricatable90Popup::OpenSaveFileDialog()
{
  QString path = "./data/";
  QUrl url = QFileDialog::getSaveFileName(this, tr("Save File"), path, tr("*.txt"));

  if (url.isEmpty())
    return;

  std::string file_path = url.toString().toStdString();
  ui_.jqogl_widget_3d_->OrderSaveData(file_path);
}


void Fabricatable90Popup::OpenImportFileDialog()
{
  QString path = "./data/";
  QUrl url = QFileDialog::getOpenFileName(this, tr("Load Mesh File"), path, tr("*.obj *.stl"));

  if (url.isEmpty())
    return;

  std::string file_path = url.toString().toStdString();
  ui_.jqogl_widget_2d_->OrderImportMesh(file_path);
}


void Fabricatable90Popup::OpenExportFileDialog()
{
  QString path = "./data/";
  QUrl url = QFileDialog::getSaveFileName(this, tr("Save Mesh File"), path, tr("*.stl"));

  if (url.isEmpty())
    return;

  std::string file_path = url.toString().toStdString();
  ui_.jqogl_widget_3d_->OrderExportMesh(file_path);
}


void Fabricatable90Popup::OrderConversion()
{
  ui_.jqogl_widget_3d_->OrderSetPlaneThick(ui_.plane_thick_spinbox_->value());
  ui_.jqogl_widget_3d_->OrderSetPreventBindGap(ui_.prevent_bind_gap_spinbox_->value());
  ui_.jqogl_widget_3d_->OrderSetBendLineWidth(ui_.bend_line_width_spinbox_->value());
  ui_.jqogl_widget_3d_->OrderSetLaminationPitch(ui_.lamination_pitch_spinbox_->value());
  ui_.jqogl_widget_3d_->OrderSetNozzleWidth(ui_.nozzle_width_spinbox_->value());
  ui_.jqogl_widget_3d_->OrderConvertion();
}


Fabricatable90Popup::Fabricatable90Popup(QWidget* parent)
  : QMainWindow(parent)
{
  ui_.setupUi(this);

  QSizePolicy size_policy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  size_policy.setHorizontalStretch(0);
  size_policy.setVerticalStretch(0);
  size_policy.setHeightForWidth(this->sizePolicy().hasHeightForWidth());

  SetJQOGLWidget2D(size_policy);
  SetJQOGLWidget3D(size_policy);

  QObject::connect(ui_.jqogl_widget_2d_, SIGNAL(SynchronizeOGLCamera(const OglForQt&)),
                   ui_.jqogl_widget_3d_, SLOT(SetOGLCamera(const OglForQt&)));
  QObject::connect(ui_.jqogl_widget_3d_, SIGNAL(SynchronizeOGLCamera(const OglForQt&)),
                   ui_.jqogl_widget_2d_, SLOT(SetOGLCamera(const OglForQt&)));

  QObject::connect(ui_.place_action_, SIGNAL(triggered()), ui_.jqogl_widget_2d_, SLOT(OrderPlacement()));
  QObject::connect(ui_.delete_action_, SIGNAL(triggered()), ui_.jqogl_widget_2d_, SLOT(OrderDeletion()));
  QObject::connect(ui_.deform_action_, SIGNAL(triggered()), ui_.jqogl_widget_2d_, SLOT(OrderDeformation()));

  QObject::connect(ui_.load_action_, SIGNAL(triggered()), this, SLOT(OpenLoadFileDialog()));
  QObject::connect(ui_.save_action_, SIGNAL(triggered()), this, SLOT(OpenSaveFileDialog()));
  QObject::connect(ui_.import_action_, SIGNAL(triggered()), this, SLOT(OpenImportFileDialog()));
  QObject::connect(ui_.export_action_, SIGNAL(triggered()), this, SLOT(OpenExportFileDialog()));

  QObject::connect(ui_.apply_button_, SIGNAL(clicked()), this, SLOT(OrderConversion()));

  QObject::connect(ToolCore::GetInstance(), SIGNAL(SetAngle90(int)), ui_.fold_angle_slider_, SLOT(setValue(int)));
  QObject::connect(ToolCore::GetInstance(), SIGNAL(SetAngle90(int)), ui_.fold_angle_spinbox_, SLOT(setValue(int)));

  QObject::connect(ui_.plane_thick_spinbox_, SIGNAL(valueChanged(double)), ui_.jqogl_widget_3d_, SLOT(OrderSetPlaneThick(double)));
  QObject::connect(ui_.prevent_bind_gap_spinbox_, SIGNAL(valueChanged(double)), ui_.jqogl_widget_3d_, SLOT(OrderSetPreventBindGap(double)));
  QObject::connect(ui_.bend_line_width_spinbox_, SIGNAL(valueChanged(double)), ui_.jqogl_widget_3d_, SLOT(OrderSetBendLineWidth(double)));
  QObject::connect(ui_.lamination_pitch_spinbox_, SIGNAL(valueChanged(double)), ui_.jqogl_widget_3d_, SLOT(OrderSetLaminationPitch(double)));
  QObject::connect(ui_.nozzle_width_spinbox_, SIGNAL(valueChanged(double)), ui_.jqogl_widget_3d_, SLOT(OrderSetNozzleWidth(double)));
}


Fabricatable90Popup::~Fabricatable90Popup()
{

}


void Fabricatable90Popup::SetJQOGLWidget2D(QSizePolicy& size_policy)
{
  delete ui_.jqogl_widget_2d_;

  ui_.jqogl_widget_2d_ = new JQOGLWidget(E_DRAW_MODE::MODE_2D, ui_.centralWidget);
  ui_.jqogl_widget_2d_->setObjectName(QString::fromUtf8("_jqogl_widget_2d"));
  size_policy.setHeightForWidth(ui_.jqogl_widget_2d_->sizePolicy().hasHeightForWidth());
  ui_.jqogl_widget_2d_->setSizePolicy(size_policy);
  ui_.jqogl_widget_2d_->setMinimumSize(QSize(800, 600));

  ui_.left_vertical_layout_->insertWidget(0, ui_.jqogl_widget_2d_);

  QObject::connect(ui_.fold_angle_slider_, SIGNAL(valueChanged(int)),
                   ui_.jqogl_widget_2d_, SLOT(OrderSimulation(int)));
}


void Fabricatable90Popup::SetJQOGLWidget3D(QSizePolicy& size_policy)
{
  delete ui_.jqogl_widget_3d_;

  ui_.jqogl_widget_3d_ = new JQOGLWidget(E_DRAW_MODE::MODE_3D, ui_.centralWidget);
  ui_.jqogl_widget_3d_->setObjectName(QString::fromUtf8("_jqogl_widget_3d"));
  size_policy.setHeightForWidth(ui_.jqogl_widget_3d_->sizePolicy().hasHeightForWidth());
  ui_.jqogl_widget_3d_->setSizePolicy(size_policy);
  ui_.jqogl_widget_3d_->setMinimumSize(QSize(600, 450));

  ui_.right_vertical_layout_->insertWidget(0, ui_.jqogl_widget_3d_);

  QObject::connect(ui_.fold_angle_slider_, SIGNAL(valueChanged(int)),
                   ui_.jqogl_widget_3d_, SLOT(OrderSimulation(int)));
}


