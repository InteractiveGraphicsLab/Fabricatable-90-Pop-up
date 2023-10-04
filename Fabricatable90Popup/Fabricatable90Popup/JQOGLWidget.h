/* JQOGLWidget.h */
/* Draw OpenGL widget classes' function implementation */
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


#pragma once

#include "E_DrawMode.h"
#include "ToolCore.h"


class JQOGLWidget :
  public QOpenGLWidget
{
  Q_OBJECT


signals:
  void SynchronizeOGLCamera(const OglForQt& ogl);


public slots:
  void SetOGLCamera(const OglForQt& ogl);


public slots:
  void OrderSimulation(int angle);
  void OrderPlacement();
  void OrderDeletion();
  void OrderDeformation();
  void OrderConvertion();


public slots:
  void OrderSetPlaneThick(double thickness);
  void OrderSetPreventBindGap(double gap);
  void OrderSetBendLineWidth(double fold_gap);
  void OrderSetLaminationPitch(double fold_thickness);
  void OrderSetNozzleWidth(double nozzle_width);


public:
  explicit JQOGLWidget(QWidget* parent = Q_NULLPTR);
  explicit JQOGLWidget(const E_DRAW_MODE& draw_mode, QWidget* parent = Q_NULLPTR);


public:
  void OrderLoadData(const std::string& file_path);
  void OrderSaveData(const std::string& file_path);
  void OrderImportMesh(const std::string& file_path);
  void OrderExportMesh(const std::string& file_path);


protected:
  void initializeGL() override;
  void resizeGL(int width, int height) override;
  void paintGL() override;

  void mousePressEvent(QMouseEvent* mouse_event) override;
  void mouseReleaseEvent(QMouseEvent* mouse_event) override;
  void mouseMoveEvent(QMouseEvent* mouse_event) override;

  void keyPressEvent(QKeyEvent* key_event) override;
  void keyReleaseEvent(QKeyEvent* key_event) override;


private:
  const E_DRAW_MODE draw_mode_;
  bool draw_mesh_plane_;
  OglForQt ogl_;
};

