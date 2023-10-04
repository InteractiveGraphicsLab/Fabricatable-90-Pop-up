/* JQOGLWidget.cpp */
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

#include "JQOGLWidget.h"


JQOGLWidget::JQOGLWidget(QWidget *parent)
	: draw_mode_(E_DRAW_MODE::DEFAULT), ogl_(this), QOpenGLWidget(parent)
{

}


JQOGLWidget::JQOGLWidget(const E_DRAW_MODE &draw_mode, QWidget *parent)
	: draw_mode_(draw_mode), ogl_(this), QOpenGLWidget(parent)
{
	if (draw_mode == E_DRAW_MODE::MODE_2D)
		this->setFocusPolicy(Qt::StrongFocus);

	draw_mesh_plane_ = true;
}


void JQOGLWidget::SetOGLCamera(const OglForQt &ogl)
{
	ogl_.SetCam(ogl.GetCamPos(), ogl.GetCamCnt(), ogl.GetCamUp());
	ogl_.RedrawWindow();
}


void JQOGLWidget::OrderSimulation(int angle)
{
	ToolCore::GetInstance()->Simulate(&ogl_, angle);
}


void JQOGLWidget::OrderPlacement()
{
	ToolCore::GetInstance()->EditMode(&ogl_, E_EDIT_MODE::PLACEMENT);
}


void JQOGLWidget::OrderDeletion()
{
	ToolCore::GetInstance()->EditMode(&ogl_, E_EDIT_MODE::DELETION);
}


void JQOGLWidget::OrderDeformation()
{
	ToolCore::GetInstance()->EditMode(&ogl_, E_EDIT_MODE::DEFORMATION);
	ToolCore::GetInstance()->InitializeDeformMode(&ogl_);
}


void JQOGLWidget::OrderLoadData(const std::string& file_path)
{
  ToolCore::GetInstance()->EditMode(&ogl_, E_EDIT_MODE::LOAD_DATA);
  ToolCore::GetInstance()->LoadData(&ogl_, file_path);
}


void JQOGLWidget::OrderSaveData(const std::string& file_path)
{
  ToolCore::GetInstance()->SaveData(&ogl_, file_path);
}


void JQOGLWidget::OrderImportMesh(const std::string &file_path)
{
	ToolCore::GetInstance()->EditMode(&ogl_, E_EDIT_MODE::IMPORT_MESH);
	ToolCore::GetInstance()->ImportMesh(&ogl_, file_path);
}


void JQOGLWidget::OrderExportMesh(const std::string &file_path)
{
	ToolCore::GetInstance()->ExportMesh(&ogl_, file_path);
}


void JQOGLWidget::OrderConvertion()
{
  ToolCore::GetInstance()->EditMode(&ogl_, E_EDIT_MODE::CONVERSION);
  ToolCore::GetInstance()->Convert(&ogl_);
}


void JQOGLWidget::OrderSetPlaneThick(double thickness)
{
	ToolCore::GetInstance()->SetPlaneThick(&ogl_, thickness);
}


void JQOGLWidget::OrderSetPreventBindGap(double gap)
{
	ToolCore::GetInstance()->SetPreventBindGap(&ogl_, gap);
}


void JQOGLWidget::OrderSetBendLineWidth(double fold_gap)
{
	ToolCore::GetInstance()->SetBendLineWidth(&ogl_, fold_gap);
}


void JQOGLWidget::OrderSetLaminationPitch(double fold_thickness)
{
	ToolCore::GetInstance()->SetLaminationPitch(&ogl_, fold_thickness);
}


void JQOGLWidget::OrderSetNozzleWidth(double nozzle_width)
{
	ToolCore::GetInstance()->SetNozzleWidth(&ogl_, nozzle_width);
}


void JQOGLWidget::initializeGL()
{
	ogl_.SetBgColor(1.0f, 1.0f, 1.0f, 0.5f);

  ogl_.SetCam(EVec3f(-2.099365, 7.457740, 12.407666),
              EVec3f(4.174096, 2.911306, 2.771953),
              EVec3f(-0.461784, 1.125189, 0.645092));

	this->setMouseTracking(true);
}


void JQOGLWidget::resizeGL(int width, int height)
{
	this->repaint();
}


void JQOGLWidget::paintGL()
{
	const int WIDTH = this->width();
	const int HEIGHT = this->height();

	ogl_.OnDrawBegin(WIDTH, HEIGHT);
	ToolCore::GetInstance()->DrawScene(&ogl_, draw_mode_, draw_mesh_plane_);
	ogl_.OnDrawEnd();
}


void JQOGLWidget::mousePressEvent(QMouseEvent* mouse_event)
{
	EVec2i p = EVec2i(mouse_event->x(), mouse_event->y());
	switch (mouse_event->button())
	{
		case Qt::LeftButton:
			ToolCore::GetInstance()->PressLeftMouseButton(&ogl_, p);
      ogl_.RedrawWindow();
      break;
		case Qt::RightButton:
			ToolCore::GetInstance()->PressRightMouseButton(&ogl_, p);
      ogl_.RedrawWindow();
			break;
		case Qt::MiddleButton:
			ToolCore::GetInstance()->PressMiddleMouseButton(&ogl_, p);
      ogl_.RedrawWindow();
			break;
		default:
			break;
	}
}


void JQOGLWidget::mouseReleaseEvent(QMouseEvent *mouse_event)
{
	EVec2i p = EVec2i(mouse_event->x(), mouse_event->y());
	switch (mouse_event->button())
	{
		case Qt::LeftButton:
			ToolCore::GetInstance()->ReleaseLeftMouseButton(&ogl_, p);
      ogl_.RedrawWindow();
      break;
		case Qt::RightButton:
			ToolCore::GetInstance()->ReleaseRightMouseButton(&ogl_, p);
      ogl_.RedrawWindow();
      break;
		case Qt::MiddleButton:
			ToolCore::GetInstance()->ReleaseMiddleMouseButton(&ogl_, p);
      ogl_.RedrawWindow();
      break;
		default:
			break;
	}
}


void JQOGLWidget::mouseMoveEvent(QMouseEvent *mouse_event)
{
	EVec2i p = EVec2i(mouse_event->x(), mouse_event->y());
	ToolCore::GetInstance()->MoveMouse(&ogl_, p);
	emit SynchronizeOGLCamera(ogl_);
  ogl_.RedrawWindow();
}


void JQOGLWidget::keyPressEvent(QKeyEvent *key_event)
{
	if ((draw_mode_ == E_DRAW_MODE::MODE_2D) && (key_event->key() == Qt::Key_Shift))
		draw_mesh_plane_ = false;

	ogl_.RedrawWindow();
}


void JQOGLWidget::keyReleaseEvent(QKeyEvent *key_event)
{
	if ((draw_mode_ == E_DRAW_MODE::MODE_2D) && (key_event->key() == Qt::Key_Shift))
		draw_mesh_plane_ = true;

	ogl_.RedrawWindow();
}



