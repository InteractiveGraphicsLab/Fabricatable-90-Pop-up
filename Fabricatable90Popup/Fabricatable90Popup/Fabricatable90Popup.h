/* Fabricatable90Popup.h */
/* Applications' main class file */
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

#include "ui_Fabricatable90Popup.h"
#include "ToolCore.h"

class Fabricatable90Popup
  : public QMainWindow
{
  Q_OBJECT

public slots:
  void OpenLoadFileDialog();
  void OpenSaveFileDialog();
  void OpenImportFileDialog();
  void OpenExportFileDialog();


public slots:
  void OrderConversion();


public:
  Fabricatable90Popup(QWidget* parent = nullptr);
  ~Fabricatable90Popup();


private:
  void SetJQOGLWidget2D(QSizePolicy& size_policy);
  void SetJQOGLWidget3D(QSizePolicy& size_policy);


private:
  Ui::Fabricatable90PopupClass ui_;
};