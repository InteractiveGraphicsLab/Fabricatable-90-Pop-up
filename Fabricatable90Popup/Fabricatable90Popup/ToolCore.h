/* ToolCore.h */
/* Applications' event manager class file */
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

#include "jmath.h"
#include "jutil.h"
#include "jmesh.h"
#include "E_DrawMode.h"
#include "E_EditMode.h"
#include "E_DeformStep.h"
#include "T_PlaceProperties.h"
#include "T_DeleteProperties.h"
#include "T_DeformPropterties.h"
#include "T_Fold.h"
#include "Popups.h"


class ToolCore
  : public QObject
{
  Q_OBJECT


signals:
  void SetAngle90(int angle);


public:
  static ToolCore* GetInstance();


public:
  void PressLeftMouseButton(OglForQt* ogl, const EVec2i& p);
  void PressRightMouseButton(OglForQt* ogl, const EVec2i& p);
  void PressMiddleMouseButton(OglForQt* ogl, const EVec2i& p);

  void ReleaseLeftMouseButton(OglForQt* ogl, const EVec2i& p);
  void ReleaseRightMouseButton(OglForQt* ogl, const EVec2i& p);
  void ReleaseMiddleMouseButton(OglForQt* ogl, const EVec2i& p);

  void MoveMouse(OglForQt* ogl, const EVec2i& p);

  void DrawScene(OglForQt* ogl, const E_DRAW_MODE& draw_mode, bool draw_mesh_plane);


public:
  void EditMode(OglForQt* ogl, const E_EDIT_MODE& edit_mode);

  void Simulate(OglForQt* ogl, int angle);
  void Simulate(OglForQt* ogl, double angle);
  void InitializeDeformMode(OglForQt* ogl);
  void ImportMesh(OglForQt* ogl, const std::string& file_path);
  void ExportMesh(OglForQt* ogl, const std::string& file_path);
  void LoadData(OglForQt* ogl, const std::string& file_path);
  void SaveData(OglForQt* ogl, const std::string& file_path);
  void Convert(OglForQt* ogl);


public:
  void SetPlaneThick(OglForQt* ogl, double thickness);
  void SetPreventBindGap(OglForQt* ogl, double gap);
  void SetBendLineWidth(OglForQt* ogl, double fold_gap);
  void SetLaminationPitch(OglForQt* ogl, double fold_thickness);
  void SetNozzleWidth(OglForQt* ogl, double nozzle_width);

private:
  ToolCore();
  void InitializeProps();


private:
  void ProgressDeformStep();

  bool HitDeformHandle(const JUtil::Ray& ray);
  void UpdateDeformHandle(double move_len);

  void DrawPlaceFoldEdge() const;
  void DrawDeformHandle() const;


private:
  bool Placable() const;
  bool Deletable() const;
  bool Deformable() const;

  bool SelectedMech() const;
  bool MovableHandle() const;



private:
  JUtil::Mouse mouse_;
  JUtil::Ray prev_ray_;

  E_EDIT_MODE edit_mode_;

  PlaceProperties place_props_;
  DeformProperties deform_props_;
  DeleteProperties delete_props_;

  bool imported_mesh_;
  bool converted_;
};

