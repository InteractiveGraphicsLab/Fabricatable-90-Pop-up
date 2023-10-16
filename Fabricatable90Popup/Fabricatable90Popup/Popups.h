/* Popups.h */
/* Pop-up Mechanism classes' function implementation */
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

#include "jutil.h"
#include "jmesh.h"
#include "E_DeformStep.h"
#include "E_FoldType.h"
#include "T_ConvertProperties.h"
#include "T_PlaceProperties.h"
#include "T_DeleteProperties.h"
#include "T_DeformPropterties.h"
#include "Component.h"


class Popups
{
public:
  static Popups* GetInstance();

public:
  void Mesh(const JMesh::Mesh& mesh) { mesh_ = mesh; }

  void PlaneThick(double thickness) { cvt_props_.thickness = thickness; }
  void PreventBindGap(double gap) { cvt_props_.gap = gap; }
  void BendLineWidth(double fold_gap) { cvt_props_.fold_gap = fold_gap; }
  void LaminationPitch(double fold_thickness) { cvt_props_.fold_thickness = fold_thickness; }
  void NozzleWidth(double nozzle_width) { cvt_props_.nozzle_width = nozzle_width; }


public:
  void Draw2D() const;
  void Draw2D(const DeleteProperties& delete_props, const EVec3f& highlight) const;
  void Draw2D(const DeformProperties& deform_props, const EVec3f& highlight) const;

  void Draw3D() const;
  void DrawMesh(bool draw_mesh_plane) const;


public:
  void UpdateComponent(double angle);
  void AddComponent(const PlaceProperties& place_props);
  void DeleteComponent(const DeleteProperties& delete_props);
  double DeformComponent(const DeformProperties& deform_props, double move_len);
  void TrimComponent();


public:
  bool CheckError();
  void FillPlane();
  void ConvertComponent();
  void TrimOutlines();
  bool SegmentSpace();

  void ResetPatch();
  void ResetMeshs();

  JMesh::Mesh ExportMesh();

public:
  void LoadData(const std::vector<std::string>& data);
  std::string SaveData();


public:
  bool HitFoldEdge(const JUtil::Ray& ray, PlaceProperties& place_props) const;
  bool HitComponent(const JUtil::Ray& ray, DeleteProperties& delete_props) const;
  bool HitComponent(const JUtil::Ray& ray, DeformProperties& deform_props) const;


private:
  Popups();


private:
  void IsConnect();
  void IsMaximumSize();
  void IsMinimumLength();


private:
  Component* root_;
  JMesh::Mesh mesh_;

  double angle_;
  ConvertProperties cvt_props_;
};

