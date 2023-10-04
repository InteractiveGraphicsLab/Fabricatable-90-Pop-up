/* Component.h */
/* Component class file */
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
#include "E_FoldType.h"
#include "E_FaceType.h"
#include "T_ConvertProperties.h"
#include "T_Fold.h"
#include "T_Rect3D.h"
#include "T_FaceCoordSystem.h"
#include "Patch.h"


class Component
{
public:
  Component(double width, double v_height, double h_height);
  Component(const Fold& depend_fold, const EVec2d& scale);
  Component(int id, const Fold& depend_fold, const std::vector<Fold>& folds);

  Component(const Component& src);

  ~Component() {};


public:
  int Id() const { return id_; };

  std::vector<std::vector<Component*>> Children() const { return children_; }
  std::vector<Component*> Children(const E_FOLD_TYPE& type) const { return children_[static_cast<int>(type)]; }
  std::vector<Component*> Children(int i) const { return children_[i]; }
  Component* Child(int i, int j) const { return children_[i][j]; }

  Fold  FoldEdge(const E_FOLD_TYPE& type) const { return folds_[static_cast<int>(type)]; }
  Fold  FoldEdge(int i) const { return folds_[i]; }
  Fold* FoldEdgePtr(const E_FOLD_TYPE& type) { return &(folds_[static_cast<int>(type)]); }

  Rect3D Rect(const E_FACE_TYPE& type) const { return rects_[static_cast<int>(type)]; }
  Rect3D Rect(int i) const { return rects_[i]; }

  bool Err() const { return err_; }
  void Err(bool err) { err_ = err; }


public:
  double Left() const;
  double Right() const;


public:
  bool IsExistChilren() const;


public:
  void Draw2D(double angle, const EVec3f& highlight = EVec3f::Zero()) const;
  void Draw3D(double angle, const ConvertProperties& cvt_props) const;

public:
  void Update(const Fold& depend_fold, double angle);

  void IsConnect(const Rect3D& parent_rect, const ConvertProperties& cvt_props);
  void IsMaximumSize(const Fold& root_fold);
  void IsMinimumLength(double nozzle_width);

  void TrimOAFace
  (
    const std::vector<Rect3D>& other_rects,
    const std::vector<Rect3D>& child_rects,
    const ConvertProperties& cvt_props
  );

  void AddChild(const E_FOLD_TYPE& fold_type, Component* mech);
  void DeleteChild(const E_FOLD_TYPE& fold_type, int mech_id);
  double DeformOAFace(const E_DEFORM_HANDLE& handle, double move_len, const ConvertProperties& cvt_props);

  Rect3D ComponentRect() const;


public:
  void FillPlane(const ConvertProperties& cvt_props);

  void CutPrintableOAFace
  (
    const JMesh::Mesh& mesh,
    const Fold& depend_fold,
    const ConvertProperties& cvt_props
  );

  void GenerateInflatedPatch
  (
    JMesh::Mesh& outlines,
    const ConvertProperties& cvt_props
  );

  void TrimOutlines
  (
    const E_FACE_TYPE& face_type,
    const JMesh::Mesh& outlines,
    const ConvertProperties& cvt_props
  );

  bool Segment
  (
    JMesh::Mesh& higher_space,
    const JMesh::Mesh& mesh,
    const E_FACE_TYPE& type,
    const ConvertProperties& cvt_props,
    const EVec3d& max_pos
  );

  void ResetPatch();
  void ResetMeshs();

  void AssembleOA(JMesh::Mesh& oa, const ConvertProperties& cvt_props) const;


public:
  bool HitFoldEdge
  (
    const JUtil::Ray& ray,
    double radius,
    const ConvertProperties& cvt_props,
    double& dist,
    E_FOLD_TYPE& fold_type,
    EVec3d& start,
    EVec3d& end
  ) const;

  bool Hit(const JUtil::Ray& ray, double& dist) const;

  std::vector<JMesh::Arrow> GetHandles(const JMesh::ArrowSize& arrow_size, int slice) const;


public:
  static void MechIdGenerater(int id) { s_mech_id_generator_ = id; }
  static int GenerateMechId() { return s_mech_id_generator_++; }

  static bool Trimmed() { return s_trimmed_; }
  static void Trimmed(bool trimmed) { s_trimmed_ = trimmed; }

  static EVec3f GenerateColor(int id);

private:
  void Initialize(const Fold& depend_fold, const EVec2d& scale);


private:
  void UpdateFolds(const Fold& depend_fold, double angle);
  void UpdateRects(double angle);


private:
  static int s_mech_id_generator_;
  static bool s_trimmed_;

  const int id_;
  const EVec3f color_;

  bool err_;

  std::vector<std::vector<Component*>> children_;
  std::vector<Patch> faces_;
  std::vector<Fold> folds_;
  std::vector<Rect3D> rects_;

  bool trimmed_;
  std::vector<Fold> prev_folds_;
  std::vector<Rect3D> prev_rects_;
};

