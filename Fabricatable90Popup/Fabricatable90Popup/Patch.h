/* Patch.h */
/* Planer patch class file */
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
#include "jdraw.h"
#include "jmesh.h"
#include "jcsg.h"
#include "E_FoldType.h"
#include "E_FaceType.h"
#include "E_DeformHandle.h"
#include "T_ConvertProperties.h"
#include "T_Fold.h"
#include "T_Rect3D.h"
#include "T_FaceCoordSystem.h"
#include "FabricatablePatch.h"


class Patch
{
public:
  static FaceCoordSystem Axis(const E_FACE_TYPE& face_type, double angle);


public:
  Patch();
  Patch(double width, double height, const E_FACE_TYPE& face_type, bool mount = false);
  Patch(const Fold& depend_fold, const E_FACE_TYPE& face_type, const EVec2d& scale);

  Patch(const Patch& src);
  void operator=(const Patch& src);


public:
  E_FACE_TYPE Type() const { return type_; }
  bool LoadedMesh() const { return fab_state_.LoadedMesh(); }


public:
  std::vector<EVec3d> ModelCoords(const EVec3d& depend_org, double angle) const;

  void Deform(const E_DEFORM_HANDLE& handle, double move_len, double width, double height);

  void TrimChildren
  (
    const Rect3D& rect,
    const std::vector<Rect3D>& child_rects,
    const ConvertProperties& cvt_props,
    bool mount_face
  );

public:
  void FillPlane(const Rect3D& rect, const ConvertProperties& cvt_props);

  void CutPrintableOAFace
  (
    const JMesh::Mesh& mesh,
    const Rect3D& rect,
    const Fold& depend_fold,
    const ConvertProperties& cvt_props
  );

  void GenerateInflatedPatch
  (
    const Rect3D& rect,
    const EVec3d& depend_org,
    const ConvertProperties& cvt_props
  );

  JMesh::Mesh InflatedPatch() const { return fab_state_.InflatedPatch(); }

  void TrimOutlines
  (
    const JMesh::Mesh& clip_plane,
    const EVec3d& depend_org,
    const std::vector<EVec3d>& rect,
    const ConvertProperties& cvt_props
  );

  bool SegmentSweepRegion
  (
    const JMesh::Mesh& higher_space,
    const JMesh::Mesh& mesh,
    const Rect3D& rect,
    const ConvertProperties& cvt_props,
    const EVec3d& max_pos,
    const std::vector<std::pair<double, double>>& concave_fold_list
  );

  void ResetPatch();
  void ResetMeshs();

  JMesh::Mesh Assemble(const EVec3d& depend_org) const;


public:
  void Draw2D
  (
    const EVec3d& depend_org,
    const EVec3f& color,
    double angle
  ) const;

  void Draw3D
  (
    const EVec3d& depend_org,
    const EVec3f& color,
    double angle,
    const ConvertProperties& cvt_props
  ) const;


public:
  bool Hit(const JUtil::Ray& ray, const EVec3d& depend_org, double& dist) const;


private:
  void Initialize(double width, double height);
  void Initialize(const Fold& depend_fold, const EVec2d& scale);


private:
  E_FACE_TYPE type_;
  std::vector<EVec2d> rel_pattern_coords_;
  EMatXi idx_array_;

  FabricatablePatch fab_state_;
};

