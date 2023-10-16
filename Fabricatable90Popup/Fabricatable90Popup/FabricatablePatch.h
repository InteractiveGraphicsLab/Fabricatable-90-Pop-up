/* FabricatablePatch.h */
/* Fabricatable patch class file */
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
#include "jdraw.h"
#include "jcsg.h"
#include "E_FaceType.h"
#include "T_ConvertProperties.h"
#include "T_Fold.h"
#include "T_FaceCoordSystem.h"
#include "T_Rect3D.h"


class FabricatablePatch
{
public:
  FabricatablePatch();
  FabricatablePatch(const E_FACE_TYPE& type, bool mount_face);

  FabricatablePatch(const FabricatablePatch& src);

  ~FabricatablePatch() {}


public:
  bool Mount() const { return mount_face_; }
  bool LoadedMesh() const { return loaded_mesh_; }
  void LoadedMesh(bool loaded_mesh) { loaded_mesh_ = loaded_mesh; }

  JMesh::Mesh InflatedPatch() const { return inflated_plane_; }


public:
  void Plane
  (
    const JMesh::Mesh& plane,
    const FaceCoordSystem& face_system,
    const EVec3d& depend_org,
    const std::vector<EVec3d>& rect,
    const ConvertProperties& cvt_props
  );

  void FillPlane
  (
    const JMesh::Mesh& plane,
    const FaceCoordSystem& face_system,
    const EVec3d& depend_org,
    const std::vector<EVec3d>& rect,
    const ConvertProperties& cvt_props
  );


public:
  void Patch(const JMesh::Mesh& patch);
  void ResetMeshs();


public:
  void Draw3D(const FaceCoordSystem& face_system, const EVec3d& depend_org, const EVec3f& color) const;


public:
  void CutPatch
  (
    const JMesh::Mesh& mesh,
    const std::vector<EVec3d>& rect,
    const Fold& depend_fold,
    const ConvertProperties& cvt_props
  );

  void GenerateInflatedPatch
  (
    const std::vector<EVec3d> rect,
    const FaceCoordSystem& face_system,
    const EVec3d& depend_org,
    const ConvertProperties& cvt_props
  );

  void TrimOutlines
  (
    const JMesh::Mesh& outlines,
    const FaceCoordSystem& face_system,
    const EVec3d& depend_org,
    const std::vector<EVec3d>& rect,
    const ConvertProperties& cvt_props
  );


  bool SegmentSweepRegion
  (
    const JMesh::Mesh& higher_space,
    const JMesh::Mesh& mesh,
    const std::vector<EVec3d>& rect,
    const FaceCoordSystem& face_system,
    const ConvertProperties& cvt_props,
    const EVec3d& max_pos,
    const std::vector<std::pair<double, double>>& base_fold_list
  );


public:
  JMesh::Mesh Assemble(const FaceCoordSystem& face_system, const EVec3d& depend_org) const;


private:
  JMesh::Mesh InfatePatch
  (
    const std::vector<EVec3d> rect,
    const FaceCoordSystem& face_system,
    const EVec3d& depend_org,
    const ConvertProperties& cvt_props
  );


private:
  JMesh::Mesh GenerateConvex(const std::vector<EVec3d>& rect, const ConvertProperties& cvt_props) const;
  JMesh::Mesh GenerateConvexCutBox(const std::vector<EVec3d>& rect, const ConvertProperties& cvt_props) const;

  JMesh::Mesh GenerateConcaveCutBox90(const std::vector<EVec3d>& rect, const Fold& depend_fold, const ConvertProperties& cvt_props) const;
  JMesh::Mesh GenerateConcaveCutBox180(const std::vector<EVec3d>& rect, const ConvertProperties& cvt_props) const;


private:
  E_FACE_TYPE type_;

  bool mount_face_;
  bool loaded_mesh_;

  JMesh::Mesh rel_plane_;
  JMesh::Mesh rel_original_plane_;

  JMesh::Mesh inflated_plane_;
  cv::Mat projected_img_;

  JMesh::Mesh rel_outside_;

  JMesh::Mesh rel_patch_;
};

