/* Patch.cpp */
/* Planer patch classes' function implementation */
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


#include "Patch.h"


using std::vector;
using std::pair;


Patch::Patch()
{
  type_ = E_FACE_TYPE::DEFAULT;
  rel_pattern_coords_ = vector<EVec2d>();
  idx_array_.resize(1, 1);
  printable_state_ = FabricatablePatch();
}


Patch::Patch(const double width, const double height, const E_FACE_TYPE& face_type, bool mount)
{
  type_ = face_type;
  Initialize(width, height);
  printable_state_ = FabricatablePatch(face_type, mount);
}


Patch::Patch(const Fold& depend_fold, const E_FACE_TYPE& face_type, const EVec2d& scale)
{
  type_ = face_type;
  Initialize(depend_fold, scale);
  printable_state_ = FabricatablePatch(face_type, false);
}


Patch::Patch(const Patch& src)
{
  this->type_ = src.type_;
  this->rel_pattern_coords_ = src.rel_pattern_coords_;
  this->idx_array_ = src.idx_array_;
  this->printable_state_ = src.printable_state_;
}


void Patch::operator=(const Patch& src)
{
  this->type_ = src.type_;
  this->rel_pattern_coords_ = src.rel_pattern_coords_;
  this->idx_array_ = src.idx_array_;
  this->printable_state_ = src.printable_state_;
}


vector<EVec3d> Patch::ModelCoords(const EVec3d& depend_org, double angle) const
{
  EVec3d dir = EVec3d::Zero();
  if (type_ == E_FACE_TYPE::VTYPE)
    dir[1] = 1.0;
  else if (type_ == E_FACE_TYPE::HTYPE)
    dir[2] = 1.0;

  vector<EVec3d> model_coords(rel_pattern_coords_.size());
//#pragma omp parallel for
  for (int i = 0; i < model_coords.size(); i++)
  {
    model_coords[i] = depend_org + (rel_pattern_coords_[i][0] * EVec3d::UnitX()) + (rel_pattern_coords_[i][1] * dir);
    model_coords[i][2] -= model_coords[i][1] * cos(J_PI - angle);
    model_coords[i][1] *= sin(J_PI - angle);
  }

  return model_coords;
}


void Patch::Deform(const E_DEFORM_HANDLE& handle, double move_len, double width, double height)
{
  EVec2d max_pos(DBL_MIN, DBL_MIN);
  vector<int> move_index_x(0);
  vector<int> move_index_y(0);
  for (int i = 0; i < rel_pattern_coords_.size(); i++)
  {
    if (rel_pattern_coords_[i][0] > max_pos[0])
    {
      max_pos[0] = rel_pattern_coords_[i][0];
      move_index_x = vector<int>(0);
      move_index_x.push_back(i);
    }
    else if (rel_pattern_coords_[i][0] == max_pos[0])
    {
      move_index_x.push_back(i);
    }

    if (rel_pattern_coords_[i][1] > max_pos[1])
    {
      max_pos[1] = rel_pattern_coords_[i][1];
      move_index_y = vector<int>(0);
      move_index_y.push_back(i);
    }
    else if (rel_pattern_coords_[i][1] == max_pos[1])
    {
      move_index_y.push_back(i);
    }
  }

  for (int i = 0; i < rel_pattern_coords_.size(); i++)
  {
    if ((handle == E_DEFORM_HANDLE::LEFT) || (handle == E_DEFORM_HANDLE::RIGHT))
    {
      for (int j : move_index_x)
      {
        if (i == j)
          rel_pattern_coords_[i][0] += move_len;
        else if (rel_pattern_coords_[i][0] > max_pos[0] + move_len)
          rel_pattern_coords_[i][0] = max_pos[0] + move_len;
      }

    }
    else if (((handle == E_DEFORM_HANDLE::UP) && (type_ == E_FACE_TYPE::VTYPE)) ||
             ((handle == E_DEFORM_HANDLE::DEPTH) && (type_ == E_FACE_TYPE::HTYPE)))
    {
      for (int j : move_index_y)
      {
        if (i == j)
          rel_pattern_coords_[i][1] += move_len;
        else if (rel_pattern_coords_[i][1] > max_pos[1] + move_len)
          rel_pattern_coords_[i][1] = max_pos[1] + move_len;
      }
    }
  }
}


void Patch::Draw2D(const EVec3d& depend_org, const EVec3f& color, double angle) const
{
  vector<EVec3d> model_coords = ModelCoords(depend_org, angle);

  const EVec4f AMBI(color[0] * 0.8f, color[1] * 0.8f, color[2] * 0.8f, 1.0f);  // ä¬ã´åıîΩéÀê¨ï™
  const EVec4f DIFF(color[0], color[1], color[2], 1.0f);  // ägéUîΩéÀê¨ï™
  const EVec4f SPEC(0.2f, 0.2f, 0.2f, 1.0f);  // ãæñ îΩéÀê¨ï™
  const float  SHIN[1] = { 64.0f };  // ãæñ îΩéÀÇÃã≠ìx

  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, DIFF.data());
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, AMBI.data());
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPEC.data());
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, SHIN);

  glBegin(GL_TRIANGLES);
  JDraw::SetEVec3fToColor3fv(color);
  JDraw::SetStdVectorToGLTriangle(model_coords, idx_array_);
  glEnd();
  glFlush();
}


void Patch::Draw3D
(
  const EVec3d& depend_org,
  const EVec3f& color,
  double angle,
  const ConvertProperties& cvt_props
) const
{
  printable_state_.Draw3D(Axis(type_, angle), depend_org, color);
}


inline static double sRectW(const Rect3D& rect)
{
  return fabs(rect.p[0][0] - rect.p[1][0]);
}

inline static double sRectH(const Rect3D& rect)
{
  return fabs(rect.p[0][2] - rect.p[3][2]);
}


void Patch::FillPlane(const Rect3D& rect, const ConvertProperties& cvt_props)
{
  EVec3d x = EVec3d::UnitX();
  EVec3d y = -EVec3d::UnitY();
  EVec3d z = EVec3d::UnitZ();
  EVec3d org = (type_ == E_FACE_TYPE::VTYPE) ? rect.p[3] : rect.p[0];
  JMesh::Rectangular rectangular_param =
  { org, x, y, z, sRectW(rect), cvt_props.thickness, sRectH(rect) };
  JMesh::Mesh rectangular = JMesh::RectangularMesh(rectangular_param);

  if (!printable_state_.Mount())
  {
    int bevel_idx1 = (type_ == E_FACE_TYPE::VTYPE) ? 2 : 3;
    EVec3d pos1 = rectangular.V(bevel_idx1) +
      ((type_ == E_FACE_TYPE::VTYPE) ? EVec3d(0.0, 0.0, cvt_props.thickness) : EVec3d(0.0, 0.0, -cvt_props.thickness));
    rectangular.V(bevel_idx1, pos1);

    int bevel_idx2 = bevel_idx1 + 4;
    EVec3d pos2 = rectangular.V(bevel_idx2) +
      ((type_ == E_FACE_TYPE::VTYPE) ? EVec3d(0.0, 0.0, cvt_props.thickness) : EVec3d(0.0, 0.0, -cvt_props.thickness));
    rectangular.V(bevel_idx2, pos2);
  }

  printable_state_.FillPlane(rectangular, Axis(type_, J_PI), rect.p[0], rect.p, cvt_props);
}


void Patch::CutPrintableOAFace
(
  const JMesh::Mesh& parent_mesh,
  const Rect3D& rect,
  const Fold& depend_fold,
  const ConvertProperties& cvt_props
)
{
  printable_state_.CutPlane(parent_mesh, rect.p, depend_fold, cvt_props);
}


void Patch::GenerateInflatedPatch
(
  JMesh::Mesh& clip_plane,
  const Rect3D& rect,
  const EVec3d& depend_org,
  const ConvertProperties& cvt_props
)
{
  printable_state_.GenerateInflatedPatch(clip_plane, rect.p, Axis(type_, J_PI), depend_org, cvt_props);
}


void Patch::TrimOutlines
(
  const JMesh::Mesh& clip_plane,
  const EVec3d& depend_org,
  const ConvertProperties& cvt_props
)
{
  printable_state_.TrimOutlines(clip_plane, Axis(type_, J_PI), depend_org, cvt_props);
}


bool Patch::SegmentSweepRegion
(
  const JMesh::Mesh& higher_space,
  const JMesh::Mesh& mesh,
  const Rect3D& rect,
  const ConvertProperties& cvt_props,
  const EVec3d& max_pos,
  const vector<pair<double, double>>& concave_fold_list
)
{
  return printable_state_.SegmentSweepRegion(higher_space, mesh, rect.p, Axis(type_, J_PI / 2.0), cvt_props, max_pos, concave_fold_list);
}


JMesh::Mesh Patch::Assemble(const EVec3d& depend_org) const
{
  return printable_state_.Assemble(Axis(type_, J_PI), depend_org);
}


void Patch::ResetPatch()
{
  printable_state_.Patch(JMesh::Mesh());
}


void Patch::ResetMeshs()
{
  printable_state_.ResetMeshs();
}


bool Patch::Hit(const JUtil::Ray& ray, const EVec3d& depend_org, double& dist) const
{
  EMatXd model = JUtil::ToEMatXd(ModelCoords(depend_org, (J_PI / 2.0)));
  if (JMesh::HitPlane(ray, model, idx_array_, dist))
    return true;

  return false;
}


FaceCoordSystem Patch::Axis(const E_FACE_TYPE& face_type, double angle)
{
  FaceCoordSystem axis = { EVec3d::UnitX(), EVec3d::UnitY(), EVec3d::UnitZ() };

  if (face_type == E_FACE_TYPE::VTYPE)
  {
    axis.v_dir = EVec3d(0.0, sin(angle), cos(angle));
    axis.n_dir = EVec3d(0.0, sin(angle - (J_PI / 2.0)), cos(angle - (J_PI / 2.0)));
  }
  else if (face_type == E_FACE_TYPE::HTYPE)
  {
    axis.v_dir = EVec3d::UnitZ();
    axis.n_dir = -EVec3d::UnitY();
  }
  else
  {
    std::cout << "E_FACE_TYPE Error: OAFace::Axis(E_FACE_TYPE "
      << static_cast<int>(face_type) << ", double " << angle << ")\n";
    return { JUtil::ErrEVec3d(), JUtil::ErrEVec3d(), JUtil::ErrEVec3d() };
  }

  return axis;
}


void Patch::Initialize(double width, double height)
{
  rel_pattern_coords_ =
  {
    EVec2d(0.0, 0.0), EVec2d(width, 0.0),
    EVec2d(width, height), EVec2d(0.0, height)
  };
}


void Patch::Initialize(const Fold& depend_fold, const EVec2d& scale)
{
  double width = scale[0] * depend_fold.u;
  double height = scale[1] * ((type_ == E_FACE_TYPE::VTYPE) ? depend_fold.v : depend_fold.w);

  rel_pattern_coords_ =
  {
    EVec2d(0.0, 0.0), EVec2d(width, 0.0),
    EVec2d(width, height), EVec2d(0.0, height)
  };
}


inline static EVec3i sAlignNormal(const EVec2d& p0, const EVec2d& p1, const EVec2d& p2, const EVec3i& idx)
{
  EVec3d pos0 = EVec3d(p0[0], 0.0, p0[1]);
  EVec3d pos1 = EVec3d(p1[0], 0.0, p1[1]);
  EVec3d pos2 = EVec3d(p2[0], 0.0, p2[1]);

  EVec3d normal = JMath::ToNormalDir(pos0, pos1, pos2);
  if (normal[1] > 0.0)
  {
    EVec3i swapped = EVec3i(idx[2], idx[1], idx[0]);
    return swapped;
  }

  return idx;
}


void Patch::TrimChildren
(
  const Rect3D& rect,
  const std::vector<Rect3D>& child_rects,
  const ConvertProperties& cvt_props,
  bool mount_face
)
{
  EVec3d x = EVec3d::UnitX();
  EVec3d y = -EVec3d::UnitY();
  EVec3d z = EVec3d::UnitZ();
  EVec3d org = (type_ == E_FACE_TYPE::VTYPE) ? rect.p[3] : rect.p[0];
  JMesh::Rectangular rectangular_param =
  { org, x, y, z, sRectW(rect), cvt_props.thickness, sRectH(rect) };

  JMesh::Mesh rectangular = JMesh::RectangularMesh(rectangular_param);

  if (!mount_face)
  {
    int bevel_idx1 = (type_ == E_FACE_TYPE::VTYPE) ? 2 : 3;
    EVec3d pos1 = rectangular.V(bevel_idx1) +
      ((type_ == E_FACE_TYPE::VTYPE) ? EVec3d(0.0, 0.0, cvt_props.thickness) : EVec3d(0.0, 0.0, -cvt_props.thickness));
    rectangular.V(bevel_idx1, pos1);
    int bevel_idx2 = bevel_idx1 + 4;
    EVec3d pos2 = rectangular.V(bevel_idx2) +
      ((type_ == E_FACE_TYPE::VTYPE) ? EVec3d(0.0, 0.0, cvt_props.thickness) : EVec3d(0.0, 0.0, -cvt_props.thickness));
    rectangular.V(bevel_idx2, pos2);
  }

  for (auto c_rect : child_rects)
  {
    //EVec3d c_org = c_rect.p[3] + EVec3d(0.0, J_CSG_OFFSET, 0.0);
    EVec3d c_org = c_rect.p[3] + EVec3d(-cvt_props.gap, J_CSG_OFFSET, 0.0);
    JMesh::Rectangular c_rectangular_param =
    { c_org, x, y, z, sRectW(c_rect) + 2.0 * cvt_props.gap, (2.0 * J_CSG_OFFSET) + cvt_props.thickness, sRectH(c_rect) };
    JMesh::Mesh c_rectangular = JMesh::RectangularMesh(c_rectangular_param);
    rectangular = JCSG::Subtract(rectangular, c_rectangular);
  }

  printable_state_.Plane(rectangular, Axis(type_, J_PI), rect.p[0], rect.p, cvt_props);
  printable_state_.LoadedMesh(false);

  rel_pattern_coords_ = vector<EVec2d>(rectangular.V().rows());
  vector<EVec3i> idx_array = vector<EVec3i>(0);
  for (int i = 0; i < rectangular.F().rows(); i++)
  {
    EVec3i idx = rectangular.F(i);
    if (rectangular.V(idx[0], 1) >= 0.0 &&
        rectangular.V(idx[1], 1) >= 0.0 &&
        rectangular.V(idx[2], 1) >= 0.0)
    {
      EVec2d p0 = EVec2d(rectangular.V(idx[0], 0) - org[0], rectangular.V(idx[0], 2) - org[2]);
      EVec2d p1 = EVec2d(rectangular.V(idx[1], 0) - org[0], rectangular.V(idx[1], 2) - org[2]);
      EVec2d p2 = EVec2d(rectangular.V(idx[2], 0) - org[0], rectangular.V(idx[2], 2) - org[2]);

      rel_pattern_coords_[idx[0]] = p0;
      rel_pattern_coords_[idx[1]] = p1;
      rel_pattern_coords_[idx[2]] = p2;

      idx_array.push_back(sAlignNormal(p0, p1, p2, idx));
    }
  }

  idx_array_ = JUtil::ToEMatXi(idx_array);

  if (type_ == E_FACE_TYPE::VTYPE)
  {
    for (int i = 0; i < rel_pattern_coords_.size(); i++)
      rel_pattern_coords_[i][1] = abs(rel_pattern_coords_[i][1] - sRectH(rect));
  }

}


