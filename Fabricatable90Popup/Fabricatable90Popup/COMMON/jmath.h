/* jmesh.h */
/* Mathmatical calculate functions */
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

#include "tmath.h"
#include <cmath>
#include <string>

#ifndef J_PI
#define J_PI 3.1415926
#endif  // J_PI


typedef Eigen::MatrixXd EMatXd;
typedef Eigen::MatrixXi EMatXi;

namespace JMath
{
  inline EVec3d RotateAroundAxis
  (
    const EVec3d& pos,
    const EVec3d& axis,
    const double& angle
  )
  {
    EMat3d axis_angle = EMat3d::Identity();
    axis_angle = Eigen::AngleAxisd(angle, axis);
    EVec3d rotated_pos = axis_angle * pos;

    return rotated_pos;
  }


  inline std::pair<EVec3f, EVec3f> CalcInsecPlane
  (
    const EVec3f& c1,
    const double& r1,
    const EVec3f& c2,
    const double& r2
  )
  {
    if ((c1 - c2).norm() > 0.0f)
    {
      double d = (c1 - c2).norm();	// d
      double cos_theta = ((d * d) + (r1 * r1) - (r2 * r2)) / (2.0 * d * r1);

      EVec3f J = c1 + (r1 * cos_theta) * (c2 - c1).normalized();
      EVec3f normal = (c2 - c1).normalized();

      return { J, normal };
    }

    return { c1, EVec3f::Zero() };
  }


  inline EVec3f CalcPlanesIntersectionLine
  (
    const EVec3f& C1, const float& r1,
    const EVec3f& C2, const float& r2,
    const EVec3f& C3, const float& r3
  )
  {
    std::pair<EVec3f, EVec3f> P = CalcInsecPlane(C1, r1, C2, r2);
    std::pair<EVec3f, EVec3f> Q = CalcInsecPlane(C1, r1, C3, r3);

    float offset_P = -(P.first[0] * P.second[0]) - (P.first[1] * P.second[1]) - (P.first[2] * P.second[2]);
    float offset_Q = -(Q.first[0] * Q.second[0]) - (Q.first[1] * Q.second[1]) - (Q.first[2] * Q.second[2]);

    double Px = (double)P.second[0];
    double Py = (double)P.second[1];
    double Pz = (double)P.second[2];
    double Pd = (double)offset_P;

    double Qx = (double)Q.second[0];
    double Qy = (double)Q.second[1];
    double Qz = (double)Q.second[2];
    double Qd = (double)offset_Q;

    double p = (Pz * Qd) - (Pd * Qz);
    double q = (Py * Qd) - (Pd * Qy);
    double r = (Py * Qz) - (Pz * Qy);
    double s = (Px * Qd) - (Pd * Qx);
    double t = (Px * Qz) - (Pz * Qx);
    double u = (Px * Qy) - (Py * Qx);

    EVec3d V = EVec3d(r, -t, u).normalized();
    EVec3d B = EVec3d(0.0, 0.0, 0.0);

    if (r != 0.0) { B = EVec3d(0.0, (p / r), -(q / r)); }
    else if (t != 0.0) { B = EVec3d((p / t), 0.0, -(s / t)); }
    else if (u != 0.0) { B = EVec3d((q / u), -(s / u), 0.0); }
    else { return EVec3f(0.0f, 0.0f, 0.0f); }

    EVec3d C = C1.cast<double>();
    double radius = r1;
    if ((r1 >= r2) && (r1 >= r3)) { C = C1.cast<double>(); radius = r1; }
    else if ((r2 >= r1) && (r2 >= r3)) { C = C2.cast<double>(); radius = r2; }
    else if ((r3 >= r1) && (r3 >= r2)) { C = C3.cast<double>(); radius = r3; }

    double a = V.squaredNorm();
    double b = 2.0 * ((B - C).dot(V));
    double c = (B - C).squaredNorm() - pow(radius, 2);

    double t1 = (-b + sqrt(pow(b, 2) - (4.0 * a * c))) / (2.0 * a);
    double t2 = (-b - sqrt(pow(b, 2) - (4.0 * a * c))) / (2.0 * a);

    if ((t1 < 0.0) && (t2 < 0.0)) { return EVec3f(0.0f, 0.0f, 0.0f); }

    EVec3f pt, dir;

    if ((V[1] >= 0.0))
    {
      if (t1 >= t2)
      {
        pt = B.cast<float>();
        dir = (B + (t1 * V)).cast<float>();
      }
      else if (t2 > t1)
      {
        pt = B.cast<float>();
        dir = (B + (t2 * V)).cast<float>();
      }
      else
      {
        pt = EVec3f(0.0f, 0.0f, 0.0f);
        dir = EVec3f(0.0f, 0.0f, 0.0f);
      }
    }
    else
    {
      if (t1 <= t2)
      {
        pt = B.cast<float>();
        dir = (B + (t1 * V)).cast<float>();
      }
      else if (t2 < t1)
      {
        pt = B.cast<float>();
        dir = (B + (t2 * V)).cast<float>();
      }
      else
      {
        pt = EVec3f(0.0f, 0.0f, 0.0f);
        dir = EVec3f(0.0f, 0.0f, 0.0f);
      }
    }

    return dir;
  }


  inline std::pair<EVec3f, float> CalcProjectedHingePt
  (
    const EVec3f& c1,
    const float   r1,
    const EVec3f& c2,
    const float   r2,
    const float   gap
  )
  {
    float d = 0.0f;

    if ((c1 - c2).norm() > gap)
      d = (c1 - c2).norm() - gap;

    float a = 0.0f;

    if (d > 0.0f)
      a = ((r1 * r1) + (d * d) - (r2 * r2)) / (2.0f * d);

    float m = 0.0f;
    EVec3f M = c1;

    if (r1 > a)
      m = sqrtf((r1 * r1) - (a * a));

    if (a > 0.0f)
      M = c1 + (a * (c2 - c1).normalized());

    return std::make_pair(M, m);
  }


  inline EVec3f CalcCenterPt(const EVec3f& pt1, const EVec3f& pt2)
  {
    EVec3f dir_pt1_to_pt2 = (pt2 - pt1).normalized();
    double distance = (pt2 - pt1).norm();
    EVec3f center_pt = pt1 + dir_pt1_to_pt2 * (distance / 2);

    return center_pt;
  }


  inline double CalcSubtendAngle(const EVec3d& v1, const EVec3d& v2)
  {
    double inner_product = (v1.dot(v2));
    double v1_norm = v1.norm();
    double v2_norm = v2.norm();

    double cos_theta = inner_product / (v1_norm * v2_norm);
    double theta = acosf(cos_theta);

    return theta;
  }


  inline EVec3f MovePtAlongDir
  (
    const EVec3f& origin,
    const EVec3f& pt,
    const EVec3f& dir
  )
  {
    EVec3f move = pt - origin;
    float move_dist = dir.dot(move);
    EVec3f moved_pt = dir * move_dist;

    return origin + moved_pt;
  }


  inline EVec2f Project3DPtTo2DPlane
  (
    const EVec3f& pt,
    const EVec3f& plane_origin,
    const EVec3f& plane_dir1,
    const EVec3f& plane_dir2
  )
  {
    EVec3f prev_project = pt - plane_origin;
    EVec2f projected_pt = EVec2f(prev_project.dot(plane_dir1), prev_project.dot(plane_dir1));

    return projected_pt;
  }


  inline float CalcPythagoreanHypo(float height, float bottom)
  {
    return sqrtf(powf(height, 2.0f) + powf(bottom, 2.0f));
  }

  inline float CalcPythagoreanAnother(float hypo, float bottom)
  {
    return sqrtf(powf(hypo, 2.0f) - powf(bottom, 2.0f));
  }


  inline float CalcDirectionMoveDist
  (
    const EVec3f& dir,
    const EVec3f& start_pt,
    const EVec3f& moved_pt
  )
  {
    EVec3f moved = moved_pt - start_pt;
    return moved.dot(dir);
  }


  inline EVec3f CalcInternalDivisionPt
  (
    float m,
    float n,
    const EVec3f& org,
    const EVec3f& p1,
    const EVec3f& p2
  )
  {
    EVec3f a = p1 - org;
    EVec3f b = p2 - org;

    return (n * a + m * b) / (m + n);
  }


  inline EVec3f CalcBarycenter(const std::vector<EVec3f>& pos_list)
  {
    EVec3f ret = EVec3f::Zero();
    for (EVec3f pos : pos_list)
      ret += pos;

    return ret / (float)(pos_list.size());
  }


  inline float Rad(int angle)
  {
    return ((float)angle / 180.0f) * J_PI;
  }

  inline int Deg(float radian)
  {
    return (int)((180.0f / J_PI) * radian);
  }


  inline EVec3f ToNormalDir(const EVec3f& org, const EVec3f& pt0, const EVec3f& pt1)
  {
    EVec3f v0 = pt0 - org;
    EVec3f v1 = pt1 - org;

    return v0.cross(v1).normalized();
  }

  inline EVec3d ToNormalDir(const EVec3d& org, const EVec3d& pt0, const EVec3d& pt1)
  {
    EVec3d v0 = pt0 - org;
    EVec3d v1 = pt1 - org;

    return v0.cross(v1).normalized();
  }


  inline EVec3f SolveEquation
  (
    const EVec3f& P1, float& s, const EVec3f& d1,
    const EVec3f& P2, float& t, const EVec3f& d2
  )
  {
    float  d = d1.dot(d2);
    EVec3f C = P1 - P2;
    float  b = C.dot(d2);
    EVec3f A = d1 - (d * d2);

    if (A.norm() == 0.0f)
    {
      s = t = 0.0f;
      return EVec3f::Zero();
    }

    s = (b * d2 - C).dot(A) / A.squaredNorm();
    t = b + s * d;

    EVec3f Q = P1 + s * d1;

    return Q;
  }


  inline EMatXd RotateMat(const EMatXd& mat, const EVec3d& axis)
  {
    EMatXd rotated_mat;
    rotated_mat.resize(mat.rows(), mat.cols());

    EVec3d x_axis = axis.cross(EVec3d(1.0, 0.0, 0.0));
    EVec3d y_axis = axis.cross(EVec3d(0.0, 1.0, 0.0));
    EVec3d z_axis = axis.cross(EVec3d(0.0, 0.0, 1.0));

    double x_radian = CalcSubtendAngle(axis, EVec3d(1.0, 0.0, 0.0));
    double y_radian = CalcSubtendAngle(axis, EVec3d(0.0, 1.0, 0.0));
    double z_radian = CalcSubtendAngle(axis, EVec3d(0.0, 0.0, 1.0));

    for (int i = 0; i < mat.rows(); ++i)
    {
      EVec3d pos = EVec3d(mat(i, 0), mat(i, 1), mat(i, 2));

      pos = RotateAroundAxis(pos, x_axis, x_radian);
      pos = RotateAroundAxis(pos, y_axis, y_radian);
      pos = RotateAroundAxis(pos, z_axis, z_radian);

      rotated_mat(i, 0) = pos[0];
      rotated_mat(i, 1) = pos[1];
      rotated_mat(i, 2) = pos[2];
    }

    return rotated_mat;
  }


  inline EMatXd RotateMatAroundAxis(const EMatXd& pts, const EVec3d& axis, const double& radian)
  {
    EMatXd rotated_mat;
    rotated_mat.resize(pts.rows(), pts.cols());

    for (int i = 0; i < pts.rows(); ++i)
    {
      EVec3d pos = EVec3d(pts(i, 0), pts(i, 1), pts(i, 2));
      EVec3d rotated_pos = RotateAroundAxis(pos, axis, radian);

      rotated_mat(i, 0) = rotated_pos[0];
      rotated_mat(i, 1) = rotated_pos[1];
      rotated_mat(i, 2) = rotated_pos[2];
    }

    return rotated_mat;
  }


  inline EMatXd MoveMat(const EMatXd& pts, const EVec3d& move)
  {
    EMat4d move_matrix;
    move_matrix << 1.0, 0.0, 0.0, move[0],
      0.0, 1.0, 0.0, move[1],
      0.0, 0.0, 1.0, move[2],
      0.0, 0.0, 0.0, 1.0;

    EMatXd moved_pts;
    moved_pts.resize(pts.rows(), pts.cols());

    for (int i = 0; i < pts.rows(); ++i)
    {
      EVec4d pt = EVec4d(pts(i, 0), pts(i, 1), pts(i, 2), 1.0);
      EVec4d moved_pt = move_matrix * pt;

      moved_pts(i, 0) = moved_pt[0];
      moved_pts(i, 1) = moved_pt[1];
      moved_pts(i, 2) = moved_pt[2];
    }

    return moved_pts;
  }


  inline EMatXd ScaleMat(const EMatXd& pts, const Eigen::Vector3d& scale)
  {
    EMat4d scale_matrix;
    scale_matrix << scale[0], 0.0, 0.0, 0.0,
      0.0, scale[1], 0.0, 0.0,
      0.0, 0.0, scale[2], 0.0,
      0.0, 0.0, 0.0, 1.0;

    EMatXd scaled_pts;
    scaled_pts.resize(pts.rows(), pts.cols());

    for (int i = 0; i < pts.rows(); i++)
    {
      EVec4d pt = EVec4d(pts(i, 0), pts(i, 1), pts(i, 2), 1.0);
      EVec4d scaled_pt = scale_matrix * pt;

      scaled_pts(i, 0) = scaled_pt[0];
      scaled_pts(i, 1) = scaled_pt[1];
      scaled_pts(i, 2) = scaled_pt[2];
    }

    return scaled_pts;
  }


  inline EMatXd TransformCoodinateToRelative
  (
    const EVec3f& org,
    const EVec3f& v0,
    const EVec3f& v1,
    const EVec3f& v2,
    const EMatXd& V
  )
  {
    EMatXd V_return;
    V_return.resize(V.rows(), V.cols());

    EVec3d O = org.cast<double>();
    EVec3d A = v0.cast<double>();
    EVec3d B = v1.cast<double>();
    EVec3d C = v2.cast<double>();

    for (int i = 0; i < V.rows(); ++i)
    {
      EVec3d pt = EVec3d(V(i, 0), V(i, 1), V(i, 2));
      EVec3d rel_pt = pt - O;
      V_return(i, 0) = rel_pt.dot(A);
      V_return(i, 1) = rel_pt.dot(B);
      V_return(i, 2) = rel_pt.dot(C);
    }

    return V_return;
  }


  inline EMatXd TransformCoodinateToStandard
  (
    const EVec3f& org,
    const EVec3f& v0,
    const EVec3f& v1,
    const EVec3f& v2,
    const EMatXd& V
  )
  {
    EMatXd V_return;
    V_return.resize(V.rows(), V.cols());

    EVec3d O = org.cast<double>();
    EVec3d A = v0.cast<double>();
    EVec3d B = v1.cast<double>();
    EVec3d C = v2.cast<double>();

    for (int i = 0; i < V.rows(); ++i)
    {
      EVec3d rel_pt = EVec3d(V(i, 0), V(i, 1), V(i, 2));
      EVec3d pt = O + (rel_pt[0] * A) + (rel_pt[1] * B) + (rel_pt[2] * C);
      V_return(i, 0) = pt[0];
      V_return(i, 1) = pt[1];
      V_return(i, 2) = pt[2];
    }

    return V_return;
  }


  inline double Dist2D(const EVec2d& p1, const EVec2d& p2)
  {
    return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
                (p1[1] - p2[1]) * (p1[1] - p2[1]));
  }


  inline EVec2d Row2D(const EMatXd& V, int idx)
  {
    return EVec2d(V(idx, 0), V(idx, 1));
  }


  inline int SearchFarPosIdx2D(const EMatXd& V, const EVec2d& org)
  {
    int far_idx = -1;
    double far_dist = -1.0;
    for (int i = 0; i < V.rows(); i++)
    {
      double dist = Dist2D(org, EVec2d(V(i, 0), V(i, 1)));
      if (far_dist < dist)
      {
        far_idx = i;
        far_dist = dist;
      }
    }

    return far_idx;
  }


  inline int SearchFarPosIdxInRange2D(const EMatXd& V, const std::vector<int> range, const EVec2d& org)
  {
    int far_idx = -1;
    double far_dist = -1.0;
    for (int i = 0; i < range.size(); i++)
    {
      double dist = Dist2D(org, EVec2d(V(range[i], 0), V(range[i], 1)));
      if (far_dist < dist)
      {
        far_idx = range[i];
        far_dist = dist;
      }
    }

    return far_idx;
  }


  inline double CalcSubtendAngle2D(const EVec2d& v1, const EVec2d& v2)
  {
    double inner_product = (v1.dot(v2));
    double v1_norm = v1.norm();
    double v2_norm = v2.norm();

    double cos_theta = inner_product / (v1_norm * v2_norm);
    double theta = acosf(cos_theta);

    return theta;
  }


  inline double Calc3VertexAngle2D(const EVec2d& org, const EVec2d& v1, const EVec2d& v2)
  {
    return CalcSubtendAngle2D(v1 - org, v2 - org);
  }


  inline bool InsideTriangle2D(const EVec2d& p, const EVec2d& a, const EVec2d& b, const EVec2d& c)
  {
    EVec2d ab = b - a;  EVec2d bp = p - b;
    EVec2d bc = c - b;  EVec2d cp = p - c;
    EVec2d ca = a - c;  EVec2d ap = p - a;

    double c1 = ab[0] * bp[1] - ab[1] * bp[0];
    double c2 = bc[0] * cp[1] - bc[1] * cp[0];
    double c3 = ca[0] * ap[1] - ca[1] * ap[0];

    if (((c1 < 0.0) && (c2 < 0.0) && (c3 < 0.0)) ||
        ((c1 > 0.0) && (c2 > 0.0) && (c3 > 0.0)))
    {
      return true;
    }

    return false;
  }


  inline bool InsideTriangle2D(const EMatXd& V, int i1, int i2, int i3)
  {
    for (int i = 0; i < V.rows(); i++)
    {
      if (i != i1 && i != i2 && i != i3)
      {
        if (InsideTriangle2D(Row2D(V, i), Row2D(V, i1), Row2D(V, i2), Row2D(V, i3)))
          return true;
      }
    }

    return false;
  }


  inline bool InsideTriangle2D(const EMatXd& V, const EVec3i& triangle)
  {
    for (int i = 0; i < V.rows(); i++)
    {
      if (i != triangle[0] && i != triangle[1] && i != triangle[2])
      {
        if (InsideTriangle2D(Row2D(V, i), Row2D(V, triangle[0]), Row2D(V, triangle[1]), Row2D(V, triangle[2])))
          return true;
      }
    }

    return false;
  }


  inline bool FaceForwardTriangle2D(const EVec2d& a, const EVec2d& b, const EVec2d& c)
  {
    EVec2d ab = b - a;
    EVec2d ac = c - a;

    if ((ab[0] * ac[1] - ab[1] * ac[0]) > 0.0)
      return true;

    return false;
  }

  inline bool FaceBackwardTriangle2D(const EVec2d& a, const EVec2d& b, const EVec2d& c)
  {
    EVec2d ab = b - a;
    EVec2d ac = c - a;

    if ((ab[0] * ac[1] - ab[1] * ac[0]) < 0.0)
      return true;

    return false;
  }


  inline bool FaceForwardTriangle2D(const EMatXd& V, const EVec3i& triangle)
  {
    EVec2d a = EVec2d(V(triangle[0], 0), V(triangle[0], 1));
    EVec2d b = EVec2d(V(triangle[1], 0), V(triangle[1], 1));
    EVec2d c = EVec2d(V(triangle[2], 0), V(triangle[2], 1));

    return FaceForwardTriangle2D(a, b, c);
  }


  inline bool FaceBackwardTriangle2D(const EMatXd& V, const EVec3i& triangle)
  {
    EVec2d a = EVec2d(V(triangle[0], 0), V(triangle[0], 1));
    EVec2d b = EVec2d(V(triangle[1], 0), V(triangle[1], 1));
    EVec2d c = EVec2d(V(triangle[2], 0), V(triangle[2], 1));

    return FaceBackwardTriangle2D(a, b, c);
  }
}





