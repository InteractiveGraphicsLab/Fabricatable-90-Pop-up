/* FabricatablePatch.cpp */
/* Fabricatable patch classes' function implementation */
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

#include "FabricatablePatch.h"


using std::vector;
using std::pair;


FabricatablePatch::FabricatablePatch()
{
	type_ = E_FACE_TYPE::DEFAULT;
	mount_face_ = false;
	loaded_mesh_ = false;
  skip_ = false;
  rel_plane_ = JMesh::Mesh();
  rel_original_plane_ = JMesh::Mesh();
  rel_fill_plane_ = JMesh::Mesh();
  rel_original_fill_plane_ = JMesh::Mesh();
  rel_outside_ = JMesh::Mesh();
  rel_patch_ = JMesh::Mesh();
}


FabricatablePatch::FabricatablePatch(const E_FACE_TYPE& type, bool mount_face)
{
	type_ = type;
	mount_face_ = mount_face;
	loaded_mesh_ = false;
  skip_ = false;
  rel_plane_ = JMesh::Mesh();
  rel_original_plane_ = JMesh::Mesh();
  rel_fill_plane_ = JMesh::Mesh();
  rel_original_fill_plane_ = JMesh::Mesh();
  inflated_plane_ = JMesh::Mesh();
  rel_outside_ = JMesh::Mesh();
  rel_patch_ = JMesh::Mesh();
}


FabricatablePatch::FabricatablePatch(const FabricatablePatch& src)
{
	this->type_ = src.type_;
	this->mount_face_ = src.mount_face_;
	this->loaded_mesh_ = src.loaded_mesh_;
  this->skip_ = src.skip_;
  this->rel_plane_ = src.rel_plane_;
  this->rel_original_plane_ = src.rel_original_plane_;
  this->rel_fill_plane_ = src.rel_fill_plane_;
  this->rel_original_fill_plane_ = src.rel_original_fill_plane_;
  this->inflated_plane_ = src.inflated_plane_;
  this->rel_outside_ = src.rel_outside_;
  this->rel_patch_ = src.rel_patch_;
}


void FabricatablePatch::Plane
(
  const JMesh::Mesh& plane,
  const FaceCoordSystem& face_system,
  const EVec3d& depend_org,
  const vector<EVec3d>& rect,
  const ConvertProperties& cvt_props
)
{
  JMesh::Mesh concave_cut_box = GenerateConcaveCutBox180(rect, cvt_props);
  rel_plane_ = JCSG::Subtract(plane, concave_cut_box);

  JMesh::RelSystem rel_system = { depend_org, face_system.u_dir, face_system.v_dir, face_system.n_dir };
  rel_plane_ = JMesh::RelMesh(rel_plane_, rel_system);
  rel_original_plane_ = rel_plane_;
}


void FabricatablePatch::FillPlane
(
  const JMesh::Mesh& plane,
  const FaceCoordSystem& face_system,
  const EVec3d& depend_org,
  const vector<EVec3d>& rect,
  const ConvertProperties& cvt_props
)
{
  JMesh::Mesh concave_cut_box = GenerateConcaveCutBox180(rect, cvt_props);
  rel_fill_plane_ = JCSG::Subtract(plane, concave_cut_box);

  JMesh::RelSystem rel_system = { depend_org, face_system.u_dir, face_system.v_dir, face_system.n_dir };
  rel_fill_plane_ = JMesh::RelMesh(rel_fill_plane_, rel_system);
  rel_original_fill_plane_ = rel_fill_plane_;
}


void FabricatablePatch::Patch(const JMesh::Mesh& patch)
{
  rel_patch_ = patch;
  rel_plane_ = JMesh::Mesh();
  rel_outside_ = JMesh::Mesh();
  inflated_plane_ = JMesh::Mesh();
}


void FabricatablePatch::ResetMeshs()
{
  if (!mount_face_)
    rel_plane_= JMesh::Mesh();

  rel_original_plane_ = JMesh::Mesh();
  rel_fill_plane_ = JMesh::Mesh();
  rel_original_fill_plane_ = JMesh::Mesh();
  inflated_plane_ = JMesh::Mesh();
  rel_outside_ = JMesh::Mesh();
}


void FabricatablePatch::Draw3D(const FaceCoordSystem& face_system, const EVec3d& depend_org, const EVec3f& color) const
{
  JMesh::RelSystem rel_system = { depend_org, face_system.u_dir, face_system.v_dir, face_system.n_dir };

  if (rel_patch_.Exist())
    JMesh::StdMesh(rel_patch_, rel_system).DrawPlane(color);
  else
    JMesh::StdMesh(rel_plane_, rel_system).DrawPlane(color);
}


static JMesh::RelSystem sGetRelSystemAngle90(const E_FACE_TYPE& type, const EVec3d& org)
{
	JMesh::RelSystem rel_system;
	if (type == E_FACE_TYPE::VTYPE)
		rel_system = { org, EVec3d::UnitX(), EVec3d::UnitY(), EVec3d::UnitZ() };
	else if (type == E_FACE_TYPE::HTYPE)
		rel_system = { org, EVec3d::UnitX(), EVec3d::UnitZ(), -EVec3d::UnitY() };

	return rel_system;
}


void FabricatablePatch::CutPlane
(
	const JMesh::Mesh& mesh,
	const vector<EVec3d>& rect,
  const Fold& depend_fold,
	const ConvertProperties& cvt_props
)
{
  std::cout << "  Generate ";
  if (type_ == E_FACE_TYPE::VTYPE)
    std::cout << "V ";
  else if (type_ == E_FACE_TYPE::HTYPE)
    std::cout << "H ";
  std::cout << "Patch Mesh";

  JMesh::RelSystem rel_system = sGetRelSystemAngle90(type_, rect[0]);
  if (mount_face_)
	{
    rel_plane_ = JMesh::StdMesh(rel_plane_, rel_system);
    rel_original_plane_ = JMesh::StdMesh(rel_original_plane_, rel_system);
    rel_fill_plane_ = JMesh::StdMesh(rel_fill_plane_, rel_system);
    rel_original_fill_plane_ = JMesh::StdMesh(rel_original_fill_plane_, rel_system);
		loaded_mesh_ = true;
	}
  else
  {
    JMesh::Mesh convex = GenerateConvex(rect, cvt_props);
    JMesh::Mesh convex_cut_box = GenerateConvexCutBox(rect, cvt_props);

    rel_original_plane_ = JMesh::StdMesh(rel_original_plane_, rel_system);
    rel_original_fill_plane_ = JMesh::StdMesh(rel_original_fill_plane_, rel_system);

    igl::copyleft::cgal::CSGTree original_plane_tree =
    {{{rel_original_plane_.V(), rel_original_plane_.F()}, {convex.V(), convex.F()}, "u"},
     {convex_cut_box.V(), convex_cut_box.F()}, "m" };
    std::cout << ".";

    igl::copyleft::cgal::CSGTree original_fill_plane_tree =
    {{{rel_original_fill_plane_.V(), rel_original_fill_plane_.F()}, {convex.V(), convex.F()}, "u"},
      {convex_cut_box.V(), convex_cut_box.F()}, "m" };
    std::cout << ".";

    rel_original_plane_ = JMesh::Mesh(original_plane_tree.cast_V<EMatXd>(), original_plane_tree.F());
    rel_original_fill_plane_ = JMesh::Mesh(original_fill_plane_tree.cast_V<EMatXd>(), original_fill_plane_tree.F());

    rel_plane_ = JCSG::Intersect(mesh, rel_original_plane_);
    std::cout << ".";
    rel_fill_plane_ = JCSG::Intersect(mesh, rel_original_fill_plane_);
    std::cout << ".";

    loaded_mesh_ = true;
	}

	JMesh::Mesh concave_cut_box = GenerateConcaveCutBox90(rect, depend_fold, cvt_props);
  rel_plane_ = JCSG::Subtract(rel_plane_, concave_cut_box);
  std::cout << ".";
  rel_original_plane_ = JCSG::Subtract(rel_original_plane_, concave_cut_box);
  std::cout << ".";
  rel_fill_plane_ = JCSG::Subtract(rel_fill_plane_, concave_cut_box);
  std::cout << ".";
  rel_original_fill_plane_ = JCSG::Subtract(rel_original_fill_plane_, concave_cut_box);
  std::cout << ".";

	rel_plane_ = JMesh::RelMesh(rel_plane_, rel_system);
  rel_original_plane_ = JMesh::RelMesh(rel_original_plane_, rel_system);
  rel_fill_plane_ = JMesh::RelMesh(rel_fill_plane_, rel_system);
  rel_original_fill_plane_ = JMesh::RelMesh(rel_original_fill_plane_, rel_system);

  std::cout << " Done" << std::endl;
}


void FabricatablePatch::GenerateInflatedPatch
(
	JMesh::Mesh& clip_plane,
	const vector<EVec3d> rect,
	const FaceCoordSystem& face_system,
	const EVec3d& depend_org,
	const ConvertProperties& cvt_props
)
{
	if (mount_face_)
		return;

	if (loaded_mesh_)
	{
    bool not_inflate_plane = !inflated_plane_.Exist();

    if (not_inflate_plane)
      std::cout << "Generate Inflate Patch";

    JMesh::Mesh inflate_plane = InfatePlane(rect, face_system, depend_org, cvt_props);
    clip_plane = JCSG::Union(clip_plane, inflate_plane);

    if (not_inflate_plane)
      std::cout << ".";

    if (not_inflate_plane)
      std::cout << "Done" << std::endl;
	}
}


void FabricatablePatch::TrimOutlines
(
	const JMesh::Mesh& outlines,
	const FaceCoordSystem& face_system,
	const EVec3d& depend_org,
	const ConvertProperties& cvt_props
)
{
	if (loaded_mesh_)
	{
    std::cout << "Trimming Patch";
		JMesh::RelSystem rel_system = { depend_org, face_system.u_dir, face_system.v_dir, face_system.n_dir };

		JMesh::Mesh plane = JMesh::StdMesh(rel_fill_plane_, rel_system);
		plane = JCSG::Subtract(plane, outlines);
    std::cout << ".";

    rel_plane_ = JMesh::RelMesh(plane, rel_system);
    std::cout << "Done" << std::endl;
	}
}


static JMesh::Mesh sCalcSpace(const EVec3d& max_pos, const vector<EVec3d>& rect, const E_FACE_TYPE& type, const ConvertProperties& cvt_props)
{
  EVec3d org = (type == E_FACE_TYPE::HTYPE) ?
    EVec3d(rect[0][0], max_pos[1], rect[0][2] + cvt_props.fold_gap) :
    EVec3d(rect[3][0], rect[3][1] - (0.5 * cvt_props.fold_gap), rect[3][2]);

  EVec3d x =  EVec3d::UnitX();
  EVec3d y = -EVec3d::UnitY();
  EVec3d z =  EVec3d::UnitZ();

  double w = (rect[1] - rect[0]).norm();
  double h = (type == E_FACE_TYPE::HTYPE) ?
    fabs(max_pos[1] - rect[0][1]) : ((rect[3] - rect[0]).norm() - (1.5 * cvt_props.fold_gap));
  double d = (type == E_FACE_TYPE::HTYPE) ?
    ((rect[3] - rect[0]).norm() - (1.5 * cvt_props.fold_gap)) : fabs(max_pos[2] - rect[0][2]);

  JMesh::Rectangular space = { org, x, y, z, w, h, d };
  return JMesh::RectangularMesh(space);
}

static vector<int> sGetConcaveFoldIndexList(const EMatXd& V, const vector<pair<double, double>>& concave_fold_list)
{
  vector<int> index_list(0);

  for (int i = 0; i < V.rows(); i++)
  {
    for (int j = 0; j < concave_fold_list.size(); j++)
    {
      if (((concave_fold_list[j].first - (1e-5)) <= V(i, 0)) &&
          ((concave_fold_list[j].second + (1e-5)) >= V(i, 0)))
      {

        if ((fabs(V(i, 1)) <= 1e-5) && (fabs(V(i, 2)) <= 1e-5))
          index_list.push_back(i);
      }
    }
  }

  return index_list;
}

static vector<int> sGetConvexFoldIndexList(const EMatXd& V, double y_max)
{
  vector<int> index_list(0);

  for (int i = 0; i < V.rows(); i++)
  {
    if ((fabs(V(i, 1) - y_max) <= 1e-5) && (fabs(V(i, 2)) <= 1e-5))
      index_list.push_back(i);
  }

  return index_list;
}

static bool sIsConnectFolds(int v_idx, const vector<int>& convex_list)
{
  for (int convex_idx : convex_list)
  {
    if (convex_idx == v_idx)
      return true;
  }

  return false;
}

static int sSearchVLabel0(const vector<int>& labels)
{
  for (int i = 0; i < labels.size(); i++)
  {
    if (labels[i] == 0)
      return i;
  }

  return -1;
}

static JMesh::Mesh sRemoveIsolateMesh(const JMesh::Mesh& mesh, double y_max, const vector<pair<double, double>>& concave_fold_list)
{
  vector<int> v_labels(mesh.V().rows(), -1);
  vector<int> f_labels(mesh.F().rows(), -1);
  vector<EVec3i> f(0);

  vector<int> concave_list = sGetConcaveFoldIndexList(mesh.V(), concave_fold_list);
  vector<int> convex_list = sGetConvexFoldIndexList(mesh.V(), y_max);

  if (concave_list.size() == 0 || convex_list.size() == 0)
  {
    if (concave_list.size() == 0)
      std::cout << "h";
    if (convex_list.size() == 0)
      std::cout << "s";

    return JMesh::Mesh();
  }
  bool connected = true;
  for (int k = 0; k < concave_list.size(); k++)
  {
    int connected_count = 0;
    if (v_labels[concave_list[k]] == 1)
      continue;

    int v_idx = concave_list[k];
    while (v_idx != -1)
    {
      for (int i = 0; i < f_labels.size(); i++)
      {
        if (f_labels[i] != -1)
          continue;

        if ((mesh.F(i, 0) == v_idx) || (mesh.F(i, 1) == v_idx) || (mesh.F(i, 2) == v_idx))
        {
          for (int j = 0; j < 3; j++)
          {
            if (v_labels[mesh.F(i, j)] == -1)
              v_labels[mesh.F(i, j)] = 0;
          }

          f_labels[i] = 1;
          f.push_back(mesh.F(i));
        }
      }

      if (sIsConnectFolds(v_idx, convex_list))
        connected_count++;

      v_labels[v_idx] = 1;
      v_idx = sSearchVLabel0(v_labels);
    }

    if (connected_count == 0)
      connected = false;
  }

  if (connected)
  {
    EMatXi F = JUtil::ToEMatXi(f);
    return JMesh::Mesh(mesh.V(), F);
  }
  else
  {
    return JMesh::Mesh();
  }
}

static vector<EVec3i> sGetTriangles(const JMesh::Mesh& mesh, int idx)
{
  vector<EVec3i> triangle_list(0);
  for (int i = 0; i < mesh.F().rows(); i++)
  {
    if (mesh.F(i, 0) == idx || mesh.F(i, 1) == idx || mesh.F(i, 2) == idx)
      triangle_list.push_back(EVec3i(mesh.F(i, 0), mesh.F(i, 1), mesh.F(i, 2)));
  }

  return triangle_list;
}

static bool sIsLongerMinimumConcaveFolds(const JMesh::Mesh& mesh, const ConvertProperties& cvt_props, const vector<pair<double, double>>& concave_fold_list)
{
  for (int i = 0; i < concave_fold_list.size(); i++)
  {
    EVec3d org = EVec3d(concave_fold_list[i].first - J_CSG_OFFSET, cvt_props.fold_gap, -cvt_props.fold_thickness - J_CSG_OFFSET);
    double width = concave_fold_list[i].second - concave_fold_list[i].first + 2.0 * J_CSG_OFFSET;
    double height = cvt_props.fold_gap + J_CSG_OFFSET;
    double depth = 2.0 * cvt_props.fold_thickness + 2.0 * J_CSG_OFFSET;
    JMesh::Rectangular concave_box = { org, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };
    JMesh::Mesh concave = JCSG::Intersect(mesh, JMesh::RectangularMesh(concave_box));

    if (fabs(concave.GetBoundBox().MaxPos()[0] - concave.GetBoundBox().MinPos()[0]) >= cvt_props.nozzle_width)
      return true;
  }

  return false;
}

static bool sIsLongerMinimumConvexFolds(const JMesh::Mesh& mesh, double y_max, const ConvertProperties& cvt_props)
{
  EVec3d org = EVec3d(-J_CSG_OFFSET, y_max + J_CSG_OFFSET, -cvt_props.fold_thickness - J_CSG_OFFSET);
  double width = mesh.GetBoundBox().MaxPos()[0] + 2.0 * J_CSG_OFFSET;
  double height = cvt_props.fold_gap + J_CSG_OFFSET;
  double depth = 2.0 * cvt_props.fold_thickness + 2.0 * J_CSG_OFFSET;
  JMesh::Rectangular convex_box = { org, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };
  JMesh::Mesh convex = JCSG::Intersect(mesh, JMesh::RectangularMesh(convex_box));

  if (fabs(convex.GetBoundBox().MaxPos()[0] - convex.GetBoundBox().MinPos()[0]) >= cvt_props.nozzle_width)
    return true;
  else
    return false;
}


bool FabricatablePatch::SegmentSweepRegion
(
  const JMesh::Mesh& higher_space,
  const JMesh::Mesh& mesh,
  const vector<EVec3d>& rect,
  const FaceCoordSystem& face_system,
  const ConvertProperties& cvt_props,
  const EVec3d& max_pos,
  const vector<pair<double, double>>& concave_fold_list
)
{
  std::cout << "Segment Sweep Region";
  JMesh::Mesh space = sCalcSpace(max_pos, rect, type_, cvt_props);
  JMesh::Mesh segment_region = JCSG::Subtract(space, higher_space);
  std::cout << ".";
  rel_outside_ = JCSG::Intersect(segment_region, mesh);
  std::cout << ".";

  JMesh::RelSystem rel_system = { rect[0], face_system.u_dir, face_system.v_dir, face_system.n_dir};
  rel_outside_ = JMesh::RelMesh(rel_outside_, rel_system);

  if (type_ == E_FACE_TYPE::VTYPE)
  {
    double plane_move = 1e-6 - rel_plane_.GetBoundBox().MaxPos()[2];
    double outside_move = -1e-6 - rel_outside_.GetBoundBox().MinPos()[2];
    rel_plane_.Move(EVec3d(0.0, 0.0, plane_move));
    rel_outside_.Move(EVec3d(0.0, 0.0, outside_move));
  }
  else if (type_ == E_FACE_TYPE::HTYPE)
  {
    double plane_move = -1e-6 - rel_plane_.GetBoundBox().MinPos()[2];
    double outside_move = 1e-6 - rel_outside_.GetBoundBox().MaxPos()[2];
    rel_plane_.Move(EVec3d(0.0, 0.0, plane_move));
    rel_outside_.Move(EVec3d(0.0, 0.0, outside_move));
  }

  rel_patch_ = JCSG::Union(rel_plane_, rel_outside_);
  std::cout << ".";

  double y_max = (rect[0] - rect[3]).norm();
  rel_patch_ = sRemoveIsolateMesh(rel_patch_, y_max, concave_fold_list);
  std::cout << ".";

  if (!skip_)
  {
    if (!rel_patch_.Exist() ||
        !sIsLongerMinimumConcaveFolds(rel_patch_, cvt_props, concave_fold_list) ||
        !sIsLongerMinimumConvexFolds(rel_patch_, y_max, cvt_props))
    {
      rel_plane_ = rel_original_plane_;
      rel_fill_plane_ = rel_original_fill_plane_;

      inflated_plane_ = JMesh::Mesh();
      std::cout << "Done2" << std::endl;
      return true;
    }

    //rel_fill_plane_ = JCSG::Intersect(rel_patch_, rel_original_fill_plane_);
    //rel_plane_ = rel_fill_plane_;

    rel_original_plane_ = JMesh::Mesh();
    rel_original_fill_plane_ = JMesh::Mesh();
    rel_outside_ = JMesh::Mesh();
    inflated_plane_ = JMesh::Mesh();

    //rel_patch_ = JMesh::Mesh();
    skip_ = true;
  }

  std::cout << "Done" << std::endl;
  return false;
}


JMesh::Mesh FabricatablePatch::Assemble(const FaceCoordSystem& face_system, const EVec3d& depend_org) const
{
  JMesh::RelSystem rel_system = { depend_org, face_system.u_dir, face_system.v_dir, face_system.n_dir };

  if (!mount_face_)
    return JMesh::StdMesh(rel_patch_, rel_system);
  else
    return JMesh::StdMesh(rel_plane_, rel_system);
}


static JMesh::Mesh sMakeFatProjection(const JMesh::Mesh &projection, const EVec3d &n_dir, double fat_radius, int slice)
{
	JMesh::Mesh fat_projection = projection;
	vector<EVec3d> projection_3d = JUtil::ToStdVectorEVec3d(projection.V());
	EMatXi F = projection.F();

	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			int idx1 = j;
			int idx2 = (j == (F.cols() - 1)) ? 0 : (j + 1);

			EVec3d pos1 = projection_3d[F(i, idx1)];
			EVec3d pos2 = projection_3d[F(i, idx2)];

			EVec3d v_dir = (pos2 - pos1).normalized();
			EVec3d u_dir = n_dir.cross(v_dir).normalized();

			JMesh::FatLine2D fat_line_plane = { pos1, pos2, u_dir, v_dir, fat_radius, slice };
			fat_projection.Compose((JMesh::FatLine2DMesh(fat_line_plane)));
		}
	}

	return fat_projection;
}

static cv::Mat sProjectImage(const vector<EVec3d>& rect, const JMesh::Mesh& projection,	const EVec3d& depend_org,	const ConvertProperties& cvt_props, const E_FACE_TYPE& type, const double AMPL)
{
	const int WIDTH  = int(AMPL * ((rect[3] - rect[0]).norm() + (0.5 * cvt_props.fold_gap)));
	const int HEIGHT = int(AMPL * ((rect[1] - rect[0]).norm() + (2.0 * cvt_props.gap)));
	cv::Mat projection_img = cv::Mat::zeros(WIDTH, HEIGHT, CV_8UC1);

	for (int i = 0; i < projection.F().rows(); i++)
	{
		vector<vector<cv::Point>> polygon_points(1, vector<cv::Point>(3));
		for (int j = 0; j < projection.F().cols(); j++)
		{
			int x_img = int(AMPL * fabs(projection.V(projection.F(i, j), 0) - depend_org[0] + cvt_props.gap));
			int y_img = int(AMPL * fabs(projection.V(projection.F(i, j), 2) - depend_org[2]));
			polygon_points[0][j] = cv::Point(x_img, y_img);
		}

  	cv::fillPoly(projection_img, polygon_points, 255);
	}

	return projection_img;
}

static vector<EVec2d> sToStdVectorContours(const vector<cv::Point>& cntrs, const E_FACE_TYPE& type, const EVec3d& depend_org, const ConvertProperties& cvt_props, const double AMPL)
{
	vector<EVec2d> cntrs_pos_2d = vector<EVec2d>(cntrs.size());
  for (int i = 0; i < cntrs.size(); i++)
  {
    cntrs_pos_2d[i] = EVec2d(double(cntrs[i].x), double(cntrs[i].y));
    cntrs_pos_2d[i] /= AMPL;

    if (type == E_FACE_TYPE::VTYPE)
    {
      cntrs_pos_2d[i][1] *= -1.0;
      cntrs_pos_2d[i] += EVec2d(depend_org[0] - cvt_props.gap,
                                depend_org[2]);
    }
    else if (type == E_FACE_TYPE::HTYPE)
    {
      cntrs_pos_2d[i] += EVec2d(depend_org[0] - cvt_props.gap,
                                depend_org[2]);
    }
  }

	return cntrs_pos_2d;
}

static EMatXi sGenerateContoursEdge(const EMatXd& cntrs_2d_V)
{
	EMatXi cntrs_2d_E(cntrs_2d_V.rows(), 2);
	for (int j = 0; j < cntrs_2d_E.rows(); j++)
	{
		cntrs_2d_E(j, 0) = j;
		cntrs_2d_E(j, 1) = (j != cntrs_2d_E.rows() - 1) ? (j + 1) : 0;
	}

	return cntrs_2d_E;
}

static EMatXi sAlignNormals(const EMatXd& tri_cntrs_2d_V, const EMatXi& tri_cntrs_2d_F)
{
	EVec3d p0 = EVec3d(tri_cntrs_2d_V(tri_cntrs_2d_F(0, 0), 0), 0.0, tri_cntrs_2d_V(tri_cntrs_2d_F(0, 0), 1));
	EVec3d p1 = EVec3d(tri_cntrs_2d_V(tri_cntrs_2d_F(0, 1), 0), 0.0, tri_cntrs_2d_V(tri_cntrs_2d_F(0, 1), 1));
	EVec3d p2 = EVec3d(tri_cntrs_2d_V(tri_cntrs_2d_F(0, 2), 0), 0.0, tri_cntrs_2d_V(tri_cntrs_2d_F(0, 2), 1));
	EVec3d normal = (p1 - p0).cross(p2 - p1).normalized();

	EMatXi aligned_F = tri_cntrs_2d_F;

	if (normal[1] < -J_DBL_EPSILON)
	{
		aligned_F(Eigen::seqN(0, tri_cntrs_2d_F.rows()), 0) = tri_cntrs_2d_F(Eigen::seqN(0, tri_cntrs_2d_F.rows()), 2);
		aligned_F(Eigen::seqN(0, tri_cntrs_2d_F.rows()), 2) = tri_cntrs_2d_F(Eigen::seqN(0, tri_cntrs_2d_F.rows()), 0);
	}

	return aligned_F;
}

static EMatXd sAddThickClipPlaneV(const EMatXd& tri_cntrs_2d_V, double thickness, const JMesh::Mesh& rel_outside)
{
	const int V_ROW = tri_cntrs_2d_V.rows();
	EMatXd thicked_V(V_ROW * 2, 3);

	thicked_V(Eigen::seqN(0, V_ROW), 0) = tri_cntrs_2d_V(Eigen::seqN(0, V_ROW), 0);
	thicked_V(Eigen::seqN(0, V_ROW), 1) = EMatXd::Constant(V_ROW, 1, J_CSG_OFFSET);
	thicked_V(Eigen::seqN(0, V_ROW), 2) = tri_cntrs_2d_V(Eigen::seqN(0, V_ROW), 1);

	thicked_V(Eigen::seqN(V_ROW, V_ROW), 0) = tri_cntrs_2d_V(Eigen::seqN(0, V_ROW), 0);
	thicked_V(Eigen::seqN(V_ROW, V_ROW), 1) = EMatXd::Constant(V_ROW, 1, -thickness - J_CSG_OFFSET);
	thicked_V(Eigen::seqN(V_ROW, V_ROW), 2) = tri_cntrs_2d_V(Eigen::seqN(0, V_ROW), 1);

	return thicked_V;
}

static bool sIsClockwise(const EMatXd& tri_cntrs_2d_V)
{
  double area = 0.0;
  for (int i = 0; i < tri_cntrs_2d_V.rows(); i++)
  {
    int next = (i == tri_cntrs_2d_V.rows() - 1) ? 0 : (i + 1);
    EVec2d p      = EVec2d(tri_cntrs_2d_V(i   , 0), tri_cntrs_2d_V(i   , 1));
    EVec2d p_next = EVec2d(tri_cntrs_2d_V(next, 0), tri_cntrs_2d_V(next, 1));

    area += p[0] * p_next[1] - p_next[0] * p[1];
  }

  return area < 0.0;
}

static EMatXi sAddThickClipPlaneF(const EMatXd& tri_cntrs_2d_V, const EMatXi& tri_cntrs_2d_F, const EMatXd& V, const vector<cv::Point>& cntrs)
{
	const int V_ROW = tri_cntrs_2d_V.rows();
	const int F_ROW = tri_cntrs_2d_F.rows();
	EMatXi thicked_F((F_ROW * 2) + (V_ROW * 2), 3);

	thicked_F(Eigen::seqN(0, F_ROW), Eigen::seqN(0, 3)) = tri_cntrs_2d_F;

	thicked_F(Eigen::seqN(F_ROW, F_ROW), 0) = tri_cntrs_2d_F(Eigen::seqN(0, F_ROW), 2) + EMatXi::Constant(F_ROW, 1, V_ROW);
	thicked_F(Eigen::seqN(F_ROW, F_ROW), 1) = tri_cntrs_2d_F(Eigen::seqN(0, F_ROW), 1) + EMatXi::Constant(F_ROW, 1, V_ROW);
	thicked_F(Eigen::seqN(F_ROW, F_ROW), 2) = tri_cntrs_2d_F(Eigen::seqN(0, F_ROW), 0) + EMatXi::Constant(F_ROW, 1, V_ROW);

	bool clockwise = sIsClockwise(tri_cntrs_2d_V);
	for (int i = 0; i < V_ROW; i++)
	{
		int idx = i * 2;

		thicked_F((F_ROW * 2) + idx, 0) = (clockwise) ? (i + V_ROW) : i;
		thicked_F((F_ROW * 2) + idx, 1) = (i != (V_ROW - 1)) ? (i + 1) : 0;
		thicked_F((F_ROW * 2) + idx, 2) = (clockwise) ? i : (i + V_ROW);

		int j = (i != (V_ROW - 1)) ? (i + 1) : 0;

		thicked_F((F_ROW * 2) + idx + 1, 0) = (clockwise) ? (i + V_ROW) : j;
		thicked_F((F_ROW * 2) + idx + 1, 1) = (i != (V_ROW - 1)) ? (i + 1 + V_ROW) : V_ROW;
		thicked_F((F_ROW * 2) + idx + 1, 2) = (clockwise) ? j : (i + V_ROW);
	}

	return thicked_F;
}


JMesh::Mesh FabricatablePatch::InfatePlane
(
  const vector<EVec3d> rect,
	const FaceCoordSystem& face_system,
	const EVec3d& depend_org,
	const ConvertProperties& cvt_props
)
{
  if (inflated_plane_.Exist())
    return inflated_plane_;

	if (loaded_mesh_ && !mount_face_)
	{
		JMesh::RelSystem rel_system = { depend_org, face_system.u_dir, face_system.v_dir, face_system.n_dir };
    JMesh::Mesh plane = JMesh::StdMesh(rel_plane_, rel_system);

    vector<EVec3d> projected_plane_3d = JMesh::ProjectMesh(plane, rel_system);
    std::cout << ".";

		JMesh::Mesh projected_plane = JMesh::Mesh(JUtil::ToEMatXd(projected_plane_3d), plane.F());
		projected_plane = sMakeFatProjection(projected_plane, face_system.n_dir, cvt_props.gap, 8);
    std::cout << ".";

		const double AMPL = 100.0;
		cv::Mat projection_img = sProjectImage(rect, projected_plane, depend_org, cvt_props, type_, AMPL);
    std::cout << ".";

		vector<vector<cv::Point>> cntrs;
		vector<cv::Vec4i> hierarchy;
		cv::findContours(projection_img, cntrs, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    std::cout << ".";

		EMatXd V(0, 3);
		EMatXi F(0, 3);
		JMesh::Mesh inflated(V, F);
		for (int i = 0; i < cntrs.size(); i++)
		{
			vector<EVec2d> cntrs_pos_2d = sToStdVectorContours(cntrs[i], type_, depend_org, cvt_props, AMPL);

			EMatXd cntrs_2d_H(0, 2);
			EMatXd cntrs_2d_V = JUtil::ToEMatXd(cntrs_pos_2d);
			EMatXi cntrs_2d_E = sGenerateContoursEdge(cntrs_2d_V);

			EMatXd tri_cntrs_2d_V;
			EMatXi tri_cntrs_2d_F;
			igl::triangle::triangulate(cntrs_2d_V, cntrs_2d_E, cntrs_2d_H, "a1000.0", tri_cntrs_2d_V, tri_cntrs_2d_F);
			tri_cntrs_2d_F = sAlignNormals(tri_cntrs_2d_V, tri_cntrs_2d_F);

			EMatXd inflated_V = sAddThickClipPlaneV(tri_cntrs_2d_V, cvt_props.thickness, rel_outside_);
			EMatXi inflated_F = sAddThickClipPlaneF(tri_cntrs_2d_V, tri_cntrs_2d_F, inflated_V, cntrs[i]);
			inflated.Compose(JMesh::Mesh(inflated_V, inflated_F));
		}
    std::cout << ".";

    inflated_plane_ = inflated;
		return inflated;
	}
}


JMesh::Mesh FabricatablePatch::GenerateConvex(const vector<EVec3d>& rect, const ConvertProperties& cvt_props) const
{
  JMesh::Rectangular convex;
  double width = rect[1][0] - rect[0][0];
  if (type_ == E_FACE_TYPE::VTYPE)
  {
    EVec3d pos = rect[3] - EVec3d(0.0, 0.0, cvt_props.fold_thickness);
    double height = cvt_props.fold_gap / 2.0;
    convex = { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, cvt_props.fold_thickness };
  }
  else
  {
    EVec3d pos = rect[3] - EVec3d(0.0, 0.0, (cvt_props.fold_gap / 2.0));
    double depth = cvt_props.fold_gap / 2.0;
    convex = { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, cvt_props.fold_thickness, depth };
  }

  return JMesh::RectangularMesh(convex);
}


JMesh::Mesh FabricatablePatch::GenerateConvexCutBox(const std::vector<EVec3d>& rect, const ConvertProperties& cvt_props) const
{
  EVec3d pos = JUtil::ErrEVec3d();
  double width = rect[2][0] - rect[3][0];
  double height = -1.0;
  double depth = -1.0;

  if (type_ == E_FACE_TYPE::VTYPE)
  {
    pos = rect[3] + EVec3d(0.0, J_CSG_OFFSET, -(cvt_props.thickness + J_CSG_OFFSET));
    height = (cvt_props.fold_gap / 2.0) + J_CSG_OFFSET;
    depth = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
  }
  else if (type_ == E_FACE_TYPE::HTYPE)
  {
    pos = rect[3] - EVec3d(0.0, cvt_props.fold_thickness, cvt_props.fold_gap / 2.0);
    height = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
    depth = cvt_props.fold_gap / 2.0 + J_CSG_OFFSET;
  }

  JMesh::Rectangular convex_cut_box =
  { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };

  return JMesh::RectangularMesh(convex_cut_box);
}


static JMesh::Mesh sMountBaseFoldCutBox(const vector<EVec3d>& rect, const E_FACE_TYPE& type, const ConvertProperties& cvt_props)
{
  EVec3d pos = JUtil::ErrEVec3d();
  double width = rect[1][0] - rect[0][0] + 2.0 * J_CSG_OFFSET;
  double height = -1.0;
  double depth = -1.0;
  if (type == E_FACE_TYPE::VTYPE)
  {
    pos = rect[0] - EVec3d(J_CSG_OFFSET, -(cvt_props.fold_gap / 2.0), cvt_props.thickness + J_CSG_OFFSET);
    height = cvt_props.fold_gap / 2.0 + J_CSG_OFFSET;
    depth = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
  }
  else if (type == E_FACE_TYPE::HTYPE)
  {
    pos = rect[0] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, J_CSG_OFFSET);
    height = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
    depth = cvt_props.fold_gap / 2.0 + J_CSG_OFFSET;
  }

  JMesh::Rectangular mount_base_fold_cut_box =
  { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };

  return JMesh::RectangularMesh(mount_base_fold_cut_box);
}

static JMesh::Mesh sBaseFoldCutBox(const vector<EVec3d>& rect, const E_FACE_TYPE& type, const ConvertProperties& cvt_props)
{
  EVec3d pos = JUtil::ErrEVec3d();
  double width = rect[1][0] - rect[0][0] + 2.0 * J_CSG_OFFSET;
  double height = -1.0;
  double depth = -1.0;

  if (type == E_FACE_TYPE::VTYPE)
  {
    pos = rect[0] + EVec3d(-J_CSG_OFFSET, cvt_props.fold_gap, -(cvt_props.thickness + J_CSG_OFFSET));
    height = cvt_props.fold_gap + J_CSG_OFFSET;
    depth = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
  }
  else if (type == E_FACE_TYPE::HTYPE)
  {
    pos = rect[0] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, J_CSG_OFFSET);
    height = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
    depth = cvt_props.fold_gap + J_CSG_OFFSET;
  }

  JMesh::Rectangular base_fold_cut_box =
  { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };

  return JMesh::RectangularMesh(base_fold_cut_box);
}


JMesh::Mesh FabricatablePatch::GenerateConcaveCutBox90(const vector<EVec3d>& rect, const Fold& depend_fold,	const ConvertProperties& cvt_props) const
{
  JMesh::Mesh cut_box = JMesh::Mesh();
  if (mount_face_)
  {
    cut_box = sMountBaseFoldCutBox(rect, type_, cvt_props);
  }
  else
  {
    cut_box = sBaseFoldCutBox(rect, type_, cvt_props);

    if ((static_cast<int>(depend_fold.type) != static_cast<int>(type_)) &&
        (depend_fold.type != E_FOLD_TYPE::GTYPE))
    {
      EVec3d pos = JUtil::ErrEVec3d();
      double width = J_CSG_OFFSET;
      double height = -1.0;
      double depth = -1.0;

      double parent_left = depend_fold.org[0];
      double parent_right = depend_fold.org[0] + depend_fold.u;
      if (parent_left > rect[0][0])
      {
        pos = rect[0] - EVec3d(J_CSG_OFFSET, 0.0, 0.0);
        width += parent_left - rect[0][0];
      }
      else if (parent_right < rect[1][0])
      {
        pos = EVec3d(parent_right, rect[1][1], rect[1][2]);
        width += rect[1][0] - parent_right;
      }

      if (type_ == E_FACE_TYPE::VTYPE)
      {
        pos += EVec3d(0.0, cvt_props.fold_gap, -(cvt_props.thickness + J_CSG_OFFSET));
        height = cvt_props.fold_gap + J_CSG_OFFSET;
        depth = cvt_props.thickness + 2.0 * J_CSG_OFFSET;
      }
      else if (type_ == E_FACE_TYPE::HTYPE)
      {
        pos += EVec3d(0.0, J_CSG_OFFSET, -J_CSG_OFFSET);
        height = cvt_props.thickness + 2.0 * J_CSG_OFFSET;
        depth = cvt_props.fold_gap + J_CSG_OFFSET;
      }

      JMesh::Rectangular overhang =
      { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };
      cut_box = JCSG::Union(cut_box, JMesh::RectangularMesh(overhang));
    }
  }

  return cut_box;
}


JMesh::Mesh FabricatablePatch::GenerateConcaveCutBox180(const vector<EVec3d>& rect, const ConvertProperties& cvt_props) const
{
  if (mount_face_)
  {
    EVec3d pos = rect[0] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, (cvt_props.fold_gap / 2.0));

    double width = rect[1][0] - rect[0][0] + 2.0 * J_CSG_OFFSET;
    double height = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
    double depth = cvt_props.fold_gap;

    JMesh::Rectangular bend_line_cut_box =
    { pos, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width, height, depth };

    return JMesh::RectangularMesh(bend_line_cut_box);
  }
  else
  {
    EVec3d pos1;
    if (type_ == E_FACE_TYPE::VTYPE)
      pos1 = rect[0] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, cvt_props.fold_gap);
    else if (type_ == E_FACE_TYPE::HTYPE)
      pos1 = rect[0] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, J_CSG_OFFSET);

    double width1 = fabs(rect[1][0] - rect[0][0]) + 2.0 * J_CSG_OFFSET;
    double height1 = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
    double depth1 = cvt_props.fold_gap + J_CSG_OFFSET;

    JMesh::Rectangular concave_cut_box1 =
    { pos1, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width1, height1, depth1 };

    EVec3d pos2;
    if (type_ == E_FACE_TYPE::VTYPE)
      pos2 = rect[3] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, J_CSG_OFFSET);
    else if (type_ == E_FACE_TYPE::HTYPE)
      pos2 = rect[3] - EVec3d(J_CSG_OFFSET, cvt_props.fold_thickness, cvt_props.fold_gap / 2.0);

    double width2 = fabs(rect[1][0] - rect[0][0]) + 2.0 * J_CSG_OFFSET;
    double height2 = cvt_props.thickness - cvt_props.fold_thickness + J_CSG_OFFSET;
    double depth2 = cvt_props.fold_gap / 2.0 + J_CSG_OFFSET;

    JMesh::Rectangular concave_cut_box2 =
    { pos2, EVec3d::UnitX(), -EVec3d::UnitY(), EVec3d::UnitZ(), width2, height2, depth2 };

    return JCSG::Union(JMesh::RectangularMesh(concave_cut_box1), JMesh::RectangularMesh(concave_cut_box2));
  }
}

