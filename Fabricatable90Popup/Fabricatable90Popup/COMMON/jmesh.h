/* jmesh.h */
/* 3D Mesh class */
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
#include "jmath.h"
#include "jdraw.h"
#include "jutil.h"
#include <iostream>
#include <queue>


namespace JMesh
{
	class BoundBox
	{
	private:
		EVec3d m_min_pos;
		EVec3d m_max_pos;

	public:
		BoundBox()
		{
			m_min_pos = EVec3d::Zero();
			m_max_pos = EVec3d::Zero();
		}

		BoundBox(const EMatXd &V) { CalcBoundBox(V); }

		BoundBox(const BoundBox &src)
		{
			this->m_min_pos = src.m_min_pos;
			this->m_max_pos = src.m_max_pos;
		}

		BoundBox &operator=(const BoundBox &src)
		{
			this->m_min_pos = src.m_min_pos;
			this->m_max_pos = src.m_max_pos;

			return *this;
		}

    ~BoundBox() { }

		void CalcBoundBox(const EMatXd &V)
		{
			m_min_pos = EVec3d(DBL_MAX, DBL_MAX, DBL_MAX);
			m_max_pos = EVec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);

			for (int i = 0; i < V.rows(); i++)
			{
				for (int j = 0; j < V.cols(); j++)
				{
					if (m_min_pos[j] > V(i, j))
						m_min_pos[j] = V(i, j);

					if (m_max_pos[j] < V(i, j))
						m_max_pos[j] = V(i, j);
				}
			}
		}

		std::vector<EVec3d> GetBoundBox() const
		{
			std::vector<EVec3d> vtxs = std::vector<EVec3d>
			{
				EVec3d(m_min_pos[0], m_min_pos[1], m_min_pos[2]),
				EVec3d(m_max_pos[0], m_min_pos[1], m_min_pos[2]),
				EVec3d(m_min_pos[0], m_max_pos[1], m_min_pos[2]),
				EVec3d(m_min_pos[0], m_min_pos[1], m_max_pos[2]),
				EVec3d(m_max_pos[0], m_max_pos[1], m_max_pos[2]),
				EVec3d(m_min_pos[0], m_max_pos[1], m_max_pos[2]),
				EVec3d(m_max_pos[0], m_min_pos[1], m_max_pos[2]),
				EVec3d(m_max_pos[0], m_max_pos[1], m_min_pos[2])
			};

			return vtxs;
		}

		EVec3d MinPos() const {	return m_min_pos;	}
		EVec3d MaxPos() const {	return m_max_pos;	}

		void Draw(const EVec3f &color) const
		{
			std::vector<EVec3d> vtxs = GetBoundBox();
			const int LINE = 2;
			const int LINE_NUM = 12;
			const int INDEX[LINE_NUM][LINE] =
			{
				{0, 1}, {0, 2}, {0, 3},
				{2, 5}, {1, 6}, {1, 7},
				{3, 5}, {3, 6}, {2, 7},
				{4, 5}, {4, 6}, {4, 7}
			};

			glBegin(GL_LINES);
			JDraw::SetEVec3fToColor3f(color);
			for (int i = 0; i < LINE_NUM; i++)
				for (int j = 0; j < LINE; j++)
					JDraw::SetEVec3fToVertex3f(vtxs[INDEX[i][j]].cast<float>());
			glEnd();
		}
	};


	class Mesh
	{
	private:
		EMatXd m_V;
		EMatXi m_F;

		BoundBox m_bound_box;

		bool m_exist;

	public:
		Mesh()
		{
			m_V = EMatXd();
			m_F = EMatXi();
			m_bound_box = BoundBox();
			m_exist = false;
		}

		Mesh(const EMatXd &V, const EMatXi &F, const std::string &center = "n")
		{
			m_V = V;
			m_F = F;
			m_bound_box = BoundBox(m_V);
			m_exist = true;
			if (center == "c")
				Centering();
		}

		Mesh(const Mesh &src)
		{
			this->m_V = src.m_V;
			this->m_F = src.m_F;
			this->m_bound_box = src.m_bound_box;
			this->m_exist = src.m_exist;
		}

		Mesh &operator=(const Mesh &src)
		{
			this->m_V = src.m_V;
			this->m_F = src.m_F;
			this->m_bound_box = src.m_bound_box;
			this->m_exist = src.m_exist;

			return *this;
		}

    ~Mesh() { }

		EMatXf Vf()            const { return m_V.cast<float>(); }
    EMatXd V ()            const { return m_V; }
    EVec3d V(int i)        const { return EVec3d(m_V(i, 0), m_V(i, 1), m_V(i, 2));}
    double V(int i, int j) const { return m_V(i, j); }
		void V(int i, const EVec3d &val) { m_V(i, 0) = val[0]; m_V(i, 1) = val[1]; m_V(i, 2) = val[2]; }
		void V(int i, int j, double val) { m_V(i, j) = val; }

    EMatXi F()             const {	return m_F; }
    EVec3i F(int i)        const {	return EVec3i(m_F(i, 0), m_F(i, 1), m_F(i, 2)); }
		int    F(int i, int j) const {	return m_F(i, j);	}
		void F(int i, int j, int val)	{ m_F(i, j) = val; }

		BoundBox GetBoundBox() const { return BoundBox(m_V); }

    bool Exist() const { return m_exist; }

		void Centering()
		{
			float x_center = (m_bound_box.MaxPos()[0] + m_bound_box.MinPos()[0]) / 2.0f;
			float y_offset = m_bound_box.MinPos()[1];
			float z_center = (m_bound_box.MaxPos()[2] + m_bound_box.MinPos()[2]) / 2.0f;

			if (x_center != 0.0f) Move(EVec3d(-x_center, 0.0f, 0.0f));
			if (y_offset != 0.0f) Move(EVec3d(0.0f, -y_offset, 0.0f));
			if (z_center != 0.0f) Move(EVec3d(0.0f, 0.0f, -z_center));

			m_bound_box = BoundBox(m_V);
		}

		void Move(const EVec3d &m)
		{
			for (int i = 0; i < m_V.rows(); i++)
			{
				m_V(i, 0) = m_V(i, 0) + (double)m[0];
				m_V(i, 1) = m_V(i, 1) + (double)m[1];
				m_V(i, 2) = m_V(i, 2) + (double)m[2];
			}

			m_bound_box = BoundBox(m_V);
		}

		void Scale(const EVec3d &s)
		{
			for (int i = 0; i < m_V.rows(); i++)
			{
				m_V(i, 0) = m_V(i, 0) * (double)s[0];
				m_V(i, 1) = m_V(i, 1) * (double)s[1];
				m_V(i, 2) = m_V(i, 2) * (double)s[2];
			}

			m_bound_box = BoundBox(m_V);
		}

    void Rotate(const EVec3d& axis, double angle)
    {
      double cos_t = cos(angle);
      double sin_t = sin(angle);

      EMat3d R;
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);

      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);

      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);
      R(0, 1) = cos_t + axis[0] * axis[0] * (1.0 - cos_t);

    }

		void Compose(const Mesh& mesh)
		{
			if (!mesh.Exist())
				return;

			m_exist = true;

			const int THIS_V_ROW = this->m_V.rows();
			const int THIS_F_ROW = this->m_F.rows();
			const int MESH_V_ROW = mesh.m_V.rows();
			const int MESH_F_ROW = mesh.m_F.rows();

			EMatXd compose_V(THIS_V_ROW + MESH_V_ROW, 3);
			EMatXi compose_F(THIS_F_ROW + MESH_F_ROW, 3);

			compose_V(Eigen::seqN(0, THIS_V_ROW), Eigen::seqN(0, 3)) = this->m_V;
			compose_F(Eigen::seqN(0, THIS_F_ROW), Eigen::seqN(0, 3)) = this->m_F;

			compose_V(Eigen::seqN(THIS_V_ROW, MESH_V_ROW), Eigen::seqN(0, 3)) = mesh.m_V;
			compose_F(Eigen::seqN(THIS_F_ROW, MESH_F_ROW), Eigen::seqN(0, 3)) = mesh.m_F + EMatXi::Constant(MESH_V_ROW, 3, THIS_V_ROW);

			this->m_V = compose_V;
			this->m_F = compose_F;
		}

		void Draw(const EVec3f &color = EVec3f(0.5f, 0.5f, 0.8f)) const
		{
			glBegin(GL_TRIANGLES);
			JDraw::SetEVec3fToColor3f(color);
			JDraw::SetEMatXdtoGLTriangle(m_V, m_F);
			glEnd();

			glBegin(GL_LINES);
			JDraw::SetEVec3fToColor3f(EVec3f(0.0f, 0.0f, 0.0f));
			JDraw::SetEMatXdtoGLLine(m_V, m_F);
			glEnd();
		}

    void DrawPlane(const EVec3f& color = EVec3f(0.5f, 0.5f, 0.8f), bool alpha = false) const
		{
      const EVec4f AMBI(color[0] * 0.8f, color[1] * 0.8f, color[2] * 0.8f, 0.4f);  // 環境光反射成分
      const EVec4f DIFF(color[0], color[1], color[2], 0.6f);  // 拡散反射成分
      const EVec4f SPEC(0.3f, 0.3f, 0.3f, 0.0f);  // 鏡面反射成分
      const float  SHIN[1] = { 64.0f };  // 鏡面反射の強度

      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, DIFF.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, AMBI.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPEC.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, SHIN);

      glBegin(GL_TRIANGLES);
      //if (alpha)
      //{
      //  std::queue<int> que;

      //  for (int i = 0; i < m_F.rows(); ++i)
      //  {
      //    int idx0 = m_F(i, 0);
      //    int idx1 = m_F(i, 1);
      //    int idx2 = m_F(i, 2);

      //    EVec3f p0 = EVec3f(m_V(idx0, 0), m_V(idx0, 1), m_V(idx0, 2));
      //    EVec3f p1 = EVec3f(m_V(idx1, 0), m_V(idx1, 1), m_V(idx1, 2));
      //    EVec3f p2 = EVec3f(m_V(idx2, 0), m_V(idx2, 1), m_V(idx2, 2));

      //    EVec3f n = ((p1 - p0).cross(p2 - p1)).normalized();

      //    if (fabs(n[0] - 1.0f) <= 1e-4f || fabs(n[1] + 1.0f) <= 1e-4f || fabs(n[2] + 1.0f) <= 1e-4f)
      //    {
      //      glNormal3fv(n.data());
      //      glVertex3fv(p0.data());
      //      glVertex3fv(p1.data());
      //      glVertex3fv(p2.data());
      //    }
      //    else
      //    {
      //      que.push(i);
      //    }
      //  }

      //  while (!que.empty())
      //  {
      //    int idx = que.front();
      //    que.pop();

      //    int idx0 = m_F(idx, 0);
      //    int idx1 = m_F(idx, 1);
      //    int idx2 = m_F(idx, 2);

      //    EVec3f p0 = EVec3f(m_V(idx0, 0), m_V(idx0, 1), m_V(idx0, 2));
      //    EVec3f p1 = EVec3f(m_V(idx1, 0), m_V(idx1, 1), m_V(idx1, 2));
      //    EVec3f p2 = EVec3f(m_V(idx2, 0), m_V(idx2, 1), m_V(idx2, 2));

      //    EVec3f n = ((p1 - p0).cross(p2 - p1)).normalized();

      //    glNormal3fv(n.data());
      //    glVertex3fv(p0.data());
      //    glVertex3fv(p1.data());
      //    glVertex3fv(p2.data());
      //  }
      //}
      //else
      //{
        JDraw::SetEMatXdtoGLTriangle(m_V, m_F, alpha);
      //}
      glEnd();
		}

    void DrawPlaneN(const EVec3f& color = EVec3f(0.5f, 0.5f, 0.8f)) const
    {
      glBegin(GL_TRIANGLES);
      JDraw::SetEVec3fToColor3f(color);
      JDraw::SetEMatXdtoGLTriangle(m_V, m_F);
      glEnd();
    }

		void DrawLine(const EVec3f &color = EVec3f(0.0f, 0.0f, 0.0f)) const
		{
			glBegin(GL_LINES);
			JDraw::SetEVec3fToColor3f(color);
			JDraw::SetEMatXdtoGLLine(m_V, m_F);
			glEnd();
		}
	};


	inline Mesh ComposeMesh(const Mesh &mesh1, const Mesh &mesh2)
	{
		const int V1_ROW = mesh1.V().rows();
		const int F1_ROW = mesh1.F().rows();
		const int V2_ROW = mesh2.V().rows();
		const int F2_ROW = mesh2.F().rows();

		EMatXd compose_mesh_V(V1_ROW + V2_ROW, 3);
		EMatXi compose_mesh_F(F1_ROW + F2_ROW, 3);

		compose_mesh_V(Eigen::seqN(0, V1_ROW), Eigen::seqN(0, 3)) = mesh1.V();
		compose_mesh_V(Eigen::seqN(V1_ROW, V2_ROW), Eigen::seqN(0, 3)) = mesh2.V();
		compose_mesh_F(Eigen::seqN(0, F1_ROW), Eigen::seqN(0, 3)) = mesh1.F();
		compose_mesh_F(Eigen::seqN(F1_ROW, F2_ROW), Eigen::seqN(0, 3)) = mesh2.F() + EMatXi::Constant(V2_ROW, 3, V1_ROW);

		return Mesh(compose_mesh_V, compose_mesh_F);
	}


	struct RelSystem
	{
		EVec3d org;
		EVec3d x;
		EVec3d y;
		EVec3d z;
	};


	inline Mesh RelMesh(const Mesh &mesh, const RelSystem &rel_system)
	{
		EMatXd V = mesh.V();

		for (int i = 0; i < mesh.V().rows(); i++)
		{
			EVec3d v = mesh.V(i);
			V(i, 0) = (v - rel_system.org).dot(rel_system.x);
			V(i, 1) = (v - rel_system.org).dot(rel_system.y);
			V(i, 2) = (v - rel_system.org).dot(rel_system.z);
		}

		return Mesh(V, mesh.F());
	}


	inline Mesh StdMesh(const Mesh &mesh, const RelSystem &rel_system)
	{
		EMatXd V = mesh.V();

		for (int i = 0; i < mesh.V().rows(); i++)
		{
			EVec3d v = rel_system.org +
				(mesh.V(i, 0) * rel_system.x) +
				(mesh.V(i, 1) * rel_system.y) +
				(mesh.V(i, 2) * rel_system.z);

			V(i, 0) = v[0];
			V(i, 1) = v[1];
			V(i, 2) = v[2];
		}

		return Mesh(V, mesh.F());
	}


	inline std::vector<EVec3d> ProjectMesh(const Mesh &mesh, const RelSystem &rel_system)
	{
		Mesh rel_mesh = RelMesh(mesh, rel_system);

		std::vector<EVec3d> projection(mesh.V().rows());
		for (int i = 0; i < mesh.V().rows(); i++)
			projection[i] = EVec3d(rel_mesh.V(i, 0), rel_mesh.V(i, 1), 0.0);

		for (int i = 0; i < projection.size(); i++)
			projection[i] = rel_system.org + projection[i][0] * rel_system.x + projection[i][1] * rel_system.y;

		return projection;
	}


	struct Rectangular
	{
		EVec3d pos;
		EVec3d width_dir;
		EVec3d height_dir;
		EVec3d depth_dir;
		double width;
		double height;
		double depth;
	};


	struct Sphere
	{
		EVec3d pos;
		double radius;
		int slice;
	};


	struct Cylinder
	{
		EVec3d top;
		EVec3d dir;
		double radius;
		double length;
		int slice;
	};


	struct Cone
	{
		EVec3d base;
		EVec3d dir;
		double radius;
		double length;
		int slice;
	};


	struct ArrowSize
	{
		double bar_radius;
		double bar_length;
		double tip_radius;
		double tip_length;
	};


	struct Arrow
	{
		EVec3d pos;
		EVec3d dir;
		ArrowSize size;
		int slice;
	};


	struct Circle
	{
		EVec3d pos;
		EVec3d u_dir;
		EVec3d v_dir;
		double radius;
		int slice;
	};


	struct FatLine2D
	{
		EVec3d pos1;
		EVec3d pos2;
		EVec3d u_dir;
		EVec3d v_dir;
		double radius;
		int slice;
	};


	inline Mesh RectangularMesh(const Rectangular &rectangular, const std::string &c = "n")
	{
		EMatXd V(8, 3);
		EMatXi F(12, 3);

		for (int i = 0; i < V.rows(); ++i)
		{
			EVec3d p = rectangular.pos;

			if (i > 3)
				p += rectangular.width * rectangular.width_dir;

			if ((i % 4 == 1) || (i % 4 == 3))
				p += rectangular.depth * rectangular.depth_dir;

			if ((i % 4 == 2) || (i % 4 == 3))
				p += rectangular.height * rectangular.height_dir;

			V(i, 0) = p[0];
			V(i, 1) = p[1];
			V(i, 2) = p[2];
		}

		F(0, 0) = 1; F(0, 1) = 0; F(0, 2) = 3;
		F(1, 0) = 0; F(1, 1) = 2; F(1, 2) = 3;
		F(2, 0) = 2; F(2, 1) = 6; F(2, 2) = 3;
		F(3, 0) = 6; F(3, 1) = 7; F(3, 2) = 3;
		F(4, 0) = 4; F(4, 1) = 0; F(4, 2) = 5;
		F(5, 0) = 5; F(5, 1) = 0; F(5, 2) = 1;
		F(6, 0) = 1; F(6, 1) = 3; F(6, 2) = 5;
		F(7, 0) = 3; F(7, 1) = 7; F(7, 2) = 5;
		F(8, 0) = 7; F(8, 1) = 6; F(8, 2) = 5;
		F(9, 0) = 5; F(9, 1) = 6; F(9, 2) = 4;
		F(10, 0) = 6; F(10, 1) = 2; F(10, 2) = 0;
		F(11, 0) = 4; F(11, 1) = 6; F(11, 2) = 0;

		return Mesh(V, F, c);
	}


	inline Mesh SphereMesh(const Sphere &sphere, const std::string &c = "n")
	{
		int V_row = sphere.slice * sphere.slice * 2 * 3;
		int F_row = sphere.slice * sphere.slice * 2;

		EMatXd V;
		EMatXi F;
		V.resize(V_row, 3);
		F.resize(F_row, 3);

		int V_idx = 0;
		int F_idx = 0;

		for (int i = 0; i < sphere.slice; ++i)
		{
			for (int j = 0; j < sphere.slice; ++j)
			{
				double theta0 = (double)i / (double)sphere.slice * 2.0 * J_PI;
				double theta1 = (double)(i + 1) / (double)sphere.slice * 2.0 * J_PI;
				double phi0 = (double)j / (double)sphere.slice * J_PI - (0.5 * J_PI);
				double phi1 = (double)(j + 1) / (double)sphere.slice * J_PI - (0.5 * J_PI);

				V(V_idx, 0) = sphere.radius * cos(theta0) * cos(phi0);
				V(V_idx, 1) = sphere.radius * sin(theta0) * cos(phi0);
				V(V_idx, 2) = sphere.radius * sin(phi0);
				F(F_idx, 0) = V_idx++;

				V(V_idx, 0) = sphere.radius * cos(theta1) * cos(phi0);
				V(V_idx, 1) = sphere.radius * sin(theta1) * cos(phi0);
				V(V_idx, 2) = sphere.radius * sin(phi0);
				F(F_idx, 1) = V_idx++;

				V(V_idx, 0) = sphere.radius * cos(theta1) * cos(phi1);
				V(V_idx, 1) = sphere.radius * sin(theta1) * cos(phi1);
				V(V_idx, 2) = sphere.radius * sin(phi1);
				F(F_idx, 2) = V_idx++;

				F_idx++;

				V(V_idx, 0) = sphere.radius * cos(theta0) * cos(phi0);
				V(V_idx, 1) = sphere.radius * sin(theta0) * cos(phi0);
				V(V_idx, 2) = sphere.radius * sin(phi0);
				F(F_idx, 0) = V_idx++;

				V(V_idx, 0) = sphere.radius * cos(theta1) * cos(phi1);
				V(V_idx, 1) = sphere.radius * sin(theta1) * cos(phi1);
				V(V_idx, 2) = sphere.radius * sin(phi1);
				F(F_idx, 1) = V_idx++;

				V(V_idx, 0) = sphere.radius * cos(theta0) * cos(phi1);
				V(V_idx, 1) = sphere.radius * sin(theta0) * cos(phi1);
				V(V_idx, 2) = sphere.radius * sin(phi1);
				F(F_idx, 2) = V_idx++;

				F_idx++;
			}
		}

		V = JMath::MoveMat(V, sphere.pos);
		return Mesh(V, F, c);
	}


	inline Cylinder GetCylinder(const EVec3d &top, const EVec3d &base, double radius, int slice)
	{
		return { top, (base - top).normalized(), radius, (base - top).norm(), slice };
	}


	inline Mesh CylinderMesh(const Cylinder &cylinder, const std::string &c = "n")
	{
		const int VERTS_ROW = 2 * (cylinder.slice + 1);
		const int VERTS_COL = 3;
		const int IDXES_ROW = 2 * (cylinder.slice * 2);
		const int IDXES_COL = 3;

		EMatXd V;
		EMatXi F;
		V.resize(VERTS_ROW, VERTS_COL);
		F.resize(IDXES_ROW, IDXES_COL);

		V(0, 0) = 0.0, V(0, 1) = 0.0, V(0, 2) = 0.0;

		for (int i = 1; i < (VERTS_ROW / 2); ++i)
		{
			double t = 2.0 * (double)J_PI / cylinder.slice * (double)i;
			V(i, 0) = cylinder.radius * cos(t);
			V(i, 1) = 0.0;
			V(i, 2) = cylinder.radius * sin(t);
		}

		V(VERTS_ROW / 2, 0) = 0.0, V(VERTS_ROW / 2, 1) = -cylinder.length, V(VERTS_ROW / 2, 2) = 0.0;

		for (int i = (VERTS_ROW / 2) + 1; i < VERTS_ROW; ++i)
		{
			double t = 2.0 * (double)J_PI / cylinder.slice * (double)i;
			V(i, 0) = cylinder.radius * cos(t);
			V(i, 1) = -cylinder.length;
			V(i, 2) = cylinder.radius * sin(t);
		}

		for (int i = 0; i < IDXES_ROW; ++i)
		{
			if ((0 <= i) && (i < cylinder.slice))
			{
				F(i, 0) = i + 1;
				F(i, 1) = 0;
				F(i, 2) = (i == (cylinder.slice - 1) ? 1 : (i + 2));
			}
			else if ((cylinder.slice <= i) && (i < (IDXES_ROW - cylinder.slice)))
			{
				int j = (int)floor((double)(i - cylinder.slice) / 2.0);

				F(i, 0) = j + 1;
				F(i, 1) = (i == ((IDXES_ROW - cylinder.slice) - 2)) ? 1 : (j + 2);
				F(i, 2) = j + 2 + cylinder.slice;

				i++;

				F(i, 0) = (i == ((IDXES_ROW - cylinder.slice) - 1) ? (cylinder.slice + 2) : (j + 3 + cylinder.slice));
				F(i, 1) = F(i - 1, 2);
				F(i, 2) = F(i - 1, 1);
			}
			else
			{
				int k = i - (IDXES_ROW - cylinder.slice);

				F(i, 0) = (i == IDXES_ROW - 1) ? (cylinder.slice + 2) : (cylinder.slice + 3 + k);
				F(i, 1) = VERTS_ROW / 2;
				F(i, 2) = cylinder.slice + 2 + k;
			}
		}

		EVec3d rotate_axis = cylinder.dir.cross(EVec3d::UnitY()).normalized();
		double angle = JMath::CalcSubtendAngle(cylinder.dir, EVec3d::UnitY());
		V = JMath::RotateMatAroundAxis(V, rotate_axis, J_PI - angle);
		V = JMath::MoveMat(V, cylinder.top);

		return Mesh(V, F, c);
	}


	inline Mesh ConeMesh(const Cone &cone, const std::string &c = "n")
	{
		const int VERTS_ROW = cone.slice + 2;
		const int VERTS_COL = 3;
		const int IDXES_ROW = cone.slice * 2;
		const int IDXES_COL = 3;

		EMatXd V;
		EMatXi F;
		V.resize(VERTS_ROW, VERTS_COL);
		F.resize(IDXES_ROW, IDXES_COL);

		V(0, 0) = 0.0, V(0, 1) = 0.0, V(0, 2) = 0.0;

		for (int i = 1; i < (VERTS_ROW - 1); ++i)
		{
			double t = 2.0 * J_PI / (double)cone.slice * (double)i;
			V(i, 0) = cone.radius * cos(t);
			V(i, 1) = 0.0;
			V(i, 2) = cone.radius * sin(t);
		}

		V(VERTS_ROW - 1, 0) = 0.0, V(VERTS_ROW - 1, 1) = -cone.length, V(VERTS_ROW - 1, 2) = 0.0;

		for (int i = 0; i < IDXES_ROW; ++i)
		{
			if (0 <= i && i < cone.slice)
			{
				F(i, 0) = i + 1;
				F(i, 1) = 0;
				F(i, 2) = ((i + 1) == cone.slice) ? 1 : (i + 2);
			}
			else
			{
				F(i, 0) = i + 1 - cone.slice;
				F(i, 1) = cone.slice + 1;
				F(i, 2) = ((i + 1 - cone.slice) == VERTS_ROW - 2) ? 1 : (i + 2 - cone.slice);
			}
		}

		EVec3d rotate_axis = cone.dir.cross(EVec3d::UnitY()).normalized();
		double angle = JMath::CalcSubtendAngle(cone.dir, EVec3d::UnitY());
		V = JMath::RotateMatAroundAxis(V, rotate_axis, J_PI - angle);
		V = JMath::MoveMat(V, cone.base);

		return Mesh(V, F, c);
	}


	inline Mesh ArrowMesh(const Arrow &arrow, const std::string &c = "n")
	{
		EVec3d base = arrow.pos;
		EVec3d top = base + (arrow.size.bar_length * arrow.dir);

		Cylinder bar = GetCylinder(base, top, arrow.size.bar_radius, arrow.slice);
		Cone tip = { top, arrow.dir, arrow.size.tip_radius, arrow.size.tip_length, arrow.slice };

		return ComposeMesh(CylinderMesh(bar, c), ConeMesh(tip, c));
	}


	inline Mesh CircleMesh(const Circle& circle, const std::string& c = "n")
	{
		const int V_ROW = circle.slice + 1;
		const int V_COL = 3;
		const int F_ROW = circle.slice;
		const int F_COL = 3;

		EMatXd V(V_ROW, V_COL);
		EMatXi F(F_ROW, F_COL);

		for (int j = 0; j < V_COL; j++)
			V(0, j) = circle.pos[j];

		for (int i = 1; i < V_ROW; i++)
		{
			double t = 2.0 * (double)J_PI / circle.slice * (double)i;
			double u = circle.radius * cos(t);
			double v = circle.radius * sin(t);

			EVec3d p = circle.pos + (u * circle.u_dir) + (v * circle.v_dir);

			for (int j = 0; j < V_COL; j++)
				V(i, j) = p[j];
		}

		for (int i = 0; i < F_ROW; i++)
		{
			F(i, 0) = i + 1;
			F(i, 1) = 0;
			F(i, 2) = (i == (F_ROW - 1) ? 1 : (i + 2));
		}

		return Mesh(V, F, c);
	}


	inline Mesh FatLine2DMesh(const FatLine2D& fat_line_2d, const std::string& c = "n")
	{
		const int V_ROW = fat_line_2d.slice + 4;
		const int V_COL = 3;
		const int F_ROW = fat_line_2d.slice + 4;
		const int F_COL = 3;

		EMatXd V(V_ROW, V_COL);
		EMatXi F(F_ROW, F_COL);

		for (int j = 0; j < V_COL; j++)
			V(0, j) = fat_line_2d.pos1[j];

		for (int i = 1; i < (V_ROW / 2); i++)
		{
			double t = 2.0 * (double)J_PI / (double)fat_line_2d.slice * (double)(i - 1);
			double u = fat_line_2d.radius * cos(t);
			double v = fat_line_2d.radius * sin(t);

			EVec3d p = fat_line_2d.pos1 + (u * fat_line_2d.u_dir) - (v * fat_line_2d.v_dir);

			for (int j = 0; j < V_COL; j++)
				V(i, j) = p[j];
		}

		for (int j = 0; j < V_COL; j++)
			V((V_ROW / 2), j) = fat_line_2d.pos2[j];

		for (int i = (V_ROW / 2) + 1; i < V_ROW; i++)
		{
			double t = 2.0 * (double)J_PI / (double)fat_line_2d.slice * (double)(i - (V_ROW / 2) - 1);
			double u = fat_line_2d.radius * cos(t);
			double v = fat_line_2d.radius * sin(t);

			EVec3d p = fat_line_2d.pos2 + (u * fat_line_2d.u_dir) + (v * fat_line_2d.v_dir);

			for (int j = 0; j < V_COL; j++)
				V(i, j) = p[j];
		}

		int i = 0;
		for (i = 0; i < (F_ROW - 4) / 2; i++)
		{
			F(i, 0) = i + 1;
			F(i, 1) = 0;
			F(i, 2) = i + 2;
		}

		i = (F_ROW - 4) / 2;
		F(i, 0) = 1;
		F(i, 1) = V_ROW / 2 + 1;
		F(i, 2) = V_ROW / 2;

		i++;
		F(i, 0) = 1;
		F(i, 1) = V_ROW / 2;
		F(i, 2) = 0;

		i++;
		F(i, 0) = 0;
		F(i, 1) = V_ROW / 2;
		F(i, 2) = V_ROW - 1;

		i++;
		F(i, 0) = 0;
		F(i, 1) = V_ROW - 1;
		F(i, 2) = V_ROW / 2 - 1;

		int j = ++i;
		for (; i < F_ROW; i++)
		{
			F(i, 0) = V_ROW / 2 + 1 + (i - j);
			F(i, 1) = V_ROW / 2;
			F(i, 2) = V_ROW / 2 + 2 + (i - j);
		}

		return Mesh(V, F);
	}


	/////////////////////////////////////////////////////////////////////////////////////


	inline bool HitSphere(const JUtil::Ray &ray, const Sphere &sphere, double &dist)
	{
		double a = ray.dir.squaredNorm();
		double b = 2.0 * ((ray.pos - sphere.pos).dot(ray.dir));
		double c = (ray.pos - sphere.pos).squaredNorm() - pow(sphere.radius, 2);

		double t1 = (-b + sqrt(pow(b, 2) - (4.0 * a * c))) / (2.0 * a);
		double t2 = (-b - sqrt(pow(b, 2) - (4.0 * a * c))) / (2.0 * a);

		dist = t1 < t2 ? t1 : t2;

		if (dist >= 0.0)
			return true;
		else
			return false;
	}


	inline bool HitCylinder(const JUtil::Ray &ray, const Cylinder &cylinder, double &dist)
	{
		EVec3d top = cylinder.top;
		EVec3d base = cylinder.top + cylinder.length * cylinder.dir;
		EVec3d cylinder_dir = (top - base).normalized();
		double cylinder_len = (top - base).norm();

		double d = ray.dir.dot(cylinder_dir);
		EVec3d v = base - ray.pos;

		double s = ((d * double(v.dot(ray.dir))) - double(v.dot(cylinder_dir))) / (1.0 - (d * d));
		double t = (double(v.dot(ray.dir)) - (d * double(v.dot(cylinder_dir)))) / (1.0 - (d * d));

		EVec3d ray_vector = ray.pos + t * ray.dir;
		EVec3d cylinder_vector = base + s * cylinder_dir;
		EVec3d normal = cylinder_vector - ray_vector;

		if ((0 <= s && s <= cylinder_len) && (normal.norm() < cylinder.radius))
		{
			dist = t;
			return true;
		}

		dist = DBL_MAX;
		return false;
	}


	inline bool HitArrow(const JUtil::Ray &ray, const Arrow &arrow, double &dist)
	{
		EVec3d base = arrow.pos;
		EVec3d top = base + ((arrow.size.bar_length + arrow.size.tip_length) * arrow.dir);

		Cylinder cylinder = GetCylinder(base, top, arrow.size.bar_radius, arrow.slice);

		return HitCylinder(ray, cylinder, dist);
	}


  inline EMat3d ToEMat3d(const EVec3d& a, const EVec3d& b, const EVec3d& c)
  {
    EMat3d A;
    A << a[0], b[0], c[0],
         a[1], b[1], c[1],
         a[2], b[2], c[2];

    return A;
  }


	inline bool HitPolygon(const JUtil::Ray &ray, const EMat3d &V, double &dist)
	{
		constexpr double EPSILON = 1e-6;

		EVec3d e1 = EVec3d(V(1, 0) - V(0, 0), V(1, 1) - V(0, 1), V(1, 2) - V(0, 2));
		EVec3d e2 = EVec3d(V(2, 0) - V(0, 0), V(2, 1) - V(0, 1), V(2, 2) - V(0, 2));
		EVec3d r  = EVec3d(ray.pos[0] - V(0, 0),	ray.pos[1] - V(0, 1),	ray.pos[2] - V(0, 2));

    EMat3d A;
    A << e1[0], e2[0], -ray.dir[0],
         e1[1], e2[1], -ray.dir[1],
         e1[2], e2[2], -ray.dir[2];

    double detA = A.determinant();

    if (detA == 0.0)
      return false;

    EMat3d B = ToEMat3d(r, e2, -ray.dir);
    EMat3d C = ToEMat3d(e1, r, -ray.dir);
    EMat3d D = ToEMat3d(e1, e2, r);

    double u = B.determinant() / detA;
    double v = C.determinant() / detA;

    if (u >= 0.0 && v >= 0.0 && (u + v <= 1.0))
    {
      double t = D.determinant() / detA;
      if (0.0 < t && t < dist)
      {
        dist = t;
        return true;
      }
    }

    return false;
	}


	inline bool HitPlane(const JUtil::Ray &ray,	const EMatXd &V, const EMatXi &F, double &dist)
	{
		double min_dist = dist;
#pragma omp parallel for
		for (int i = 0; i < F.rows(); i++)
		{
			EMat3d v;
			v << V(F(i, 0), 0), V(F(i, 0), 1), V(F(i, 0), 2),
				V(F(i, 1), 0), V(F(i, 1), 1), V(F(i, 1), 2),
				V(F(i, 2), 0), V(F(i, 2), 1), V(F(i, 2), 2);

			double d = DBL_MAX;
			if (HitPolygon(ray, v, d))
			{
				if (d < min_dist)
				{
					min_dist = d;
					break;
				}
			}
		}

		if (min_dist < dist)
		{
			dist = min_dist;
			return true;
		}

		return false;
	}


	inline EMatXi Triangulation2D(const EMatXd &V, bool inverse_face)
	{
		std::queue<EVec3i> idx_que;

		std::vector<int> not_visit = JUtil::Enumrate(V.rows());
		EVec3i triangle = EVec3i(-1, -1, -1);
		while (not_visit.size() > 3)
		{
			if (triangle[0] == -1)
				triangle[0] = JMath::SearchFarPosIdxInRange2D(V, not_visit, EVec2d(0.0, 0.0));

			int far_find = JUtil::Findi(not_visit, triangle[0]);

			if (inverse_face)
			{
				if (triangle[1] == -1)
					triangle[1] = (triangle[0] == not_visit[not_visit.size() - 1]) ? (not_visit[0]) : (not_visit[far_find + 1]);
				if (triangle[2] == -1)
					triangle[2] = (triangle[0] == not_visit[0]) ? (not_visit[not_visit.size() - 1]) : (not_visit[far_find - 1]);
			}
			else
			{
				if (triangle[1] == -1)
					triangle[1] = (triangle[0] == not_visit[0]) ? (not_visit[not_visit.size() - 1]) : (not_visit[far_find - 1]);
				if (triangle[2] == -1)
					triangle[2] = (triangle[0] == not_visit[not_visit.size() - 1]) ? (not_visit[0]) : (not_visit[far_find + 1]);
			}

			if (JMath::InsideTriangle2D(V, triangle) ||
					(JMath::FaceBackwardTriangle2D(V, triangle) && inverse_face) ||
					(JMath::FaceForwardTriangle2D(V, triangle) && !inverse_face))
			{
				triangle = EVec3i(triangle[1], -1, triangle[0]);
				continue;
			}

			idx_que.push(triangle);
			not_visit.erase(not_visit.begin() + far_find);
			triangle = EVec3i(-1, -1, -1);
		}

		triangle[0] = JMath::SearchFarPosIdxInRange2D(V, not_visit, EVec2d(0.0, 0.0));
		int far_find = JUtil::Findi(not_visit, triangle[0]);
		if (inverse_face)
		{
			if (triangle[1] == -1)
				triangle[1] = (triangle[0] == not_visit[not_visit.size() - 1]) ? (not_visit[0]) : (not_visit[far_find + 1]);
			if (triangle[2] == -1)
				triangle[2] = (triangle[0] == not_visit[0]) ? (not_visit[not_visit.size() - 1]) : (not_visit[far_find - 1]);
		}
		else
		{
			if (triangle[1] == -1)
				triangle[1] = (triangle[0] == not_visit[0]) ? (not_visit[not_visit.size() - 1]) : (not_visit[far_find - 1]);
			if (triangle[2] == -1)
				triangle[2] = (triangle[0] == not_visit[not_visit.size() - 1]) ? (not_visit[0]) : (not_visit[far_find + 1]);
		}
		idx_que.push(triangle);

		return JUtil::ToEMatXi(idx_que);
	}


	inline std::string ToStringMesh(const Mesh &mesh)	{ return JUtil::ToStringEMatXdEMatXi(mesh.V(), mesh.F());	}

	inline std::string ToStringBoundBox(const BoundBox &bound_box)
	{
		std::string s = "Min Pos : " + JUtil::ToStringEVec3d(bound_box.MinPos()) + "\n"
			+ "Max Pos : " + JUtil::ToStringEVec3d(bound_box.MaxPos()) + "\n------\n";

		for (const EVec3d &v : bound_box.GetBoundBox())
			s += JUtil::ToStringEVec3d(v) + "\n";

		return s;
	}

	inline void Debug(const Mesh &mesh)
	{
		std::cout << "Mesh V:" << mesh.V().rows() << " F:" << mesh.F().rows() << "\n-----\n";
		std::cout << ToStringMesh(mesh) << "\n";
	}

	inline void Debug(const BoundBox &bound_box)
	{
		std::cout << ToStringBoundBox(bound_box) << "\n";
	}

}











