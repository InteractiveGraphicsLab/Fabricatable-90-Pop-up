/* jcsg.h */
/* Easy-to-use for libigl mesh_boolean operation */
/*
* Copyright(C) 2023 Junpei Fujikawa (ma22121@shibaura-it.ac.jp)
*
* This program is free software; you can redistribute it and /or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program uses libigl library.
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

#include "jmesh.h"
#include "igl/MeshBooleanType.h"
#include "igl/copyleft/cgal/mesh_boolean.h"
#include "igl/copyleft/cgal/CSGTree.h"


#ifndef J_CSG_OFFSET
#define J_CSG_OFFSET 0.001
#endif // !CSG_OFFSET



namespace JCSG
{
  inline JMesh::Mesh CSG(const JMesh::Mesh& m1, const JMesh::Mesh& m2, const std::string& ope)
  {
    if (ope == "u")
    {
      EMatXd V; EMatXi F;
      igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_UNION, V, F);
      return JMesh::Mesh(V, F);
    }
    else if (ope == "i")
    {
      EMatXd V; EMatXi F;
      igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_INTERSECT, V, F);
      return JMesh::Mesh(V, F);
    }
    else if (ope == "m")
    {
      EMatXd V; EMatXi F;
      igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_MINUS, V, F);
      return JMesh::Mesh(V, F);
    }
    else if (ope == "r")
    {
      EMatXd V; EMatXi F;
      igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_RESOLVE, V, F);
      return JMesh::Mesh(V, F);
    }
    else if (ope == "x")
    {
      EMatXd V; EMatXi F;
      igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_XOR, V, F);
      return JMesh::Mesh(V, F);
    }
    else
    {
      std::cout << "CSG OPERATOR IS NONE!!\n";
      exit(1);
    }
  }


  inline JMesh::Mesh Union(const JMesh::Mesh &m1, const JMesh::Mesh &m2)
  {
    if (!m1.Exist() && m2.Exist())
      return m2;
    else if (m1.Exist() && !m2.Exist())
      return m1;

    EMatXd V; EMatXi F;
    igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_UNION, V, F);
    return JMesh::Mesh(V, F);
  }


  inline JMesh::Mesh Intersect(const JMesh::Mesh &m1, const JMesh::Mesh &m2)
  {
		if (!m1.Exist() && m2.Exist())
			return m2;
		else if (m1.Exist() && !m2.Exist())
			return m1;

    EMatXd V; EMatXi F;
    igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_INTERSECT, V, F);
    return JMesh::Mesh(V, F);
  }


  inline JMesh::Mesh Subtract(const JMesh::Mesh &m1, const JMesh::Mesh &m2)
  {
		if (!m1.Exist() && m2.Exist())
		{
			std::cout << "Error : Minuend Mesh is not Exist!\n";
			exit(1);
		}
		else if (m1.Exist() && !m2.Exist())
		{
			return m1;
		}

    EMatXd V; EMatXi F;
    igl::copyleft::cgal::mesh_boolean(m1.V(), m1.F(), m2.V(), m2.F(), igl::MESH_BOOLEAN_TYPE_MINUS, V, F);
    return JMesh::Mesh(V, F);
  }
}




