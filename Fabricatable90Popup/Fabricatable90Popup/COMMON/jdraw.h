/* jdraw.h */
/* Easy-to-use for glew draw function */
/*
* Copyright(C) 2023 Junpei Fujikawa (ma22121@shibaura-it.ac.jp)
*
* This program is free software; you can redistribute it and /or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program uses glew library.
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
#include "GL/glew.h"


namespace JDraw
{
	inline void SetEVec3fToVertex3f(const EVec3f &v)
	{
		glVertex3f(v[0], v[1], v[2]);
	}


	inline void SetEVec3fToColor3f(const EVec3f &color)
	{
		glColor3f(color[0], color[1], color[2]);
	}


  inline void SetEVec4fToColor4f(const EVec4f& color)
  {
    glColor4f(color[0], color[1], color[2], color[3]);
  }


	inline void SetEVec3fToColor3fv(const EVec3f& color)
	{
		glColor3fv(color.data());
	}


	inline void SetEVec3dToColor3fv(const EVec3d& color)
	{
		SetEVec3fToColor3fv(color.cast<float>());
	}


	inline void SetEVec3fToVertex3fv(const EVec3f& vtx)
	{
		glVertex3fv(vtx.data());
	}


	inline void SetEVec3dToVertex3fv(const EVec3d& vtx)
	{
		SetEVec3fToVertex3fv(vtx.cast<float>());
	}


	inline void SetStdVectorToGLTriangle(const std::vector<EVec3d> &V, const std::vector<EVec3i>& F)
	{
		if (!V.empty() && !F.empty())
		{
			for (int i = 0; i < F.size(); ++i)
			{
				SetEVec3dToVertex3fv(V[F[i][0]]);
				SetEVec3dToVertex3fv(V[F[i][1]]);
				SetEVec3dToVertex3fv(V[F[i][2]]);
			}
		}
	}


	inline void SetStdVectorToGLTriangle(const std::vector<EVec3d> &V, const EMatXi& F)
	{
		if (!V.empty() && F.nonZeros())
		{
      for (int i = 0; i < F.rows(); ++i)
			{
        EVec3d p0 = V[F(i, 0)];
        EVec3d p1 = V[F(i, 1)];
        EVec3d p2 = V[F(i, 2)];

        EVec3f n = -((p1 - p0).cross(p2 - p0)).cast<float>();
        n.normalize();

        glNormal3fv(n.data());
        SetEVec3dToVertex3fv(p0);
				SetEVec3dToVertex3fv(p1);
        SetEVec3dToVertex3fv(p2);
			}
		}
	}


	inline void SetStdVectorToGLLine(const std::vector<EVec3d> &V, const EMatXi &F)
	{
		if (!V.empty() && F.nonZeros())
		{
			for (int i = 0; i < F.rows(); ++i)
			{
				SetEVec3dToVertex3fv(V[F(i, 0)]);
				SetEVec3dToVertex3fv(V[F(i, 1)]);

				SetEVec3dToVertex3fv(V[F(i, 1)]);
				SetEVec3dToVertex3fv(V[F(i, 2)]);

				SetEVec3dToVertex3fv(V[F(i, 2)]);
				SetEVec3dToVertex3fv(V[F(i, 0)]);
			}
		}
	}


	inline void SetEMatXftoGLTriangle(const EMatXf &V, const EMatXi &F, bool alpha=false)
	{
    //if (alpha)
    //{
    //  std::queue<int> que;
    //  for (int i = 0; i < F.rows(); ++i)
    //  {
    //    int idx0 = F(i, 0);
    //    int idx1 = F(i, 1);
    //    int idx2 = F(i, 2);

    //    EVec3f p0 = EVec3f(V(idx0, 0), V(idx0, 1), V(idx0, 2));
    //    EVec3f p1 = EVec3f(V(idx1, 0), V(idx1, 1), V(idx1, 2));
    //    EVec3f p2 = EVec3f(V(idx2, 0), V(idx2, 1), V(idx2, 2));

    //    EVec3f n = ((p1 - p0).cross(p2 - p1)).normalized();

    //    if (fabs(n[0]-1.0f) <= 1e-4f || fabs(n[1] + 1.0f) <= 1e-4f || fabs(n[2] + 1.0f) <= 1e-4f)
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

    //    int idx0 = F(idx, 0);
    //    int idx1 = F(idx, 1);
    //    int idx2 = F(idx, 2);

    //    EVec3f p0 = EVec3f(V(idx0, 0), V(idx0, 1), V(idx0, 2));
    //    EVec3f p1 = EVec3f(V(idx1, 0), V(idx1, 1), V(idx1, 2));
    //    EVec3f p2 = EVec3f(V(idx2, 0), V(idx2, 1), V(idx2, 2));

    //    EVec3f n = ((p1 - p0).cross(p2 - p1)).normalized();

    //    glNormal3fv(n.data());
    //    glVertex3fv(p0.data());
    //    glVertex3fv(p1.data());
    //    glVertex3fv(p2.data());
    //  }

    //  return;
    //}

		if (V.nonZeros() && F.nonZeros())
		{
#pragma omp parallel for
      for (int i = 0; i < F.rows(); ++i)
			{
        int idx0 = F(i, 0);
        int idx1 = F(i, 1);
        int idx2 = F(i, 2);

        EVec3f p0 = EVec3f(V(idx0, 0), V(idx0, 1), V(idx0, 2));
        EVec3f p1 = EVec3f(V(idx1, 0), V(idx1, 1), V(idx1, 2));
        EVec3f p2 = EVec3f(V(idx2, 0), V(idx2, 1), V(idx2, 2));

        EVec3f n = ((p1 - p0).cross(p2 - p1)).normalized();

        glNormal3fv(n.data());
        glVertex3fv(p0.data());
        glVertex3fv(p1.data());
        glVertex3fv(p2.data());
			}
		}
	}

	inline void SetEMatXdtoGLTriangle(const EMatXd &V, const EMatXi &F, bool alpha=false)
	{
		SetEMatXftoGLTriangle(V.cast<float>(), F, alpha);
	}


	inline void SetEMatXftoGLLine(const EMatXf &V, const EMatXi &F)
	{
		if (V.nonZeros() && F.nonZeros())
		{
			for (int i = 0; i < F.rows(); ++i)
			{
				for (int j = 0; j < F.cols(); ++j)
				{

					if (j != 2)
					{
						glVertex3f(V(F(i, j    ), 0), V(F(i, j    ), 1), V(F(i, j    ), 2));
						glVertex3f(V(F(i, j + 1), 0), V(F(i, j + 1), 1), V(F(i, j + 1), 2));
					}
					else
					{
						glVertex3f(V(F(i, j), 0), V(F(i, j), 1), V(F(i, j), 2));
						glVertex3f(V(F(i, 0), 0), V(F(i, 0), 1), V(F(i, 0), 2));
					}
				}
			}
		}
	}

	inline void SetEMatXdtoGLLine(const EMatXd &V, const EMatXi &F)
	{
		SetEMatXftoGLLine(V.cast<float>(), F);
	}
}
