/* jmesh.h */
/* Utility functions */
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
#include <queue>


#ifndef J_DBL_EPSILON
#define J_DBL_EPSILON 1e-6
#endif


namespace JUtil
{
	struct Mouse
	{
		bool left;
		bool right;
		bool middle;
	};


	struct Ray
	{
		EVec3d pos;
		EVec3d dir;
	};


	inline Ray GetRay(const OglForQt &ogl, const EVec2i &p)
	{
		EVec3f ray_pos, ray_dir;
		ogl.GetCursorRay(p, ray_pos, ray_dir);

		Ray ray = { EVec3d::Zero(), EVec3d::Zero() };
		ray.pos = ray_pos.cast<double>();
		ray.dir = ray_dir.cast<double>();

		return ray;
	}


	inline std::vector<int> Enumrate(int n)
	{
		std::vector<int> num_array(n);
		for (int i = 0; i < n; i++)
			num_array[i] = i;
		return num_array;
	}


	inline int Findi(const std::vector<int> &list, int val)
	{
		std::vector<int>::const_iterator itr = std::find(list.begin(), list.end(), val);
		if (itr == list.end())
			return -1;

		int idx = std::distance(list.begin(), itr);
		return idx;
	}


	inline std::vector<std::string> Split(const std::string &str, const std::string &sep)
	{
		int sep_len = sep.length();
		std::vector<std::string> split = std::vector<std::string>();

		if (sep_len == 0)
		{
			split.push_back(str);
			return split;
		}

		int offset = 0;
		while (1)
		{
			int sep_pos = str.find(sep, offset);
			if (sep_pos == std::string::npos)
			{
				split.push_back(str.substr(offset));
				break;
			}
			split.push_back(str.substr(offset, sep_pos - offset));
			offset = sep_pos + sep_len;
		}

		return split;
	}


	////////// Error Data Function /////////


	inline EVec2d ErrEVec2d()
	{
		return EVec2d(DBL_MAX, DBL_MAX);
	}


	inline EVec3d ErrEVec3d()
	{
		return EVec3d(DBL_MAX, DBL_MAX, DBL_MAX);
	}


	inline EVec3f ErrEVec3f()
	{
		return EVec3f(FLT_MAX, FLT_MAX, FLT_MAX);
	}


	inline EMatXd ErrEMatXd()
	{
		EMatXd err_mat;
		err_mat.resize(0, 0);
		return err_mat;
	}


	inline EMatXi ErrEMatXi()
	{
		EMatXi err_mat;
		err_mat.resize(0, 0);
		return err_mat;
	}


	////////// ToString Function //////////


	inline std::string ToStringEVec2d(const EVec2d &v)
	{
		std::string x = std::to_string(v[0]);
		std::string y = std::to_string(v[1]);

		return "(" + x + ", " + y + ")";
	}

  inline std::string ToStringEVec2dNonBrackets(const EVec2d& v)
  {
    std::string x = std::to_string(v[0]);
    std::string y = std::to_string(v[1]);

    return x + "," + y;
  }

  inline EVec2d ParseEVec2d(const std::string& v)
  {
    std::vector<std::string> v_str = Split(v, ",");
    double x = std::stod(v_str[0]);
    double y = std::stod(v_str[1]);
    return EVec2d(x, y);
  }

	inline std::string ToStringEVec3d(const EVec3d &v)
	{
		std::string x = std::to_string(v[0]);
		std::string y = std::to_string(v[1]);
		std::string z = std::to_string(v[2]);

		return "(" + x + ", " + y + ", " + z + ")";
	}


	inline std::string ToStringEVec3f(const EVec3f &v)
	{
		std::string x = std::to_string(v[0]);
		std::string y = std::to_string(v[1]);
		std::string z = std::to_string(v[2]);

		return "(" + x + ", " + y + ", " + z + ")";
	}


	inline std::string ToStringStdVectorEVec3d(const std::vector<EVec3d> &vtxs)
	{
		std::string s = "";

		for (int i = 0; i < vtxs.size(); i++)
		{
			s += "[" + std::to_string(i) + "] : (";
			s += std::to_string(vtxs[i][0]) + ", "
				+ std::to_string(vtxs[i][1]) + ", "
				+ std::to_string(vtxs[i][2]);
			s += ")\n";
		}

		return s;
	}


	inline std::string ToStringStdVectorEVec2d(const std::vector<EVec2d> &vtxs)
	{
		std::string s = "";

		for (int i = 0; i < vtxs.size(); i++)
		{
			s += "[" + std::to_string(i) + "] : (";
			s += std::to_string(vtxs[i][0]) + ", "
				+ std::to_string(vtxs[i][1]);
			s += ")\n";
		}

		return s;
	}


	inline std::string ToStringEMatXd(const EMatXd &v)
	{
		std::string s = "";

		for (int i = 0; i < v.rows(); i++)
		{
			s += "[" + std::to_string(i) + "] : (";
			for (int j = 0; j < v.cols(); j++)
			{
				s += std::to_string(v(i, j));
				s += (j != v.cols() - 1) ? ", " : ")\n";
			}
		}

		return s;
	}

	inline std::string ToStringEMatXi(const EMatXi &f)
	{
		std::string s = "";

		for (int i = 0; i < f.rows(); i++)
		{
			s += "[" + std::to_string(i) + "] : (";
			s += std::to_string(f(i, 0)) + ", "
				+ std::to_string(f(i, 1)) + ", "
				+ std::to_string(f(i, 2));
			s += ")\n";
		}

		return s;
	}

	inline std::string ToStringEMatXdEMatXi(const EMatXd &V, const EMatXi &F)
	{
		std::string s = "";

		if (F.rows() > 100)
		{
			for (int i = 0; i < F.rows(); i++)
			{
				if ((i < 10) || (i > F.rows() - 11))
				{
					for (int j = 0; j < F.cols(); j++)
					{
						s += "[" + std::to_string(F(i, j)) + "] : (";
						for (int k = 0; k < V.cols() - 1; k++)
							s += std::to_string(V(F(i, j), k)) + ", ";

						s += std::to_string(V(F(i, j), V.cols() - 1)) + ")\n";
					}
					s += "\n";
				}
				else if (i == F.rows() - 12)
				{
					s += ".\n.\n.\n";
				}
			}
		}
		else
		{
			for (int i = 0; i < F.rows(); i++)
			{
				for (int j = 0; j < F.cols(); j++)
				{
					s += "[" + std::to_string(F(i, j)) + "] : (";
					for (int k = 0; k < V.cols() - 1; k++)
						s += std::to_string(V(F(i, j), k)) + ", ";

					s += std::to_string(V(F(i, j), V.cols() - 1)) + ")\n";
				}
				s += "\n";
			}
		}

		return s;
	}


	////////// Debug Function //////////


	inline void Debug(const EVec2d &v)
	{
		std::cout << ToStringEVec2d(v) << "\n";
	}

	inline void Debug(const EVec3d &v)
	{
		std::cout << ToStringEVec3d(v) << "\n";
	}

	inline void Debug(const EVec3f &v)
	{
		std::cout << ToStringEVec3f(v) << "\n";
	}

	inline void Debug(const std::vector<EVec3d> vtxs)
	{
		std::cout << ToStringStdVectorEVec3d(vtxs) << "\n";
	}

	inline void Debug(const std::vector<EVec2d> vtxs)
	{
		std::cout << ToStringStdVectorEVec2d(vtxs) << "\n";
	}

	inline void Debug(const std::vector<int> v)
	{
		for (int i : v)
			std::cout << i << ", ";
		std::cout << "\n";
	}

	inline void Debug(const EMatXd &v)
	{
		std::cout << ToStringEMatXd(v) << "\n";
	}

	inline void Debug(const EMatXi &f)
	{
		std::cout << ToStringEMatXi(f) << "\n";
	}

	inline void Debug(const EMatXd &V, const EMatXi &F)
	{
		std::cout << ToStringEMatXdEMatXi(V, F) << "\n";
	}


	////////// Cast Function //////////


	// vector<EVec3f> => EMat3f
	inline EMatXf ToEMatXf(const std::vector<EVec3f> &vtxs)
	{
		EMatXf V;
		V.resize(vtxs.size(), 3);

		for (int i = 0; i < vtxs.size(); i++)
		{
			V(i, 0) = vtxs[i][0];
			V(i, 1) = vtxs[i][1];
			V(i, 2) = vtxs[i][2];
		}

		return V;
	}

	// vector<EVec3f> => EMat3d
	inline EMatXd ToEMatXd(const std::vector<EVec3d> &vtxs)
	{
		EMatXd V;
		V.resize(vtxs.size(), 3);

		for (int i = 0; i < vtxs.size(); i++)
		{
			V(i, 0) = vtxs[i][0];
			V(i, 1) = vtxs[i][1];
			V(i, 2) = vtxs[i][2];
		}

		return V;
	}

	// vector<EVec3f> => EMat3d
	inline EMatXd ToEMatXd(const std::vector<EVec3f> &vtxs)
	{
		EMatXd V;
		V.resize(vtxs.size(), 3);

		for (int i = 0; i < vtxs.size(); i++)
		{
			V(i, 0) = (double)(vtxs[i][0]);
			V(i, 1) = (double)(vtxs[i][1]);
			V(i, 2) = (double)(vtxs[i][2]);
		}

		return V;
	}

	// vector<EVec2d> => EMatXd
	inline EMatXd ToEMatXd(const std::vector<EVec2d> &vtxs)
	{
		EMatXd V;
		V.resize(vtxs.size(), 2);

		for (int i = 0; i < vtxs.size(); i++)
		{
			V(i, 0) = (double)(vtxs[i][0]);
			V(i, 1) = (double)(vtxs[i][1]);
		}

		return V;
	}

	// vector<EVec3i> => EMat3i
	inline EMatXi ToEMatXi(const std::vector<EVec3i> &idxs)
	{
		EMatXi F;
		F.resize(idxs.size(), 3);

		for (int i = 0; i < idxs.size(); i++)
		{
			F(i, 0) = idxs[i][0];
			F(i, 1) = idxs[i][1];
			F(i, 2) = idxs[i][2];
		}

		return F;
	}

	inline EMatXi ToEMatXi(std::queue<EVec3i> &idxs)
	{
		EMatXi F(idxs.size(), 3);
		int row = 0;
		while (!idxs.empty())
		{
			EVec3i idx = idxs.front();
			idxs.pop();

			for (int i = 0; i < 3; ++i)
				F(row, i) = idx[i];

			row++;
		}

		return F;
	}


	// EMatXf => vector<EVec3f>
	inline std::vector<EVec3f> ToStdVectorEVec3f(const EMatXf &V)
	{
		std::vector<EVec3f> vtxs = std::vector<EVec3f>(V.rows());

		for (int i = 0; i < V.rows(); i++)
		{
			vtxs[i][0] = V(i, 0);
			vtxs[i][1] = V(i, 1);
			vtxs[i][2] = V(i, 2);
		}

		return vtxs;
	}

	// EMatXd => vector<EVec3f>
	inline std::vector<EVec3f> ToStdVectorEVec3f(const EMatXd &V)
	{
		std::vector<EVec3f> vtxs = std::vector<EVec3f>(V.rows());

		for (int i = 0; i < V.rows(); i++)
		{
			vtxs[i][0] = (float)V(i, 0);
			vtxs[i][1] = (float)V(i, 1);
			vtxs[i][2] = (float)V(i, 2);
		}

		return vtxs;
	}

	// EMatXd => vector<EVec3f>
	inline std::vector<EVec3d> ToStdVectorEVec3d(const EMatXd& V)
	{
		std::vector<EVec3d> vtxs = std::vector<EVec3d>(V.rows());

		for (int i = 0; i < V.rows(); i++)
		{
			vtxs[i][0] = V(i, 0);
			vtxs[i][1] = V(i, 1);
			vtxs[i][2] = V(i, 2);
		}

		return vtxs;
	}

	// EMatXi => vector<EVec3i>
	inline std::vector<EVec3i> ToStdVectorEVec3i(const EMatXi &F)
	{
		std::vector<EVec3i> idxs = std::vector<EVec3i>(F.rows());

		for (int i = 0; i < F.rows(); i++)
		{
			idxs[i][0] = F(i, 0);
			idxs[i][1] = F(i, 1);
			idxs[i][2] = F(i, 2);
		}

		return idxs;
	}


	inline EMat3d ToEMat3d(const EVec3f &v0, const EVec3f &v1, const EVec3f &v2)
	{
		EMat3d mat;
		mat << v0[0], v1[0], v2[0],
			v0[1], v1[1], v2[1],
			v0[2], v1[2], v2[2];
		return mat;
	}


	////////// Equal Function //////////

	inline bool EqualEVec2d(const EVec2d &v1, const EVec2d &v2)
	{
		return (fabs(v1[0] - v2[0]) <= J_DBL_EPSILON) &&
			(fabs(v1[1] - v2[1]) <= J_DBL_EPSILON);
	}


	inline bool EqualEVec3d(const EVec3d &v1, const EVec3d &v2)
	{
		return (fabs(v1[0] - v2[0]) <= J_DBL_EPSILON) &&
			(fabs(v1[1] - v2[1]) <= J_DBL_EPSILON) &&
			(fabs(v1[2] - v2[2]) <= J_DBL_EPSILON);
	}
}
