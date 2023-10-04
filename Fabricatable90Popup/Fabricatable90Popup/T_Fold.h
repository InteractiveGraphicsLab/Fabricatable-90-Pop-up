/* T_Fold.h */
/* Fold properties */
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

#include "jutil.h"
#include "E_FoldType.h"


struct Fold
{
	E_FOLD_TYPE type = E_FOLD_TYPE::DEFAULT;
	EVec2d rel_org = JUtil::ErrEVec2d();
  EVec3d org = JUtil::ErrEVec3d();
  EVec3d org90 = JUtil::ErrEVec3d();
	double u = -1.0;
	double v = -1.0;
	double w = -1.0;
};
