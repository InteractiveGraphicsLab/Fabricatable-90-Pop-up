/* T_DeformProperties.h */
/* Properties of to be required for Place operation */
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
#include "Component.h"


struct PlaceProperties
{
	Component *parent = nullptr;
	Component *grand_parent = nullptr;
	E_FOLD_TYPE fold_type = E_FOLD_TYPE::DEFAULT;
	EVec3d start_pos = JUtil::ErrEVec3d();
	EVec3d end_pos = JUtil::ErrEVec3d();
	const double radius = 0.1;
	const int slice = 10;
};
