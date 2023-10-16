/* T_DeformProperties.h */
/* Properties of to be required for Resize and Translate operation */
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

#include "jmesh.h"
#include "E_DeformStep.h"
#include "E_FoldType.h"
#include "E_DeformHandle.h"


struct DeformProperties
{
	E_DEFORM_STEP step = E_DEFORM_STEP::DEFAULT;
	Component *deform_comp = nullptr;
	Component *parent = nullptr;
	Component *grand_parent = nullptr;
	E_FOLD_TYPE fold_type = E_FOLD_TYPE::DEFAULT;
	std::vector<JMesh::Arrow> handles = std::vector<JMesh::Arrow>();
	E_DEFORM_HANDLE handle = E_DEFORM_HANDLE::DEFAULT;
	double handle_dist = -1.0;
	const JMesh::ArrowSize handle_size = { 0.2, 1.6, 0.3, 0.4 };
	const int slice = 10;
};
