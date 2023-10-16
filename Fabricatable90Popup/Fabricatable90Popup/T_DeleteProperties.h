/* T_DeformProperties.h */
/* Properties of to be required for Delete operation */
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

#include "E_FoldType.h"
#include "Component.h"


struct DeleteProperties
{
	Component *delete_comp = nullptr;
	Component *parent = nullptr;
	Component *grand_parent = nullptr;
	E_FOLD_TYPE fold_type = E_FOLD_TYPE::DEFAULT;
};
