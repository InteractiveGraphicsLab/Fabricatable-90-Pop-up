/* E_EditMode.h */
/* Operation mode name class file */
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


enum class E_EDIT_MODE : int
{
	DEFAULT = -1,
	SIMULATE = 0,
	PLACEMENT = 1,
	DELETION = 2,
	DEFORMATION = 3,
	IMPORT_MESH = 4,
	EXPORT_MESH = 5,
	CONVERSION = 6,
  LOAD_DATA = 7,
  SAVE_DATA = 8
};
