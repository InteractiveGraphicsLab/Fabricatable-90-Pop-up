/* stdafx.h */
/* Precompile header */
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

#include <iostream>
#include <vector>
#include <queue>
#include <string>

#include "GL/glew.h"
#include <GL/GL.h>
#include <GL/GLU.h>

#include <opencv2/opencv.hpp>

#include <QWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QApplication>
#include <QtGui/QMouseEvent>
#include <QKeyEvent>
#include <QString>
#include <QUrl>
#include <QtWidgets/QFileDialog>
#include <QtOpenGLWidgets/QOpenGLWidget>
#include <QSlider>

#include <igl/MeshBooleanType.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/CSGTree.h>
#include <igl/triangle/triangulate.h>
#include <igl/writeSTL.h>
#include <igl/readOBJ.h>
#include <igl/readSTL.h>

#include "OglForQt.h"
#include "tmath.h"

