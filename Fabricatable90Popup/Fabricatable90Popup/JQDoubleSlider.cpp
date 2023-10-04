/* JQDoubleSlider.cpp */
/* Double slider classes' function implementation */
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


#include "JQDoubleSlider.h"


JQDoubleSlider::JQDoubleSlider(QWidget *parent)
	: QSlider(parent)
{
	connect(this, SIGNAL(valueChanged(int)), this, SLOT(NoticeValueChanged(int)));
}


void JQDoubleSlider::NoticeValueChanged(int value)
{
	double value_double = (double)value / 1000.0;
	emit ValueChangedDouble(value_double);
}


void JQDoubleSlider::SetValueDouble(double value)
{
	int value_int = (int)(value * 1000);
	setValue(value_int);
	emit ValueChangedDouble(value);
}
