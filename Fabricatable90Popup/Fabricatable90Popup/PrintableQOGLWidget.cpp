#include "PrintableQOGLWidget.h"

using DMode = E_DrawMode;
using Mode = E_Mode;
using TCore = ToolCore;


PrintableQOGLWidget::PrintableQOGLWidget(QWidget *parent)
	: QOpenGLWidget(parent), m_ogl(this), m_draw_mode(DMode::Printable)
{

}


void PrintableQOGLWidget::initializeGL()
{
	m_ogl.SetBgColor(0.2f, 0.2f, 0.2f, 0.5f);
	m_ogl.SetCam(EVec3f(16.1f, 9.9f, 18.0f),
		           EVec3f( 3.0f, 0.0f,  0.0f),
		           EVec3f( 0.0f, 1.0f,  0.0f));

	// マウスをドラッグしていない間も moseMoveEvent を発生させる
	this->setMouseTracking(true);
}


void PrintableQOGLWidget::resizeGL(int width, int height)
{
	this->repaint();
}

void PrintableQOGLWidget::paintGL()
{
	const int width = this->width();
	const int height = this->height();

	TCore::getInst()->setDrawMode(m_draw_mode);
	m_ogl.OnDrawBegin(width, height);
	ToolCore::getInst()->drawScene(&m_ogl);
	m_ogl.OnDrawEnd();
}


void PrintableQOGLWidget::mousePressEvent(QMouseEvent* mouse_event)
{
	EVec2i p = EVec2i(mouse_event->x(), mouse_event->y());

	switch (mouse_event->button())
	{
	case Qt::LeftButton:
		TCore::getInst()->lMouseButtonPress(&m_ogl, p);
		break;
	case Qt::RightButton:
		TCore::getInst()->rMouseButtonPress(&m_ogl, p);
		break;
	case Qt::MiddleButton:
		TCore::getInst()->mMouseButtonPress(&m_ogl, p);
		break;
	default:
		break;
	}
}

void PrintableQOGLWidget::mouseReleaseEvent(QMouseEvent *mouse_event)
{
	EVec2i p = EVec2i(mouse_event->x(), mouse_event->y());

	switch (mouse_event->button())
	{
	case Qt::LeftButton:
		TCore::getInst()->lMouseButtonRelease(&m_ogl, p);
		break;
	case Qt::RightButton:
		TCore::getInst()->rMouseButtonRelease(&m_ogl, p);
		break;
	case Qt::MiddleButton:
		TCore::getInst()->mMouseButtonRelease(&m_ogl, p);
		break;
	default:
		break;
	}
}

void PrintableQOGLWidget::mouseMoveEvent(QMouseEvent *mouse_event)
{
	EVec2i p = EVec2i(mouse_event->x(), mouse_event->y());
	TCore::getInst()->mouseMove(&m_ogl, p);
	emit manipulatedOGLCamera(m_ogl);
}


void PrintableQOGLWidget::setOGLCamera(const OglForQt &ogl)
{
	m_ogl.SetCam(ogl.GetCamPos(), ogl.GetCamCnt(), ogl.GetCamUp());
	m_ogl.RedrawWindow();
}

void PrintableQOGLWidget::orderToSimulate(int angle)
{
	TCore::getInst()->simulateOA(&m_ogl, angle);
}

void PrintableQOGLWidget::orderToExport()
{
	TCore::getInst()->exportMesh(&m_ogl, "./OA/test1.stl");
}