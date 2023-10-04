#pragma once

#include "E_DrawMode.h"
#include "ToolCore.h"


class PrintableQOGLWidget :
    public QOpenGLWidget
{
  Q_OBJECT
  Q_DISABLE_COPY(PrintableQOGLWidget)


private:
  const E_DrawMode m_draw_mode;

  OglForQt m_ogl;


protected:
  void initializeGL() override;
  void resizeGL(int w, int h) override;
  void paintGL() override;

  void mousePressEvent(QMouseEvent *mouse_event) override;
  void mouseReleaseEvent(QMouseEvent *mouse_event) override;
  void mouseMoveEvent(QMouseEvent *mouse_event) override;


public:
  PrintableQOGLWidget(QWidget *parent = Q_NULLPTR);


signals:
  void manipulatedOGLCamera(const OglForQt &ogl);

public slots:
  void setOGLCamera(const OglForQt &ogl);

public slots:
  void orderToSimulate(int angle);
  void orderToExport();
};

