/* ToolCore.cpp */
/* Applications' event manager classes' function implementation */
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


#include "ToolCore.h"


using std::string;
using std::vector;


ToolCore* ToolCore::GetInstance()
{
  static ToolCore core;
  return &core;
}


void ToolCore::PressLeftMouseButton(OglForQt* ogl, const EVec2i& p)
{
  mouse_.left = true;

  JUtil::Ray ray = JUtil::GetRay(*ogl, p);
  if (SelectedMech() && HitDeformHandle(ray))
  {
    prev_ray_ = JUtil::GetRay(*ogl, p);
    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
  }
  else
  {
    ogl->BtnDown_Trans(p);
  }

  ogl->RedrawWindow();
}


void ToolCore::PressRightMouseButton(OglForQt* ogl, const EVec2i& p)
{
  mouse_.right = true;
  ogl->BtnDown_Rot(p);
  ogl->RedrawWindow();
}


void ToolCore::PressMiddleMouseButton(OglForQt* ogl, const EVec2i& p)
{
  mouse_.middle = true;
  ogl->BtnDown_Zoom(p);
  ogl->RedrawWindow();
}


void ToolCore::ReleaseLeftMouseButton(OglForQt* ogl, const EVec2i& p)
{
  JUtil::Ray ray = JUtil::GetRay(*ogl, p);
  switch (edit_mode_)
  {
  case E_EDIT_MODE::PLACEMENT:
    if (Placable())
      Popups::GetInstance()->AddComponent(place_props_);
    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
    ogl->RedrawWindow();
    break;

  case E_EDIT_MODE::DELETION:
    if (Deletable())
      Popups::GetInstance()->DeleteComponent(delete_props_);
    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
    ogl->RedrawWindow();
    break;

  case E_EDIT_MODE::DEFORMATION:
    if (Deformable() && Popups::GetInstance()->HitComponent(ray, deform_props_))
      ProgressDeformStep();
    else if (MovableHandle())
      ProgressDeformStep();
    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
    ogl->RedrawWindow();
    break;

  default:
    InitializeProps();
    break;
  }

  mouse_.left = false;
  ogl->BtnUp();
  ogl->RedrawWindow();
}


void ToolCore::ReleaseRightMouseButton(OglForQt* ogl, const EVec2i& p)
{
  mouse_.right = false;
  ogl->BtnUp();
  ogl->RedrawWindow();
}


void ToolCore::ReleaseMiddleMouseButton(OglForQt* ogl, const EVec2i& p)
{
  mouse_.middle = false;
  ogl->BtnUp();
  ogl->RedrawWindow();
}


void ToolCore::MoveMouse(OglForQt* ogl, const EVec2i& p)
{
  if (!mouse_.left && !mouse_.right && !mouse_.middle)
  {
    JUtil::Ray ray = JUtil::GetRay(*ogl, p);
    switch (edit_mode_)
    {
    case E_EDIT_MODE::PLACEMENT:
      Popups::GetInstance()->HitFoldEdge(ray, place_props_);
      Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
      break;

    case E_EDIT_MODE::DELETION:
      Popups::GetInstance()->HitComponent(ray, delete_props_);
      Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
      break;

    case E_EDIT_MODE::DEFORMATION:
      if (Deformable())
        Popups::GetInstance()->HitComponent(ray, deform_props_);
      else if (SelectedMech())
        HitDeformHandle(ray);
      Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
      break;

    default:
      InitializeProps();
      break;
    }
  }
  else if (mouse_.left && !mouse_.right && !mouse_.middle && MovableHandle())
  {
    JUtil::Ray ray = JUtil::GetRay(*ogl, p);
    EVec3d temp_pos = ray.pos + deform_props_.handle_dist * ray.dir;
    EVec3d prev_pos = prev_ray_.pos + deform_props_.handle_dist * prev_ray_.dir;

    EVec3d move = temp_pos - prev_pos;
    int handle_id = static_cast<int>(deform_props_.handle);
    double move_len = move.dot(deform_props_.handles[handle_id].dir);
    move_len = Popups::GetInstance()->DeformComponent(deform_props_, move_len);
    UpdateDeformHandle(move_len);

    Popups::GetInstance()->UpdateComponent(J_PI / 2.0);
    prev_ray_ = ray;
  }
  else
  {
    ogl->MouseMove(p);
  }

  ogl->RedrawWindow();
}

void ToolCore::DrawScene(OglForQt* ogl, const E_DRAW_MODE& draw_mode, bool draw_mesh_plane)
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, EVec4f(0.3f, 0.3f, 0.3f, 1.0f).data());
  glLightfv(GL_LIGHT0, GL_DIFFUSE, EVec4f(0.8f, 0.8f, 0.8f, 1.0f).data());
  glLightfv(GL_LIGHT0, GL_POSITION, EVec4f(5.0f, 100.0f, 50.0f, 1.0f).data());

  if (draw_mode == E_DRAW_MODE::MODE_2D)
  {
    if (Placable())
    {
      Popups::GetInstance()->Draw2D();
      DrawPlaceFoldEdge();
    }
    else if (Deletable())
    {
      const EVec3f highlight = EVec3f(0.1f, 0.1f, 0.1f);
      Popups::GetInstance()->Draw2D(delete_props_, highlight);
    }
    else if (Deformable())
    {
      if (deform_props_.deform_comp != nullptr)
      {
        const EVec3f highlight = EVec3f(0.1f, 0.1f, 0.1f);
        Popups::GetInstance()->Draw2D(deform_props_, highlight);
      }
      else
      {
        Popups::GetInstance()->Draw2D();
      }
    }
    else if (SelectedMech())
    {
      const EVec3f highlight = EVec3f(0.1f, 0.1f, 0.1f);
      Popups::GetInstance()->Draw2D(deform_props_, highlight);
      DrawDeformHandle();
    }
    else
    {
      Popups::GetInstance()->Draw2D();
    }

    if (imported_mesh_)
      Popups::GetInstance()->DrawMesh(draw_mesh_plane);
  }
  else if (draw_mode == E_DRAW_MODE::MODE_3D)
  {
    Popups::GetInstance()->Draw3D();
    if (imported_mesh_ && !converted_)
      Popups::GetInstance()->DrawMesh(draw_mesh_plane);
  }
}


void ToolCore::EditMode(OglForQt* ogl, const E_EDIT_MODE& edit_mode)
{
  if (converted_ && (edit_mode != E_EDIT_MODE::SIMULATE && edit_mode != E_EDIT_MODE::EXPORT_MESH && edit_mode != E_EDIT_MODE::SAVE_DATA))
  {
    Popups::GetInstance()->ResetPatch();
    converted_ = false;
  }

  InitializeProps();
  edit_mode_ = edit_mode;
  Popups::GetInstance()->UpdateComponent(J_PI / 2.0);
  ogl->RedrawWindow();
  emit SetAngle90(90);
}


void ToolCore::Simulate(OglForQt* ogl, int angle)
{
  InitializeProps();
  edit_mode_ = E_EDIT_MODE::SIMULATE;
  Popups::GetInstance()->UpdateComponent(JMath::Rad(angle));
  ogl->RedrawWindow();
}

void ToolCore::Simulate(OglForQt* ogl, double angle)
{
  InitializeProps();
  edit_mode_ = E_EDIT_MODE::SIMULATE;
  Popups::GetInstance()->UpdateComponent(angle);
  ogl->RedrawWindow();
}


void ToolCore::InitializeDeformMode(OglForQt* ogl)
{
  deform_props_.step = E_DEFORM_STEP::SELECT_MECH;
}


void ToolCore::ImportMesh(OglForQt* ogl, const std::string& file_path)
{
  imported_mesh_ = true;
  string extension = JUtil::Split(file_path, ".").back();

  if (extension == "obj")
  {
    EMatXd V;
    EMatXi F;
    igl::readOBJ(file_path, V, F);
    JMesh::Mesh mesh = JMesh::Mesh(V, F, "c");
    mesh.Move(EVec3d(5.0, -0.6, 0.5));
    Popups::GetInstance()->Mesh(mesh);
  }
  else if (extension == "stl")
  {
    FILE* fp = fopen(file_path.c_str(), "rb");
    EMatXd V;	EMatXi F;	EMatXd N;
    igl::readSTL(fp, V, F, N);
    JMesh::Mesh mesh = JMesh::Mesh(V, F);

    Popups::GetInstance()->Mesh(mesh);
  }

  edit_mode_ = E_EDIT_MODE::DEFAULT;
  ogl->RedrawWindow();
}


void ToolCore::ExportMesh(OglForQt* ogl, const std::string& file_path)
{
  if (converted_)
  {
    InitializeProps();
    edit_mode_ = E_EDIT_MODE::EXPORT_MESH;
    Popups::GetInstance()->UpdateComponent(J_PI);

    std::cout << "\n ----- Start Export Origamic Architecture -----\n\n";
    JMesh::Mesh popups = Popups::GetInstance()->ExportMesh();
    popups.Scale(EVec3d(10.0, 10.0, 10.0));
    igl::writeSTL(file_path, popups.V(), popups.F());
    std::cout << "\n\n ----- Exported Origamic Architecture -----\n";

    Popups::GetInstance()->UpdateComponent(J_PI / 2.0);
  }

  edit_mode_ = E_EDIT_MODE::DEFAULT;
  ogl->RedrawWindow();
}


void ToolCore::LoadData(OglForQt* ogl, const std::string& file_path)
{
  Popups::GetInstance()->UpdateComponent(J_PI / 2.0);

  std::ifstream file(file_path);
  if (file.fail())
  {
    std::cout << "File Load Error" << std::endl;
    return;
  }

  vector<std::string> data;
  data.reserve(10000);
  std::string str;
  while (std::getline(file, str))
    data.push_back(str);

  Popups::GetInstance()->LoadData(data);

  edit_mode_ = E_EDIT_MODE::DEFAULT;
  ogl->RedrawWindow();
}


void ToolCore::SaveData(OglForQt* ogl, const std::string& file_path)
{
  Popups::GetInstance()->UpdateComponent(J_PI / 2.0);

  std::cout << "\n ----- Start Save Origamic Architecture -----\n\n";
  std::ofstream file;
  file.open(file_path, std::ios::out);
  std::string data = Popups::GetInstance()->SaveData();
  file << data << std::endl;
  file.close();
  std::cout << "\n ----- Saved Origamic Architecture -----\n";

  edit_mode_ = E_EDIT_MODE::DEFAULT;
  ogl->RedrawWindow();
}


void ToolCore::Convert(OglForQt* ogl)
{
  if (!Popups::GetInstance()->CheckError())
  {
    converted_ = true;

    Popups::GetInstance()->UpdateComponent(JMath::Rad(180));
    Popups::GetInstance()->FillPlane();

    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
    Popups::GetInstance()->ConvertComponent();

    Popups::GetInstance()->UpdateComponent(JMath::Rad(180));
    Popups::GetInstance()->TrimOutlines();

    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
    bool one_more = Popups::GetInstance()->SegmentSpace();

    if (one_more)
    {
      Popups::GetInstance()->UpdateComponent(JMath::Rad(180));
      Popups::GetInstance()->TrimOutlines();

      Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
      Popups::GetInstance()->SegmentSpace();
    }

    Popups::GetInstance()->UpdateComponent(JMath::Rad(90));
    Popups::GetInstance()->ResetMeshs();
  }
  else
  {
    std::cout << "Detect Error\n" << std::endl;
  }

  edit_mode_ = E_EDIT_MODE::DEFAULT;
  emit SetAngle90(90);

  ogl->RedrawWindow();
}


void ToolCore::SetPlaneThick(OglForQt* ogl, double thickness)
{
  Popups::GetInstance()->PlaneThick(thickness / 10.0);
  ogl->RedrawWindow();
  converted_ = false;
}


void ToolCore::SetPreventBindGap(OglForQt* ogl, double gap)
{
  Popups::GetInstance()->PreventBindGap(gap / 10.0);
  ogl->RedrawWindow();
  converted_ = false;
}


void ToolCore::SetBendLineWidth(OglForQt* ogl, double fold_gap)
{
  Popups::GetInstance()->BendLineWidth(fold_gap / 10.0);
  ogl->RedrawWindow();
  converted_ = false;
}


void ToolCore::SetLaminationPitch(OglForQt* ogl, double fold_thickness)
{
  Popups::GetInstance()->LaminationPitch(fold_thickness / 10.0);
  ogl->RedrawWindow();
  converted_ = false;
}


void ToolCore::SetNozzleWidth(OglForQt* ogl, double nozzle_width)
{
  Popups::GetInstance()->NozzleWidth(nozzle_width / 10.0);
  ogl->RedrawWindow();
  converted_ = false;
}


ToolCore::ToolCore()
  : QObject()
{
  InitializeProps();
  prev_ray_ = { EVec3d::Zero(), EVec3d::Zero() };
  imported_mesh_ = false;
  converted_ = false;
}


void ToolCore::InitializeProps()
{
  mouse_ = { false, false, false };
  edit_mode_ = E_EDIT_MODE::DEFAULT;

  place_props_.parent = nullptr;
  place_props_.grand_parent = nullptr;
  place_props_.fold_type = E_FOLD_TYPE::DEFAULT;
  place_props_.start_pos = JUtil::ErrEVec3d();
  place_props_.end_pos = JUtil::ErrEVec3d();

  delete_props_ = { nullptr, nullptr, nullptr, E_FOLD_TYPE::DEFAULT };

  deform_props_.step = E_DEFORM_STEP::DEFAULT;
  deform_props_.deform_comp = nullptr;
  deform_props_.parent = nullptr;
  deform_props_.grand_parent = nullptr;
  deform_props_.fold_type = E_FOLD_TYPE::DEFAULT;
  deform_props_.handles = vector<JMesh::Arrow>();
  deform_props_.handle = E_DEFORM_HANDLE::DEFAULT;
  deform_props_.handle_dist = -1.0;
}


void ToolCore::ProgressDeformStep()
{
  int step_num = static_cast<int>(deform_props_.step) + 1;

  if (step_num < static_cast<int>(E_DEFORM_STEP::FINISH))
  {
    deform_props_.step = static_cast<E_DEFORM_STEP>(step_num);
  }
  else
  {
    deform_props_.handle = E_DEFORM_HANDLE::DEFAULT;
    deform_props_.handle_dist = -1.0;
    deform_props_.step = E_DEFORM_STEP::SELECT_MECH;
    Popups::GetInstance()->TrimComponent();
  }
}


bool ToolCore::HitDeformHandle(const JUtil::Ray& ray)
{
  double min_dist = DBL_MAX;
  for (int i = 0; i < deform_props_.handles.size(); ++i)
  {
    double dist = DBL_MAX;
    if (JMesh::HitArrow(ray, deform_props_.handles[i], dist))
    {
      if (min_dist > dist)
      {
        min_dist = dist;
        deform_props_.handle = static_cast<E_DEFORM_HANDLE>(i);
      }
    }
  }

  if (min_dist < DBL_MAX)
  {
    deform_props_.handle_dist = min_dist;
    return true;
  }

  deform_props_.handle = E_DEFORM_HANDLE::DEFAULT;
  deform_props_.handle_dist = -1.0;
  return false;
}


void ToolCore::UpdateDeformHandle(double move_len)
{
  int handle_id = static_cast<int>(deform_props_.handle);
  if ((deform_props_.handle == E_DEFORM_HANDLE::LEFT) || (deform_props_.handle == E_DEFORM_HANDLE::RIGHT))
  {
    deform_props_.handles[handle_id].pos += move_len * deform_props_.handles[handle_id].dir;

    deform_props_.handles[1].pos += (move_len / 2.0) * deform_props_.handles[handle_id].dir;
    deform_props_.handles[2].pos += (move_len / 2.0) * deform_props_.handles[handle_id].dir;
  }
  else if ((deform_props_.handle == E_DEFORM_HANDLE::UP) || (deform_props_.handle == E_DEFORM_HANDLE::DEPTH))
  {
    for (JMesh::Arrow& handle : deform_props_.handles)
      handle.pos += move_len * deform_props_.handles[handle_id].dir;
  }
}


void ToolCore::DrawPlaceFoldEdge() const
{
  glDisable(GL_LIGHTING);
  if (Placable())
  {
    const EVec3f color = EVec3f(1.0f, 0.0f, 0.0f);
    JMesh::Cylinder cylinder =
      JMesh::GetCylinder(place_props_.start_pos, place_props_.end_pos, place_props_.radius, place_props_.slice);

    JMesh::CylinderMesh(cylinder).DrawPlaneN(color);
  }
  glEnable(GL_LIGHTING);
}


void ToolCore::DrawDeformHandle() const
{
  glDisable(GL_LIGHTING);

  for (int i = 0; i < deform_props_.handles.size(); ++i)
  {
    EVec3f color = EVec3f::Zero();
    int handle_id = static_cast<int>(deform_props_.handle);
    if (SelectedMech())
      color[i % 3] = (i == handle_id) ? 1.0f : 0.7f;
    else
      color[i % 3] = 0.7f;

    JMesh::ArrowMesh(deform_props_.handles[i]).DrawPlaneN(color);
  }

  glEnable(GL_LIGHTING);
}


bool ToolCore::Placable() const
{
  return (edit_mode_ == E_EDIT_MODE::PLACEMENT) &&
    (place_props_.parent != nullptr) &&
    (place_props_.grand_parent != nullptr) &&
    (place_props_.grand_parent->Id() >= 0) &&
    (place_props_.parent->Id() >= place_props_.grand_parent->Id()) &&
    (place_props_.fold_type != E_FOLD_TYPE::DEFAULT) &&
    (place_props_.start_pos[1] > 0.0 - J_DBL_EPSILON) &&
    (place_props_.end_pos[1] > 0.0 - J_DBL_EPSILON);
}


bool ToolCore::Deletable() const
{
  return (edit_mode_ == E_EDIT_MODE::DELETION) &&
    (delete_props_.delete_comp != nullptr) &&
    (delete_props_.parent != nullptr) &&
    (delete_props_.grand_parent != nullptr) &&
    (delete_props_.grand_parent->Id() >= 0) &&
    (delete_props_.parent->Id() >= delete_props_.grand_parent->Id()) &&
    (delete_props_.delete_comp->Id() > delete_props_.parent->Id()) &&
    (delete_props_.fold_type != E_FOLD_TYPE::DEFAULT);
}


bool ToolCore::Deformable() const
{
  return (edit_mode_ == E_EDIT_MODE::DEFORMATION) &&
    (deform_props_.step == E_DEFORM_STEP::SELECT_MECH);
}


bool ToolCore::SelectedMech() const
{
  return (edit_mode_ == E_EDIT_MODE::DEFORMATION) &&
    (deform_props_.step == E_DEFORM_STEP::MOVE_HANDLE) &&
    (deform_props_.deform_comp != nullptr) &&
    (deform_props_.parent != nullptr) &&
    (deform_props_.grand_parent != nullptr) &&
    (deform_props_.grand_parent->Id() >= 0) &&
    (deform_props_.parent->Id() >= deform_props_.grand_parent->Id()) &&
    (deform_props_.deform_comp->Id() > deform_props_.parent->Id()) &&
    (deform_props_.fold_type != E_FOLD_TYPE::DEFAULT) &&
    (deform_props_.handles.size() == 4);
}


bool ToolCore::MovableHandle() const
{
  int handle_id = static_cast<int>(deform_props_.handle);

  return (edit_mode_ == E_EDIT_MODE::DEFORMATION) &&
    (deform_props_.step == E_DEFORM_STEP::MOVE_HANDLE) &&
    (deform_props_.deform_comp != nullptr) &&
    (deform_props_.parent != nullptr) &&
    (deform_props_.grand_parent != nullptr) &&
    (deform_props_.grand_parent->Id() >= 0) &&
    (deform_props_.parent->Id() >= deform_props_.grand_parent->Id()) &&
    (deform_props_.deform_comp->Id() > deform_props_.parent->Id()) &&
    (deform_props_.fold_type != E_FOLD_TYPE::DEFAULT) &&
    (deform_props_.handles.size() == 4) &&
    (handle_id >= 0) &&
    (deform_props_.handle_dist > 0.0) &&
    (deform_props_.handles[handle_id].pos.nonZeros()) &&
    (deform_props_.handles[handle_id].dir.nonZeros());
}
