/* Popups.cpp */
/* Pop-up Mechanism classes' function implementation */
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

#include "Popups.h"


using std::vector;
using std::pair;
using std::tuple;
using std::queue;
using std::stack;


Popups* Popups::GetInstance()
{
  static Popups oa;
  return &oa;
}


void Popups::Draw2D() const
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->Draw2D(angle_);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}


void Popups::Draw2D(const DeleteProperties& delete_props, const EVec3f& highlight) const
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    if (comp->Id() == delete_props.delete_comp->Id())
      comp->Draw2D(angle_, highlight);
    else
      comp->Draw2D(angle_);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}


void Popups::Draw2D(const DeformProperties& deform_props, const EVec3f& highlight) const
{
  DeleteProperties props = { deform_props.deform_comp,
    deform_props.parent, deform_props.grand_parent, deform_props.fold_type };

  Draw2D(props, highlight);
}


void Popups::Draw3D() const
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->Draw3D(angle_, cvt_props_);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}


void Popups::DrawMesh(bool draw_mesh_plane) const
{
  if (draw_mesh_plane)
    mesh_.DrawPlane();
}


void Popups::UpdateComponent(double angle)
{
  angle_ = angle;

  queue<pair<Component*, Fold>> que;
  Fold fold = { E_FOLD_TYPE::DEFAULT, EVec2d::Zero(), EVec3d::Zero(), EVec3d::Zero(), 0.0, 0.0, 0.0 };
  que.push({ root_, fold });

  while (!que.empty())
  {
    Component* comp = que.front().first;
    Fold depend_fold = que.front().second;
    que.pop();

    comp->Update(depend_fold, angle_);

    for (int i = 0; i < comp->Children().size(); i++)
    {
      for (int j = 0; j < comp->Children(i).size(); j++)
        que.push({ comp->Child(i, j), comp->FoldEdge(i) });
    }
  }
}


void Popups::AddComponent(const PlaceProperties& place_props)
{
  const EVec2d scale = EVec2d(0.8, 0.35);
  Fold place_fold = place_props.parent->FoldEdge(place_props.fold_type);

  place_fold.rel_org[0] = place_props.start_pos[0] - place_fold.org[0];
  place_fold.org = place_props.start_pos;
  place_fold.u = (place_props.start_pos - place_props.end_pos).norm();

  Component* child = new Component(place_fold, scale);
  place_props.parent->AddChild(place_props.fold_type, child);

  TrimComponent();
}


void Popups::DeleteComponent(const DeleteProperties& delete_props)
{
  delete_props.parent->DeleteChild(delete_props.fold_type, delete_props.delete_comp->Id());

  queue<Component*> que;
  que.push(delete_props.delete_comp);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }

    delete comp;
  }

  TrimComponent();
}


double Popups::DeformComponent(const DeformProperties& deform_props, double move_len)
{
  Fold prev_v_fold = deform_props.deform_comp->FoldEdge(E_FOLD_TYPE::VTYPE);
  Fold prev_h_fold = deform_props.deform_comp->FoldEdge(E_FOLD_TYPE::HTYPE);
  move_len = deform_props.deform_comp->DeformPatch(deform_props.handle, move_len, cvt_props_);

  return move_len;
}


bool Popups::CheckError()
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    if (comp->Err())
      return true;

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }

  return false;
}


void Popups::FillPlane()
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->FillPlane(cvt_props_);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}


void Popups::ConvertComponent()
{
  queue<pair<Component*, Fold>> que;
  Fold default_fold = { E_FOLD_TYPE::DEFAULT, EVec2d::Zero(), EVec3d::Zero(), EVec3d::Zero(), 0.0, 0.0, 0.0 };
  que.push({ root_, default_fold });

  while (!que.empty())
  {
    Component* comp = que.front().first;
    Fold depend_fold = que.front().second;
    que.pop();

    comp->CutPrintableOAFace(mesh_, depend_fold, cvt_props_);

    for (int i = 0; i < comp->Children().size(); i++)
    {
      for (int j = 0; j < comp->Children(i).size(); j++)
        que.push({ comp->Child(i, j), comp->FoldEdge(i) });
    }
  }
}


static JMesh::Mesh sGenerateOutline(Component* comp, const E_FACE_TYPE& type, const ConvertProperties& cvt_props)
{

  int i = (comp->Id() != 0) ? (static_cast<int>(type)) : (static_cast<int>(E_FOLD_TYPE::GTYPE));

  queue<Component*> que;
  for (int j = 0; j < (comp->Children(i)).size(); j++)
    que.push(comp->Child(i, j));

  JMesh::Mesh outlines = JMesh::Mesh();
  while (!que.empty())
  {
    Component* child = que.front();
    que.pop();

    outlines = JCSG::Union(child->InflatedPatch(), outlines);
    for (vector<Component*> grand_children : child->Children())
    {
      for (Component* grand_child : grand_children)
        que.push(grand_child);
    }
  }

  return outlines;
}


void Popups::TrimOutlines()
{
  queue<Component*> que;
  stack<Component*> stk;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    if (comp->IsExistChilren())
    {
      stk.push(comp);
    }
    else
    {
      for (int i = 0; i < static_cast<int>(E_FACE_TYPE::NUM_OF_TYPES); i++)
      {
        E_FACE_TYPE type = static_cast<E_FACE_TYPE>(i);
        JMesh::Mesh outlines = JMesh::Mesh();
        comp->TrimOutlines(type, outlines, cvt_props_);
      }

      comp->GenerateInflatedPatch(cvt_props_);
    }

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }

  while (!stk.empty())
  {
    Component* comp = stk.top();
    stk.pop();

    for (int i = 0; i < static_cast<int>(E_FACE_TYPE::NUM_OF_TYPES); i++)
    {
      E_FACE_TYPE type = static_cast<E_FACE_TYPE>(i);
      JMesh::Mesh outlines = sGenerateOutline(comp, type, cvt_props_);
      comp->TrimOutlines(type, outlines, cvt_props_);
    }

    comp->GenerateInflatedPatch(cvt_props_);
  }
}


static vector<Component*> sCollectMech(Component* root)
{
  queue<Component*> que;
  que.push(root);

  vector<Component*> comp_list(0);
  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp_list.push_back(comp);

    for (int i = 0; i < comp->Children().size(); i++)
    {
      for (int j = 0; j < comp->Children(i).size(); j++)
        que.push(comp->Child(i, j));
    }
  }

  return comp_list;
}

static vector<Component*> sSortMech(const vector<Component*> comp_list, const E_FACE_TYPE& type)
{
  vector<Component*> sort_list = comp_list;
  int dim = (type == E_FACE_TYPE::HTYPE) ? 1 : 2;
  for (int k = 1; k < sort_list.size(); ++k)
  {
    for (int i = 0; i < sort_list.size() - k; ++i)
    {
      if (sort_list[i]->Rect(type).p[0][dim] < sort_list[i + 1]->Rect(type).p[0][dim])
      {
        Component* tmp = sort_list[i];
        sort_list[i] = sort_list[i + 1];
        sort_list[i + 1] = tmp;
      }
    }
  }

  return sort_list;
}


bool Popups::SegmentSpace()
{
  UpdateComponent(J_PI / 2.0);

  vector<Component*> comp_list = sCollectMech(root_);

  vector<Component*> v_sort_list = sSortMech(comp_list, E_FACE_TYPE::VTYPE);
  vector<Component*> h_sort_list = sSortMech(comp_list, E_FACE_TYPE::HTYPE);

  Fold root_fold = root_->FoldEdge(E_FOLD_TYPE::GTYPE);
  EVec3d max_pos = EVec3d(root_fold.u, root_fold.v, root_fold.w);
  JMesh::Mesh higher_space = JMesh::Mesh();

  bool one_more = false;
  for (int i = 0; i < v_sort_list.size(); i++)
  {
    if (one_more)
      v_sort_list[i]->Segment(higher_space, mesh_, E_FACE_TYPE::VTYPE, cvt_props_, max_pos);
    else
      one_more = v_sort_list[i]->Segment(higher_space, mesh_, E_FACE_TYPE::VTYPE, cvt_props_, max_pos);
  }

  higher_space = JMesh::Mesh();
  for (int i = 0; i < h_sort_list.size(); i++)
  {
    if (one_more)
      h_sort_list[i]->Segment(higher_space, mesh_, E_FACE_TYPE::HTYPE, cvt_props_, max_pos);
    else
      one_more = h_sort_list[i]->Segment(higher_space, mesh_, E_FACE_TYPE::HTYPE, cvt_props_, max_pos);
  }

  return one_more;
}



JMesh::Mesh Popups::ExportMesh()
{
  JMesh::Mesh oa = JMesh::Mesh();

  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->AssembleOA(oa, cvt_props_);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }

  return oa;
}


void Popups::ResetPatch()
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->ResetPatch();

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }

  TrimComponent();
}


void Popups::ResetMeshs()
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->ResetMeshs();

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}



static struct DataProps
{
  int id = -1;
  int parent_id = -1;
  E_FOLD_TYPE depend = E_FOLD_TYPE::DEFAULT;
  EVec2d v_rel_org = JUtil::ErrEVec2d();
  double v_u = -1.0;
  double v_v = -1.0;
  double v_w = -1.0;
  EVec2d h_rel_org = JUtil::ErrEVec2d();
  double h_u = -1.0;
  double h_v = -1.0;
  double h_w = -1.0;
};

void Popups::LoadData(const vector<std::string>& data)
{
  vector<Component*> root_children = root_->Children(E_FOLD_TYPE::GTYPE);
  for (Component* child : root_children)
    DeleteComponent({ child, root_, nullptr, E_FOLD_TYPE::GTYPE });

  vector<DataProps> data_props(data.size() - 1);

  for (int i = 0; i < data.size() - 1; ++i)
  {
    vector<std::string> elements = JUtil::Split(data[i], "/");
    for (std::string element : elements)
    {
      vector<std::string> e = JUtil::Split(element, " ");
      if (e[0] == "id")
        data_props[i].id = std::stoi(e[1]);
      else if (e[0] == "parent_id")
        data_props[i].parent_id = std::stoi(e[1]);
      else if (e[0] == "depend")
        data_props[i].depend = static_cast<E_FOLD_TYPE>(std::stoi(e[1]));
      else if (e[0] == "v")
      {
        if (e[1] == "rel_org")
          data_props[i].v_rel_org = JUtil::ParseEVec2d(e[2]);
        else if (e[1] == "u")
          data_props[i].v_u = stod(e[2]);
        else if (e[1] == "v")
          data_props[i].v_v = stod(e[2]);
        else if (e[1] == "w")
          data_props[i].v_w = stod(e[2]);
      }
      else if (e[0] == "h")
      {
        if (e[1] == "rel_org")
          data_props[i].h_rel_org = JUtil::ParseEVec2d(e[2]);
        else if (e[1] == "u")
          data_props[i].h_u = stod(e[2]);
        else if (e[1] == "v")
          data_props[i].h_v = stod(e[2]);
        else if (e[1] == "w")
          data_props[i].h_w = stod(e[2]);
      }
    }
  }

  vector<Component*> comps;
  comps.reserve(data.size());
  comps.push_back(root_);
  for (int i = 0; i < data_props.size(); ++i)
  {
    Component* parent = nullptr;
    Fold depend_fold;
    for (int j = 0; j < comps.size(); ++j)
    {
      if (data_props[i].parent_id == comps[j]->Id())
      {
        parent = comps[j];
        depend_fold = comps[j]->FoldEdge(data_props[i].depend);
        break;
      }
    }

    vector<Fold> folds =
    {
      {E_FOLD_TYPE::VTYPE, data_props[i].v_rel_org, JUtil::ErrEVec3d(), JUtil::ErrEVec3d(), data_props[i].v_u, data_props[i].v_v, data_props[i].v_w},
      {E_FOLD_TYPE::HTYPE, data_props[i].h_rel_org, JUtil::ErrEVec3d(), JUtil::ErrEVec3d(), data_props[i].h_u, data_props[i].h_v, data_props[i].h_w}
    };

    Component* comp = new Component(data_props[i].id, depend_fold, folds);
    parent->AddChild(data_props[i].depend, comp);
    comps.push_back(comp);
  }

  TrimComponent();
}


std::string Popups::SaveData()
{
  std::string data = "";

  queue<tuple<Component*, Component*, E_FOLD_TYPE>> que;
  que.push({ root_, nullptr, E_FOLD_TYPE::DEFAULT });

  while (!que.empty())
  {
    Component* comp = std::get<0>(que.front());
    Component* parent = std::get<1>(que.front());
    E_FOLD_TYPE depend = std::get<2>(que.front());
    que.pop();

    if (parent != nullptr)
    {
      data += "id " + std::to_string(comp->Id()) + "/";
      data += "parent_id " + std::to_string(parent->Id()) + "/";
      data += "depend " + std::to_string(static_cast<int>(depend)) + "/";
      data += "v rel_org " + JUtil::ToStringEVec2dNonBrackets(comp->FoldEdge(E_FOLD_TYPE::VTYPE).rel_org) + "/";
      data += "v u " + std::to_string(comp->FoldEdge(E_FOLD_TYPE::VTYPE).u) + "/";
      data += "v v " + std::to_string(comp->FoldEdge(E_FOLD_TYPE::VTYPE).v) + "/";
      data += "v w " + std::to_string(comp->FoldEdge(E_FOLD_TYPE::VTYPE).w) + "/";
      data += "h rel_org " + JUtil::ToStringEVec2dNonBrackets(comp->FoldEdge(E_FOLD_TYPE::HTYPE).rel_org) + "/";
      data += "h u " + std::to_string(comp->FoldEdge(E_FOLD_TYPE::HTYPE).u) + "/";
      data += "h v " + std::to_string(comp->FoldEdge(E_FOLD_TYPE::HTYPE).v) + "/";
      data += "h w " + std::to_string(comp->FoldEdge(E_FOLD_TYPE::HTYPE).w) + "\n";
    }

    for (int i = 0; i < comp->Children().size(); ++i)
    {
      E_FOLD_TYPE depend = static_cast<E_FOLD_TYPE>(i);
      for (int j = 0; j < comp->Children(i).size(); ++j)
        que.push({ comp->Child(i, j), comp, depend });
    }
  }

  return data;
}


static struct HitProperties
{
  Component* comp = nullptr;
  Component* parent = nullptr;
  Component* grand_parent = nullptr;
  E_FOLD_TYPE fold_type = E_FOLD_TYPE::DEFAULT;
  double dist = DBL_MAX;
};

static Component* sGetGrandParent
(
  const E_FOLD_TYPE& type,
  const tuple<Component*, Component*, Component*>& parents
)
{
  if (type == E_FOLD_TYPE::GTYPE)
    return std::get<static_cast<int>(E_FOLD_TYPE::GTYPE)>(parents);
  else if (type == E_FOLD_TYPE::VTYPE)
    return std::get<static_cast<int>(E_FOLD_TYPE::HTYPE)>(parents);
  else if (type == E_FOLD_TYPE::HTYPE)
    return std::get<static_cast<int>(E_FOLD_TYPE::VTYPE)>(parents);
}

static void sPushOAMechanismTuple
(
  int fold_idx,
  queue<tuple<Component*, Component*, Component*>>& que,
  const tuple<Component*, Component*, Component*>& grand_parents
)
{
  E_FOLD_TYPE fold_type = static_cast<E_FOLD_TYPE>(fold_idx);

  if (fold_type == E_FOLD_TYPE::VTYPE)
  {
    for (Component* child : std::get<2>(grand_parents)->Children(fold_type))
      que.push({ std::get<2>(grand_parents),
                 std::get<static_cast<int>(E_FOLD_TYPE::HTYPE)>(grand_parents), child });
  }
  else if (fold_type == E_FOLD_TYPE::HTYPE)
  {
    for (Component* child : std::get<2>(grand_parents)->Children(fold_type))
      que.push({ std::get<static_cast<int>(E_FOLD_TYPE::VTYPE)>(grand_parents),
               std::get<2>(grand_parents), child });
  }
  else if (fold_type == E_FOLD_TYPE::GTYPE)
  {
    for (Component* child : std::get<2>(grand_parents)->Children(fold_type))
      que.push({ std::get<static_cast<int>(E_FOLD_TYPE::VTYPE)>(grand_parents),
               std::get<static_cast<int>(E_FOLD_TYPE::HTYPE)>(grand_parents), child });
  }
}

static void sUpdateHitProperties
(
  Component* comp,
  const tuple<Component*, Component*, Component*>& comp_tuple,
  const E_FOLD_TYPE& fold_type,
  double dist,
  HitProperties& hit_props
)
{
  if (dist < hit_props.dist)
  {
    hit_props.comp = comp;
    hit_props.fold_type = fold_type;
    hit_props.dist = dist;
    hit_props.parent = std::get<2>(comp_tuple);
    hit_props.grand_parent = sGetGrandParent(fold_type, comp_tuple);
  }
}


bool Popups::HitFoldEdge(const JUtil::Ray& ray, PlaceProperties& place_props) const
{
  queue<tuple<Component*, Component*, Component*>> que;
  que.push({ root_, root_, root_ });

  HitProperties hit_props = { nullptr, nullptr, nullptr, E_FOLD_TYPE::DEFAULT, DBL_MAX };
  EVec3d min_dist_start_pos = JUtil::ErrEVec3d();
  EVec3d min_dist_end_pos = JUtil::ErrEVec3d();
  while (!que.empty())
  {
    tuple<Component*, Component*, Component*> comp_tuple = que.front();
    que.pop();

    double dist = DBL_MAX;
    E_FOLD_TYPE fold_type = E_FOLD_TYPE::DEFAULT;
    EVec3d start_pos = JUtil::ErrEVec3d();
    EVec3d end_pos = JUtil::ErrEVec3d();
    if (std::get<2>(comp_tuple)->HitFoldEdge(ray, place_props.radius, cvt_props_, dist, fold_type, start_pos, end_pos))
    {
      sUpdateHitProperties(std::get<2>(comp_tuple), comp_tuple, fold_type, dist, hit_props);
      min_dist_start_pos = start_pos;
      min_dist_end_pos = end_pos;
    }

    for (int i = 0; i < std::get<2>(comp_tuple)->Children().size(); i++)
      sPushOAMechanismTuple(i, que, comp_tuple);
  }

  if ((hit_props.comp != nullptr) && (hit_props.grand_parent != nullptr))
  {
    place_props.parent = hit_props.comp;
    place_props.grand_parent = hit_props.grand_parent;
    place_props.fold_type = hit_props.fold_type;
    place_props.start_pos = min_dist_start_pos;
    place_props.end_pos = min_dist_end_pos;
    return true;
  }

  place_props.parent = nullptr;
  place_props.grand_parent = nullptr;
  place_props.fold_type = E_FOLD_TYPE::DEFAULT;
  place_props.start_pos = JUtil::ErrEVec3d();
  place_props.end_pos = JUtil::ErrEVec3d();
  return false;
}


bool Popups::HitComponent(const JUtil::Ray& ray, DeleteProperties& delete_props) const
{
  queue<tuple<Component*, Component*, Component*>> que;
  que.push({ root_, root_, root_ });

  HitProperties hit_props = { nullptr, nullptr, nullptr, E_FOLD_TYPE::DEFAULT, DBL_MAX };
  while (!que.empty())
  {
    tuple<Component*, Component*, Component*> comp_tuple = que.front();
    que.pop();

    for (int i = 0; i < std::get<2>(comp_tuple)->Children().size(); i++)
    {
      for (int j = 0; j < std::get<2>(comp_tuple)->Children(i).size(); j++)
      {
        sPushOAMechanismTuple(i, que, comp_tuple);

        Component* child = std::get<2>(comp_tuple)->Child(i, j);
        double dist = DBL_MAX;
        if (child->Hit(ray, dist))
          sUpdateHitProperties(child, comp_tuple, static_cast<E_FOLD_TYPE>(i), dist, hit_props);
      }
    }
  }

  if ((hit_props.comp != nullptr) && (hit_props.parent != nullptr) && (hit_props.grand_parent != nullptr))
  {
    delete_props = { hit_props.comp, hit_props.parent, hit_props.grand_parent, hit_props.fold_type };
    return true;
  }

  delete_props = { nullptr, nullptr, nullptr, E_FOLD_TYPE::DEFAULT };
  return false;
}


bool Popups::HitComponent(const JUtil::Ray& ray, DeformProperties& deform_props) const
{
  DeleteProperties props = { nullptr, nullptr, nullptr, E_FOLD_TYPE::DEFAULT };

  if (HitComponent(ray, props))
  {
    deform_props.deform_comp = props.delete_comp;
    deform_props.parent = props.parent;
    deform_props.grand_parent = props.grand_parent;
    deform_props.fold_type = props.fold_type;

    deform_props.handles =
      deform_props.deform_comp->GetHandles(deform_props.handle_size, deform_props.slice);

    return true;
  }

  deform_props.deform_comp = nullptr;
  deform_props.parent = nullptr;
  deform_props.grand_parent = nullptr;
  deform_props.fold_type = E_FOLD_TYPE::DEFAULT;
  return false;
}


void Popups::IsConnect()
{
  queue<tuple<Component*, Rect3D, bool>> que;
  que.push({ root_, { vector<EVec3d>()}, false });

  while (!que.empty())
  {
    Component* comp = std::get<0>(que.front());
    Rect3D parent_rect = std::get<1>(que.front());
    bool root_h_child = std::get<2>(que.front());
    que.pop();

    if (comp->Id() != 0)
    {
      comp->IsConnect(parent_rect, cvt_props_);

      if (!root_h_child)
      {
        for (int i = 0; i < comp->Children().size(); i++)
        {
          for (Component* child : comp->Children(i))
            que.push({ child, comp->Rect(i), false });
        }
      }
    }
    else
    {
      for (Component* child : comp->Children(E_FOLD_TYPE::GTYPE))
      {
        que.push({ child, comp->Rect(E_FACE_TYPE::VTYPE), false });
        que.push({ child, comp->Rect(E_FACE_TYPE::HTYPE), true });
      }
    }
  }
}


void Popups::IsMaximumSize()
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    if (comp->Id() != 0)
      comp->IsMaximumSize(root_->FoldEdge(E_FOLD_TYPE::GTYPE));

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}


void Popups::IsMinimumLength()
{
  queue<Component*> que;
  que.push(root_);

  while (!que.empty())
  {
    Component* comp = que.front();
    que.pop();

    comp->IsMinimumLength(cvt_props_.nozzle_width);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push(child);
    }
  }
}


static queue<int> sSearchRoute(Component* root, int id)
{
  std::stack<int> route_stk;

  std::stack<pair<Component*, Component*>> stk;
  stk.push({ root, root });

  while (!stk.empty())
  {
    Component* comp = stk.top().first;
    Component* parent = stk.top().second;
    stk.pop();

    if (!route_stk.empty())
    {
      while (parent->Id() != route_stk.top())
        route_stk.pop();
    }

    if (comp->Id() == id)
      break;

    if (comp->IsExistChilren())
      route_stk.push(comp->Id());

    for (auto children : comp->Children())
    {
      for (auto child : children)
        stk.push({ child, comp });
    }
  }

  route_stk.push(id);
  std::stack<int> temp;
  while (!route_stk.empty())
  {
    temp.push(route_stk.top());
    route_stk.pop();
  }

  queue<int> route;
  while (!temp.empty())
  {
    route.push(temp.top());
    temp.pop();
  }

  return route;
}

static vector<Rect3D> sCalcChildRects(Component* root, Component* comp)
{
  queue<pair<Component*, bool>> que;
  que.push({ root, false });

  vector<Rect3D> child_rects = vector<Rect3D>();
  queue<int> route = sSearchRoute(root, comp->Id());
  while (!que.empty())
  {
    Component* child_comp = que.front().first;
    bool is_comp_child = que.front().second;
    que.pop();

    if (!is_comp_child)
    {
      if (child_comp->Id() == route.front())
        route.pop();
      else if (comp->Id() != child_comp->Id())
        child_rects.push_back(child_comp->ComponentRect());
    }

    for (vector<Component*> children : child_comp->Children())
    {
      if (comp->Id() == child_comp->Id())
      {
        for (Component* child : children)
          que.push({ child, true });
      }
      else
      {
        for (Component* child : children)
          que.push({ child, false });
      }
    }
  }

  return child_rects;
}

static vector<Rect3D> sChildrenRects(Component* comp)
{
  vector<Rect3D> child_rects(0);

  for (vector<Component*> children : comp->Children())
  {
    for (Component* child : children)
    {
      if (child->Err())
        break;

      child_rects.push_back(child->ComponentRect());
    }
  }

  return child_rects;
}


void Popups::TrimComponent()
{
  UpdateComponent(J_PI);
  IsConnect();
  IsMaximumSize();
  IsMinimumLength();

  queue<pair<Component*, int>> que;
  que.push({ root_, 0 });

  while (!que.empty())
  {
    Component* comp = que.front().first;
    int depth = que.front().second;
    que.pop();

    vector<Rect3D> rects = sCalcChildRects(root_, comp);
    vector<Rect3D> child_rects = sChildrenRects(comp);

    comp->TrimPatch(rects, child_rects, cvt_props_);

    for (vector<Component*> children : comp->Children())
    {
      for (Component* child : children)
        que.push({ child, depth + 1 });
    }
  }

  UpdateComponent(J_PI / 2.0);
}


Popups::Popups()
{
  // change base component size
  // root_ = new Component(width, height, depth);

  root_ = new Component(10.0000000, 7.3483635, 4.4757919); // bunny.stl

  mesh_ = JMesh::Mesh();
  angle_ = J_PI / 2.0;
  cvt_props_ = { 0.2, 0.05, 0.05, 0.02, 0.04 };
  TrimComponent();
  UpdateComponent(J_PI / 2.0);
}

