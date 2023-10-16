/* Component.cpp */
/* Component classes' function implementation */
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


#include "Component.h"


using std::vector;
using std::pair;


int Component::s_comp_id_generator_ = 0;
bool Component::s_trimmed_ = false;


inline static int sCast(const E_FOLD_TYPE& type) { return static_cast<int>(type); }
inline static int sCast(const E_FACE_TYPE& type) { return static_cast<int>(type); }


Component::Component(double width, double v_height, double h_height)
  : id_(GenerateMechId()), color_(GenerateColor(id_))
{
  err_ = false;

  faces_ = vector<Patch>
  {
    Patch(width, v_height, E_FACE_TYPE::VTYPE, true),
    Patch(width, h_height, E_FACE_TYPE::HTYPE, true)
  };

  folds_ = vector<Fold>(sCast(E_FOLD_TYPE::NUM_OF_TYPES));
  folds_[0] = { E_FOLD_TYPE::VTYPE, JUtil::ErrEVec2d(), JUtil::ErrEVec3d(), JUtil::ErrEVec3d(), -1.0, -1.0, -1.0 };
  folds_[1] = { E_FOLD_TYPE::HTYPE, JUtil::ErrEVec2d(), JUtil::ErrEVec3d(), JUtil::ErrEVec3d(), -1.0, -1.0, -1.0 };
  folds_[2] = { E_FOLD_TYPE::GTYPE, EVec2d::Zero(), EVec3d::Zero(), EVec3d::Zero(), width, v_height, h_height };

  children_ = vector<vector<Component*>>(3, vector<Component*>(0));

  rects_ = vector<Rect3D>
  {
    { vector<EVec3d>{ EVec3d::Zero(), EVec3d(width, 0.0, 0.0), EVec3d(width, v_height, 0.0), EVec3d(0.0, v_height, 0.0) }},
    { vector<EVec3d>{ EVec3d::Zero(), EVec3d(width, 0.0, 0.0), EVec3d(width, 0.0, h_height), EVec3d(0.0, 0.0, h_height) }}
  };

  trimmed_ = false;
  prev_folds_ = folds_;
  prev_rects_ = rects_;
}


Component::Component(const Fold& depend_fold, const EVec2d& scale)
  : id_(GenerateMechId()), color_(GenerateColor(id_))
{
  err_ = false;

  faces_ = vector<Patch>
  {
    Patch(depend_fold, E_FACE_TYPE::VTYPE, scale),
    Patch(depend_fold, E_FACE_TYPE::HTYPE, scale)
  };

  Initialize(depend_fold, scale);
  children_ = vector<vector<Component*>>(3, vector<Component*>(0));

  trimmed_ = false;
  prev_folds_ = folds_;
  prev_rects_ = rects_;
}


Component::Component(int id, const Fold& depend_fold, const vector<Fold>& folds)
  : id_(id), color_(GenerateColor(id_))
{
  MechIdGenerater(8);

  err_ = false;

  folds_ = vector<Fold>(2);

  for (int i = 0; i < 2; i++)
  {
    folds_[i].type = static_cast<E_FOLD_TYPE>(i);
    folds_[i].rel_org = folds[i].rel_org;
    folds_[i].u = folds[i].u;
    folds_[i].v = folds[i].v;
    folds_[i].w = folds[i].w;
  }

  rects_ = vector<Rect3D>
  {
    { vector<EVec3d>{ EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero() }},
    { vector<EVec3d>{ EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero() }}
  };

  Update(depend_fold, (J_PI / 2.0));

  for (int i = 0; i < 2; i++)
    folds_[i].org90 = folds_[i].org;

  faces_ = vector<Patch>
  {
    Patch(folds_[0].u, folds_[0].v, E_FACE_TYPE::VTYPE),
    Patch(folds_[1].u, folds_[1].w, E_FACE_TYPE::HTYPE)
  };

  children_ = vector<vector<Component*>>(3, vector<Component*>(0));

  trimmed_ = false;
  prev_folds_ = folds_;
  prev_rects_ = rects_;
}



Component::Component(const Component& src)
  : id_(src.id_), color_(src.color_)
{
  this->err_ = src.err_;
  this->faces_ = src.faces_;
  this->folds_ = src.folds_;
  this->rects_ = src.rects_;
  this->children_ = src.children_;
  this->trimmed_ = src.trimmed_;
  this->prev_folds_ = src.prev_folds_;
  this->prev_rects_ = src.prev_rects_;
}


double Component::Left() const
{
  if (id_ == 0)
    return folds_[2].org[0];
  else
    return folds_[0].org[0];
}


double Component::Right() const
{
  if (id_ == 0)
    return folds_[2].org[0] + folds_[2].u;
  else
    return folds_[0].org[0] + folds_[0].u;
}


bool Component::IsExistChilren() const
{
  for (vector<Component*> children : children_)
  {
    if (children.size() != 0)
      return true;
  }

  return false;
}


static void sDrawMountLine(double left, double right, const vector<vector<Component*>>& children)
{
  glDisable(GL_LIGHTING);
  glLineWidth(1.0f);
  glColor3f(0.0f, 0.0f, 0.0f);
  glBegin(GL_LINES);
  int fold_idx = sCast(E_FOLD_TYPE::GTYPE);
  for (int i = 0; i < children[fold_idx].size() + 1; i++)
  {
    double start_x = (i == 0) ? (left) : (children[fold_idx][i - 1]->Right());
    double end_x = (i == children[fold_idx].size()) ? (right) : (children[fold_idx][i]->Left());

    glVertex3f((float)start_x, 0.0f, 0.0f);
    glVertex3f((float)end_x, 0.0f, 0.0f);
  }
  glEnd();
  glFlush();
  glEnable(GL_LIGHTING);
}

static void sDrawConnectionFold(const vector<Rect3D>& rects)
{
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  JDraw::SetEVec3dToColor3fv(EVec3d(0.0, 0.0, 0.0));
  EVec3d p0 = rects[0].p[2] + EVec3d(0.0, 0.0, 0.0);
  EVec3d p1 = rects[0].p[3] + EVec3d(0.0, 0.0, 0.0);
  JDraw::SetEVec3dToVertex3fv(p0);
  JDraw::SetEVec3dToVertex3fv(p1);
  glEnd();
  glEnable(GL_LIGHTING);
}


void Component::Draw2D(double angle, const EVec3f& highlight) const
{
  EVec3f color = err_ ? EVec3f(1.0f, 0.0f, 0.0f) : color_;

  if (id_ == 0)
  {
    int fold_idx = sCast(E_FOLD_TYPE::GTYPE);
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].Draw2D(folds_[fold_idx].org90, color + highlight, angle);

    if (Trimmed())
      sDrawMountLine(Left(), Right(), children_);
  }
  else
  {
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].Draw2D(folds_[i].org90, color + highlight, angle);

    sDrawConnectionFold(rects_);
  }
}


void Component::Draw3D(double angle, const ConvertProperties& cvt_props) const
{
  vector<Fold>   draw_folds = (trimmed_) ? (folds_) : (prev_folds_);
  vector<Rect3D> draw_rects = (trimmed_) ? (rects_) : (prev_rects_);

  if (id_ == 0)
  {
    int fold_idx = sCast(E_FOLD_TYPE::GTYPE);
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].Draw3D(draw_folds[fold_idx].org, color_, angle, cvt_props);

    if (Trimmed())
      sDrawMountLine(Left(), Right(), children_);
  }
  else
  {
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].Draw3D(draw_folds[i].org, color_, angle, cvt_props);

    if (!faces_[0].LoadedMesh())
      sDrawConnectionFold(draw_rects);
  }
}


static Fold sUpdateFold(const Fold& depend_fold, const Fold& update_fold, const Fold& another_fold, double angle)
{
  Fold fold = update_fold;
  E_FACE_TYPE face_type = (fold.type == E_FOLD_TYPE::VTYPE) ? E_FACE_TYPE::HTYPE : E_FACE_TYPE::VTYPE;

  fold.org = depend_fold.org
    + (fold.rel_org[0] * Patch::Axis(face_type, angle).u_dir)
    + (fold.rel_org[1] * Patch::Axis(face_type, angle).v_dir);

  fold.org90 = depend_fold.org90
    + (fold.rel_org[0] * Patch::Axis(face_type, J_PI / 2.0).u_dir)
    + (fold.rel_org[1] * Patch::Axis(face_type, J_PI / 2.0).v_dir);

  if (update_fold.type == E_FOLD_TYPE::VTYPE)
    fold.w = depend_fold.w - another_fold.w;
  else
    fold.v = depend_fold.v - another_fold.v;

  return fold;
}


void Component::UpdateFolds(const Fold& depend_fold, double angle)
{
  Fold v_fold = folds_[sCast(E_FOLD_TYPE::VTYPE)];
  Fold h_fold = folds_[sCast(E_FOLD_TYPE::HTYPE)];
  folds_[sCast(E_FOLD_TYPE::VTYPE)] = sUpdateFold(depend_fold, v_fold, h_fold, angle);
  folds_[sCast(E_FOLD_TYPE::HTYPE)] = sUpdateFold(depend_fold, h_fold, v_fold, angle);
}


static Rect3D sUpdateMountRect(const Fold& fold, const E_FACE_TYPE& type, double angle)
{
  Rect3D rect = { vector<EVec3d>(4) };
  EVec3d org = fold.org;
  EVec3d offset = ((type == E_FACE_TYPE::VTYPE) ? (fold.v) : (fold.w)) * Patch::Axis(type, angle).v_dir;

  rect.p[0] = org;
  rect.p[1] = org + EVec3d(fold.u, 0.0, 0.0);
  rect.p[2] = org + EVec3d(fold.u, 0.0, 0.0) + offset;
  rect.p[3] = org + offset;

  return rect;
}

static Rect3D sUpdateRect(const Fold& fold, const E_FACE_TYPE& type, double angle)
{
  Rect3D rect = { vector<EVec3d>(4) };
  EVec3d org = fold.org;
  EVec3d offset = ((fold.type == E_FOLD_TYPE::VTYPE) ? (fold.v) : (fold.w)) * Patch::Axis(type, angle).v_dir;

  rect.p[0] = org;
  rect.p[1] = org + EVec3d(fold.u, 0.0, 0.0);
  rect.p[2] = org + EVec3d(fold.u, 0.0, 0.0) + offset;
  rect.p[3] = org + offset;

  return rect;
}


void Component::UpdateRects(double angle)
{
  rects_[sCast(E_FACE_TYPE::VTYPE)] = sUpdateRect(folds_[sCast(E_FOLD_TYPE::VTYPE)], E_FACE_TYPE::VTYPE, angle);
  rects_[sCast(E_FACE_TYPE::HTYPE)] = sUpdateRect(folds_[sCast(E_FOLD_TYPE::HTYPE)], E_FACE_TYPE::HTYPE, angle);
}


void Component::Update(const Fold& depend_fold, double angle)
{
  if (id_ == 0)
  {
    folds_[sCast(E_FOLD_TYPE::GTYPE)].org = EVec3d::Zero();
    folds_[sCast(E_FOLD_TYPE::GTYPE)].org90 = EVec3d::Zero();
    rects_[sCast(E_FACE_TYPE::VTYPE)] = sUpdateMountRect(folds_[sCast(E_FOLD_TYPE::GTYPE)], E_FACE_TYPE::VTYPE, angle);
    rects_[sCast(E_FACE_TYPE::HTYPE)] = sUpdateMountRect(folds_[sCast(E_FOLD_TYPE::GTYPE)], E_FACE_TYPE::HTYPE, angle);

    return;
  }

  UpdateFolds(depend_fold, angle);
  UpdateRects(angle);
}


void Component::IsConnect(const Rect3D& parent_rect, const ConvertProperties& cvt_props)
{
  double parent_max_x = parent_rect.p[1][0];
  double parent_min_x = parent_rect.p[0][0];
  double parent_max_z = (parent_rect.p[0][2] > parent_rect.p[3][2]) ? parent_rect.p[0][2] : parent_rect.p[3][2];
  double parent_min_z = (parent_rect.p[0][2] < parent_rect.p[3][2]) ? parent_rect.p[0][2] : parent_rect.p[3][2];

  if ((ComponentRect().p[0][0] <= parent_min_x) && (parent_max_x <= ComponentRect().p[1][0]))
    err_ = true;
  else if ((ComponentRect().p[1][0] <= parent_min_x) || (parent_max_x <= ComponentRect().p[0][0]))
    err_ = true;
  else if ((ComponentRect().p[3][2] <= parent_min_z) && (parent_max_z < ComponentRect().p[0][2]))
    err_ = true;
  else
    err_ = false;

  for (const Fold& fold : folds_)
  {
    if (fold.v < cvt_props.thickness || fold.w < cvt_props.thickness)
    {
      err_ = true;
      break;
    }
  }

  if (err_)
    std::cout << "Not connect fold line\n" << std::endl;
}


void Component::IsMaximumSize(const Fold& root_fold)
{
  if (err_)
    return;

  EVec2d max = EVec2d(root_fold.u, root_fold.w);
  EVec2d min = EVec2d(0.0, -root_fold.v);

  if (ComponentRect().p[0][0] < min[0] || max[0] < ComponentRect().p[1][0] ||
      ComponentRect().p[3][2] < min[1] || max[1] < ComponentRect().p[0][2])
    err_ = true;
  else
    err_ = false;

  if (err_)
    std::cout << "Over base card\n" << std::endl;
}


void Component::IsMinimumLength(double nozzle_width)
{
  if (err_)
    return;

  for (int i = 0; i < folds_.size(); i++)
  {
    int fold_idx = sCast(folds_[i].type);
    if (children_[fold_idx].size() > 0)
    {
      for (int j = 0; j < children_[fold_idx].size() + 1; j++)
      {
        double left = (j == 0) ? (this->Left()) : (children_[fold_idx][j - 1]->Right());
        double right = (j == children_[fold_idx].size()) ? (this->Right()) : (children_[fold_idx][j]->Left());

        bool err = false;
        if ((left < right))
        {
          err = (right - left) < nozzle_width;
          err_ = err;
        }
        else if ((j == 0) && (left > right))
        {
          left = this->Left();
          right = children_[fold_idx][j]->Right();
          err = (right - left) < nozzle_width;

          if (!(children_[fold_idx][j]->Err()))
            children_[fold_idx][j]->Err(err);
        }
        else if ((j == children_[fold_idx].size()) && (left > right))
        {
          left = children_[fold_idx][j - 1]->Left();
          right = this->Right();
          err = (right - left) < nozzle_width;

          if (!(children_[fold_idx][j - 1]->Err()))
            children_[fold_idx][j - 1]->Err(err);
        }
      }
    }
    else
    {
      double left = this->Left();
      double right = this->Right();

      if ((left < right))
        err_ = (right - left) < nozzle_width;
    }

    if (err_)
    {
      std::cout << id_ << ": Some areas are too thin\n" << std::endl;
      return;
    }
  }
}


static bool sIsIntersectX(double rect_0x, double rect_2x, double other_0x, double other_2x)
{
  return ((rect_0x <= other_0x) && (other_0x <= rect_2x) && (rect_2x <= other_2x)) ||
    ((other_0x <= rect_0x) && (rect_0x <= other_2x) && (other_2x <= rect_2x)) ||
    ((rect_0x >= other_0x) && (other_2x >= rect_2x));
}

static bool sIsIntersectZ(double rect_0z, double rect_2z, double other_0z, double other_2z, const E_FACE_TYPE& type)
{
  if (type == E_FACE_TYPE::VTYPE)
  {
    return ((rect_2z <= other_2z) && (other_2z <= rect_0z) && (rect_0z <= other_0z)) ||
      ((other_2z <= rect_2z) && (rect_2z <= other_0z) && (other_0z <= rect_0z)) ||
      ((rect_2z >= other_2z) && (other_0z >= rect_0z));
  }
  else
  {
    return ((rect_0z <= other_2z) && (other_2z <= rect_2z) && (rect_2z <= other_0z)) ||
      ((other_2z <= rect_0z) && (rect_0z <= other_0z) && (other_0z <= rect_2z)) ||
      ((rect_0z >= other_2z) && (other_0z >= rect_2z));
  }
}

static bool sIsOverlapping(const Rect3D& rect, vector<Rect3D>& other_rects, const E_FACE_TYPE& type)
{
  double rect_0x = rect.p[0][0];
  double rect_0z = rect.p[0][2];
  double rect_2x = rect.p[2][0];
  double rect_2z = rect.p[2][2];

  std::queue<int> deletes;
  for (int i = 0; i < other_rects.size(); i++)
  {
    double other_0x = other_rects[i].p[0][0];    double other_0z = other_rects[i].p[0][2];
    double other_2x = other_rects[i].p[2][0];    double other_2z = other_rects[i].p[2][2];

    if (sIsIntersectX(rect_0x, rect_2x, other_0x, other_2x) &&
        sIsIntersectZ(rect_0z, rect_2z, other_0z, other_2z, type))
    {
      deletes.push(i);
    }
  }

  if (!deletes.empty())
  {
    int delete_count = 0;
    while (!deletes.empty())
    {
      int idx = deletes.front();
      deletes.pop();

      other_rects.erase(other_rects.begin() + idx - delete_count);
      delete_count++;
    }

    return true;
  }
  else
  {
    return false;
  }
}


void Component::TrimPatch
(
  const vector<Rect3D>& other_rects,
  const vector<Rect3D>& child_rects,
  const ConvertProperties& cvt_props
)
{
  vector<Rect3D> others = other_rects;

  if (!err_)
  {
    bool v_err = sIsOverlapping(rects_[sCast(E_FACE_TYPE::VTYPE)], others, E_FACE_TYPE::VTYPE);
    bool h_err = sIsOverlapping(rects_[sCast(E_FACE_TYPE::HTYPE)], others, E_FACE_TYPE::HTYPE);

    err_ = v_err || h_err;

    if (err_)
    {
      std::cout << "Crossing component each other\n" << std::endl;
    }
  }

  vector<Rect3D> rects(others.size() + child_rects.size());
  for (int i = 0; i < rects.size(); i++)
  {
    if (i < others.size())
      rects[i] = others[i];
    else
      rects[i] = child_rects[i - others.size()];
  }

  if (id_ == 0)
  {
    faces_[sCast(E_FACE_TYPE::VTYPE)].TrimChildren(rects_[sCast(E_FACE_TYPE::VTYPE)], rects, cvt_props, true);
    faces_[sCast(E_FACE_TYPE::HTYPE)].TrimChildren(rects_[sCast(E_FACE_TYPE::HTYPE)], rects, cvt_props, true);
  }
  else
  {
    faces_[sCast(E_FACE_TYPE::VTYPE)].TrimChildren(rects_[sCast(E_FACE_TYPE::VTYPE)], rects, cvt_props, false);
    faces_[sCast(E_FACE_TYPE::HTYPE)].TrimChildren(rects_[sCast(E_FACE_TYPE::HTYPE)], rects, cvt_props, false);
  }

  trimmed_ = true;
  Trimmed(true);
}


static void sSortChildren(vector<Component*>& children)
{
  for (int i = 1; i < children.size(); i++)
  {
    for (int j = 0; j < children.size() - i; j++)
    {
      if (children[j]->Left() > children[j + 1]->Left())
      {
        Component* tmp = children[j];
        children[j] = children[j + 1];
        children[j + 1] = tmp;
      }
    }
  }
}


void Component::AddChild(const E_FOLD_TYPE& fold_type, Component* comp)
{
  if (sCast(fold_type) <= sCast(E_FOLD_TYPE::DEFAULT) ||
      sCast(fold_type) >= sCast(E_FOLD_TYPE::NUM_OF_TYPES))
  {
    std::cout << "E_FOLD_TYPE Error: OAMechanism::AddChild(E_FOLD_EDGE_TYPE "
      << sCast(fold_type) << ", OAMechanism *" << comp->id_ << ")\n";
  }

  children_[sCast(fold_type)].push_back(comp);
  sSortChildren(children_[sCast(fold_type)]);
}


void Component::DeleteChild(const E_FOLD_TYPE& fold_type, int mech_id)
{
  int fold_edge_idx = sCast(fold_type);
  int delete_idx = -1;
  for (int i = 0; i < children_[fold_edge_idx].size(); i++)
  {
    if (children_[fold_edge_idx][i]->Id() == mech_id)
      delete_idx = i;
  }

  children_[fold_edge_idx].erase(children_[fold_edge_idx].begin() + delete_idx);
}

static void sAdjustChildrenFolds(vector<vector<Component*>>& children, double move_len)
{
  for (vector<Component*> childs : children)
  {
    for (Component* child : childs)
    {
      child->FoldEdgePtr(E_FOLD_TYPE::VTYPE)->rel_org[0] += move_len;
      child->FoldEdgePtr(E_FOLD_TYPE::HTYPE)->rel_org[0] += move_len;
    }
  }
}

static void sMoveLeft(vector<Fold>& folds, double move_len, double nozzle_width)
{
  int v_type = sCast(E_FOLD_TYPE::VTYPE);
  int h_type = sCast(E_FOLD_TYPE::HTYPE);

  if (folds[v_type].u + move_len > nozzle_width)
  {
    folds[v_type].org[0] -= move_len;
    folds[v_type].rel_org[0] -= move_len;
    folds[v_type].u += move_len;
    folds[h_type].org[0] -= move_len;
    folds[h_type].rel_org[0] -= move_len;
    folds[h_type].u += move_len;
  }
}

static void sMoveUp(vector<Fold>& folds, double move_len, double thickness)
{
  int v_type = sCast(E_FOLD_TYPE::VTYPE);
  int h_type = sCast(E_FOLD_TYPE::HTYPE);

  if (folds[v_type].v + move_len > thickness)
  {
    folds[h_type].org[1] += move_len;
    folds[h_type].rel_org[1] += move_len;
    folds[h_type].v -= move_len;
    folds[v_type].v += move_len;
  }
}

static void sMoveDepth(vector<Fold>& folds, double move_len, double thickness)
{
  int v_type = sCast(E_FOLD_TYPE::VTYPE);
  int h_type = sCast(E_FOLD_TYPE::HTYPE);

  if (folds[h_type].w + move_len > thickness)
  {
    folds[v_type].org[2] += move_len;
    folds[v_type].rel_org[1] += move_len;
    folds[v_type].w -= move_len;
    folds[h_type].w += move_len;
  }
}

static void sMoveRight(vector<Fold>& folds, double move_len, double nozzle_width)
{
  int v_type = sCast(E_FOLD_TYPE::VTYPE);
  int h_type = sCast(E_FOLD_TYPE::HTYPE);

  if (folds[v_type].u + move_len > nozzle_width)
  {
    folds[v_type].u += move_len;
    folds[h_type].u += move_len;
  }
}


double Component::DeformPatch(const E_DEFORM_HANDLE& handle, double move_len, const ConvertProperties& cvt_props)
{
  if (trimmed_)
  {
    trimmed_ = false;
    prev_folds_ = folds_;
    prev_rects_ = rects_;
    Trimmed(false);
  }

  switch (handle)
  {
  case E_DEFORM_HANDLE::LEFT:
    sMoveLeft(folds_, move_len, cvt_props.nozzle_width);
    sAdjustChildrenFolds(children_, move_len);
    break;

  case E_DEFORM_HANDLE::UP:
    sMoveUp(folds_, move_len, cvt_props.thickness);
    break;

  case E_DEFORM_HANDLE::DEPTH:
    sMoveDepth(folds_, move_len, cvt_props.thickness);
    break;

  case E_DEFORM_HANDLE::RIGHT:
    sMoveRight(folds_, move_len, cvt_props.nozzle_width);
    break;

  default:
    break;
  }

  for (int i = 0; i < 2; i++)
    folds_[i].org90 = folds_[i].org;

  for (Patch& face : faces_)
  {
    double width = Left() - Right();
    double height = (face.Type() == E_FACE_TYPE::VTYPE) ? (folds_[0].v) : (folds_[0].w);
    face.Deform(handle, move_len, width, height);
  }

  return move_len;
}


Rect3D Component::ComponentRect() const
{
  Rect3D rect = { vector<EVec3d>(4) };
  rect.p[0] = rects_[0].p[0];
  rect.p[1] = rects_[0].p[1];
  rect.p[2] = rects_[1].p[1];
  rect.p[3] = rects_[1].p[0];

  return rect;
}


void Component::FillPlane(const ConvertProperties& cvt_props)
{
  faces_[sCast(E_FACE_TYPE::VTYPE)].FillPlane(rects_[sCast(E_FACE_TYPE::VTYPE)], cvt_props);
  faces_[sCast(E_FACE_TYPE::HTYPE)].FillPlane(rects_[sCast(E_FACE_TYPE::HTYPE)], cvt_props);
}



void Component::CutPrintableOAFace(const JMesh::Mesh& parent_mesh, const Fold& depend_fold, const ConvertProperties& cvt_props)
{
  std::cout << "Trimming Patch Outlines (" << id_ << ")" << std::endl;

  if (id_ == 0)
  {
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].CutPrintableOAFace(parent_mesh, rects_[i], depend_fold, cvt_props);
  }
  else
  {
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].CutPrintableOAFace(parent_mesh, rects_[i], depend_fold, cvt_props);
  }

  std::cout << std::endl;
}


void Component::GenerateInflatedPatch(const ConvertProperties& cvt_props)
{
  if (id_ != 0)
  {
    for (int i = 0; i < faces_.size(); i++)
      faces_[i].GenerateInflatedPatch(rects_[i], folds_[i].org, cvt_props);
  }
}


JMesh::Mesh Component::InflatedPatch() const
{
  if (id_ != 0)
  {
    JMesh::Mesh v_patch = faces_[static_cast<int>(E_FACE_TYPE::VTYPE)].InflatedPatch();
    JMesh::Mesh h_patch = faces_[static_cast<int>(E_FACE_TYPE::HTYPE)].InflatedPatch();
    return JCSG::Union(v_patch, h_patch);
  }
}


void Component::TrimOutlines
(
  const E_FACE_TYPE& type,
  const JMesh::Mesh& outlines,
  const ConvertProperties& cvt_props
)
{
  if (id_ == 0)
  {
    std::cout << "" << id_ << ": ";
    faces_[sCast(type)].TrimOutlines(outlines, folds_[sCast(E_FOLD_TYPE::GTYPE)].org, rects_[sCast(type)].p, cvt_props);
  }
  else
  {
    if (outlines.Exist())
      std::cout << "" << id_ << ": ";

    faces_[sCast(type)].TrimOutlines(outlines, folds_[sCast(type)].org, rects_[sCast(type)].p, cvt_props);
  }
}


static JMesh::Mesh sCalcChildSweepRegion
(
  const EVec3d& max_pos,
  const vector<Rect3D>& rects,
  const E_FACE_TYPE& type,
  const ConvertProperties& cvt_props
)
{
  EVec3d x = EVec3d::UnitX();
  EVec3d y = -EVec3d::UnitY();
  EVec3d z = EVec3d::UnitZ();

  EVec3d org1 = (type == E_FACE_TYPE::HTYPE) ?
    EVec3d(rects[1].p[0][0] - cvt_props.gap, max_pos[1] + J_CSG_OFFSET, -J_CSG_OFFSET) :
    EVec3d(rects[1].p[0][0] - cvt_props.gap, rects[1].p[0][1] + J_CSG_OFFSET, rects[1].p[0][2]);

  double w1 = (rects[0].p[0] - rects[0].p[1]).norm() + (2.0 * cvt_props.gap);
  double h1 = ((type == E_FACE_TYPE::HTYPE) ? fabs(max_pos[1] - rects[0].p[0][1]) : fabs(rects[1].p[0][1])) + J_CSG_OFFSET;
  double d1 = ((type == E_FACE_TYPE::HTYPE) ? fabs(max_pos[2]) : fabs(max_pos[2] - rects[1].p[0][2])) + J_CSG_OFFSET;

  JMesh::Mesh region1 = JMesh::RectangularMesh({ org1, x, y, z, w1, h1, d1 });

  EVec3d org2 = (type == E_FACE_TYPE::HTYPE) ?
    EVec3d(rects[1].p[0][0] - cvt_props.gap, max_pos[1] + J_CSG_OFFSET, -J_CSG_OFFSET) :
    EVec3d(rects[1].p[0][0] - cvt_props.gap, rects[1].p[0][1], -J_CSG_OFFSET);

  double w2 = (rects[0].p[0] - rects[0].p[1]).norm() + (2.0 * cvt_props.gap);
  double h2 = ((type == E_FACE_TYPE::HTYPE) ? fabs(max_pos[1]) : fabs(rects[0].p[3][1])) + J_CSG_OFFSET;
  double d2 = ((type == E_FACE_TYPE::HTYPE) ? fabs(rects[1].p[3][2]) : fabs(max_pos[2])) + J_CSG_OFFSET;

  JMesh::Mesh region2 = JMesh::RectangularMesh({ org2, x, y, z, w2, h2, d2 });
  return JCSG::Union(region1, region2);
}


vector<pair<double, double>> sGetConcaveFoldList(const vector<Component*>& children, const Fold& fold)
{
  vector<pair<double, double>> concave_fold_list(0);
  if (children.size() > 0)
  {
    double rel = fold.org[0];
    for (int i = 0; i < children.size() + 1; i++)
    {
      double left = (i == 0) ? (fold.org[0]) : (children[i - 1]->Right());
      double right = (i == children.size()) ? (fold.org[0] + fold.u) : (children[i]->Left());

      if (((i == 0) && (left > right)) || ((i == children.size()) && (left > right)))
        continue;

      pair<double, double> concave_fold = { left - rel, right - rel };
      concave_fold_list.push_back(concave_fold);
    }
  }
  else
  {
    concave_fold_list.push_back({ 0.0, fold.u });
  }

  return concave_fold_list;
}


bool Component::Segment
(
  JMesh::Mesh& higher_space,
  const JMesh::Mesh& mesh,
  const E_FACE_TYPE& type,
  const ConvertProperties& cvt_props,
  const EVec3d& max_pos
)
{
  if (id_ == 0)
    return false;

  std::cout << id_ << ": ";
  Fold fold = (type == E_FACE_TYPE::VTYPE) ? (folds_[sCast(E_FOLD_TYPE::VTYPE)]) : (folds_[sCast(E_FOLD_TYPE::HTYPE)]);
  vector<Component*> children = (type == E_FACE_TYPE::VTYPE) ? (children_[sCast(E_FOLD_TYPE::VTYPE)]) : (children_[sCast(E_FOLD_TYPE::HTYPE)]);
  vector<pair<double, double>> concave_fold_list = sGetConcaveFoldList(children, fold);

  int face = sCast(type);
  bool one_more = faces_[face].SegmentSweepRegion(higher_space, mesh, rects_[face], cvt_props, max_pos, concave_fold_list);

  higher_space = JCSG::Union(higher_space, sCalcChildSweepRegion(max_pos, rects_, type, cvt_props));
  return one_more;
}


void Component::AssembleOA(JMesh::Mesh& oa, const ConvertProperties& cvt_props) const
{
  if (id_ == 0)
  {
    JMesh::Mesh vface = faces_[sCast(E_FACE_TYPE::VTYPE)].Assemble(folds_[sCast(E_FOLD_TYPE::GTYPE)].org);
    JMesh::Mesh hface = faces_[sCast(E_FACE_TYPE::HTYPE)].Assemble(folds_[sCast(E_FOLD_TYPE::GTYPE)].org);

    igl::copyleft::cgal::CSGTree tree = { {vface.V(), vface.F()}, {hface.V(), hface.F()}, "u" };
    oa = JMesh::Mesh(tree.cast_V<EMatXd>(), tree.F());
  }
  else
  {
    JMesh::Mesh vface = faces_[sCast(E_FACE_TYPE::VTYPE)].Assemble(folds_[sCast(E_FOLD_TYPE::VTYPE)].org);
    JMesh::Mesh hface = faces_[sCast(E_FACE_TYPE::HTYPE)].Assemble(folds_[sCast(E_FOLD_TYPE::HTYPE)].org);

    igl::copyleft::cgal::CSGTree tree = { { {vface.V(), vface.F()}, {hface.V(), hface.F()}, "u" }, { oa.V(), oa.F() }, "u" };
    oa = JMesh::Mesh(tree.cast_V<EMatXd>(), tree.F());
  }

  std::cout << "(" << id_ << ") - ";
}


void Component::ResetPatch()
{
  for (int i = 0; i < faces_.size(); i++)
    faces_[i].ResetPatch();
}

void Component::ResetMeshs()
{
  for (int i = 0; i < faces_.size(); i++)
    faces_[i].ResetMeshs();
}


static struct HitFoldEdgeProperties
{
  double min_dist;
  E_FOLD_TYPE min_dist_fold_type;
  EVec3d min_dist_start_pos;
  EVec3d min_dist_end_pos;
};

static void sHitCylinder
(
  const JUtil::Ray& ray,
  const EVec3d& org,
  const E_FOLD_TYPE& fold_type,
  double top_x,
  double base_x,
  double radius,
  HitFoldEdgeProperties& hit_props
)
{
  EVec3d top = EVec3d(top_x, org[1], org[2]);
  EVec3d base = EVec3d(base_x, org[1], org[2]);

  double d = DBL_MAX;
  JMesh::Cylinder cylinder = JMesh::GetCylinder(top, base, radius, 0);
  if (JMesh::HitCylinder(ray, cylinder, d))
  {
    if (d < hit_props.min_dist)
    {
      hit_props.min_dist = d;
      hit_props.min_dist_fold_type = fold_type;
      hit_props.min_dist_start_pos = top;
      hit_props.min_dist_end_pos = base;
    }
  }
}


bool Component::HitFoldEdge
(
  const JUtil::Ray& ray,
  double radius,
  const ConvertProperties& cvt_props,
  double& dist,
  E_FOLD_TYPE& fold_type,
  EVec3d& start,
  EVec3d& end
) const
{
  HitFoldEdgeProperties hit_props = { dist, E_FOLD_TYPE::DEFAULT, JUtil::ErrEVec3d(), JUtil::ErrEVec3d() };
  for (int i = 0; i < folds_.size(); i++)
  {
    int fold_idx = sCast(folds_[i].type);
    if ((folds_[fold_idx].u * 0.8 > (3.0 * cvt_props.nozzle_width)) &&
        (folds_[fold_idx].v * 0.35 > cvt_props.thickness) && (folds_[fold_idx].v * 0.65 > cvt_props.thickness) &&
        (folds_[fold_idx].w * 0.35 > cvt_props.thickness) && (folds_[fold_idx].w * 0.65 > cvt_props.thickness))
    {
      if (children_[fold_idx].size() == 0)
      {
        sHitCylinder(ray, folds_[i].org, folds_[i].type, this->Left(), this->Right(), radius, hit_props);
        continue;
      }

      for (int j = 0; j < children_[fold_idx].size() + 1; j++)
      {
        double top_x = (j == 0) ? (this->Left()) : (children_[fold_idx][j - 1]->Right());
        double base_x = (j == children_[fold_idx].size()) ? (this->Right()) : (children_[fold_idx][j]->Left());

        sHitCylinder(ray, folds_[i].org, folds_[i].type, top_x, base_x, radius, hit_props);
      }
    }
  }

  if (hit_props.min_dist_fold_type != E_FOLD_TYPE::DEFAULT)
  {
    dist = hit_props.min_dist;
    fold_type = hit_props.min_dist_fold_type;
    start = hit_props.min_dist_start_pos;
    end = hit_props.min_dist_end_pos;
    return true;
  }

  return false;
}


bool Component::Hit(const JUtil::Ray& ray, double& dist) const
{
  double min_dist = dist;
  for (const Patch& face : faces_)
  {
    EVec3d org = EVec3d::Zero();
    if (face.Type() == E_FACE_TYPE::VTYPE)
      org = folds_[sCast(E_FOLD_TYPE::VTYPE)].org;
    else if (face.Type() == E_FACE_TYPE::HTYPE)
      org = folds_[sCast(E_FOLD_TYPE::HTYPE)].org;

    double d = DBL_MAX;
    if (face.Hit(ray, org, d))
    {
      if (d < min_dist)
        min_dist = d;
    }
  }

  if (min_dist < dist)
  {
    dist = min_dist;
    return true;
  }

  return false;
}


vector<JMesh::Arrow> Component::GetHandles(const JMesh::ArrowSize& arrow_size, int slice) const
{
  vector<JMesh::Arrow> handles(4);
  for (int i = 0; i < handles.size(); i++)
  {
    EVec3d pos = folds_[0].org + EVec3d(0.0, folds_[0].v, 0.0);
    if (0 < i && i < handles.size() - 1)
    {
      if (i == 1)
        pos += EVec3d(folds_[0].u / 2.0, arrow_size.bar_radius, 0.0);
      else
        pos += EVec3d(folds_[0].u / 2.0, 0.0, arrow_size.bar_radius);
    }
    else if (i >= handles.size() - 1)
    {
      pos += EVec3d(folds_[0].u, 0.0, 0.0);
    }

    EVec3d dir = EVec3d::Zero();
    dir[i % 3] = (i == 0) ? -1.0 : 1.0;

    handles[i] = { pos, dir, arrow_size, slice };
  }

  return handles;
}


EVec3f Component::GenerateColor(int id)
{
  if (id == 0)
    return EVec3f(0.8f, 0.8f, 0.8f);

  std::srand(id);
  return EVec3f(static_cast<float>(std::rand() % 100) / 100.0f,
                static_cast<float>(std::rand() % 100) / 100.0f,
                static_cast<float>(std::rand() % 100) / 100.0f);
}


void Component::Initialize(const Fold& depend_fold, const EVec2d& scale)
{
  folds_ = vector<Fold>(2);

  for (int i = 0; i < 2; i++)
  {
    folds_[i].type = static_cast<E_FOLD_TYPE>(i);
    folds_[i].u = scale[0] * depend_fold.u;
    folds_[i].rel_org[0] = depend_fold.rel_org[0] + (((1.0 - scale[0]) / 2.0) * depend_fold.u);

    if (folds_[i].type == E_FOLD_TYPE::VTYPE)
    {
      folds_[i].v = scale[1] * depend_fold.v;
      folds_[i].w = depend_fold.w - (scale[1] * depend_fold.w);
      folds_[i].rel_org[1] = scale[1] * depend_fold.w;
    }
    else if (folds_[i].type == E_FOLD_TYPE::HTYPE)
    {
      folds_[i].v = depend_fold.v - (scale[1] * depend_fold.v);
      folds_[i].w = scale[1] * depend_fold.w;
      folds_[i].rel_org[1] = scale[1] * depend_fold.v;
    }
  }

  Fold fold = depend_fold;
  fold.org[0] -= depend_fold.rel_org[0];

  rects_ = vector<Rect3D>
  {
    { vector<EVec3d>{ EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero() }},
    { vector<EVec3d>{ EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero(), EVec3d::Zero() }}
  };

  Update(fold, (J_PI / 2.0));

  for (int i = 0; i < 2; i++)
    folds_[i].org90 = folds_[i].org;
}





