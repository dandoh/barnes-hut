//
// Created by Dandoh on 1/11/18.
//

#ifndef BTL_QUAD_H
#define BTL_QUAD_H

#include <cmath>
#include "body.h"

struct quad {

  double mass = 0;
  // x mass and y mass represent for all bodies within this quad
  double mass_x = 0;
  double mass_y = 0;

  // from center_x|y - length / 2 to center_x|y + length / 2
  double center_x;
  double center_y;
  double length;

  // index of its child. default is -1
  int nw = -1, ne = -1, sw = -1, se = -1;

  quad() {}
  quad(double center_x, double center_y, double length) : center_x(center_x), center_y(center_y), length(length) {}

  // update mass as well as it's representative mass X and Y
  void update_mass_and_pos(body b);

  // not having any children
  bool is_external();

  // distance to a body
  double distance_to(body &b);

};

#endif //BTL_QUAD_H
