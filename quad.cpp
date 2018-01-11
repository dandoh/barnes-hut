//
// Created by Dandoh on 1/11/18.
//

#include "quad.h"
void quad::update_mass_and_pos(body b) {
  mass_x = (mass_x * mass + b.x * b.mass) / (b.mass + mass);
  mass_y = (mass_y * mass + b.y * b.mass) / (b.mass + mass);
  mass += b.mass;
}
bool quad::is_external() {
  return nw == -1 && ne == -1 && sw == -1 && se == -1;
}
double quad::distance_to(body &b) {
  double dx = mass_x - b.x;
  double dy = mass_y - b.y;
  return sqrt(dx * dx + dy * dy);
}
