//
// Created by Dandoh on 1/11/18.
//

#include "body.h"
void body::reset_force() {
  f_x = 0;
  f_y = 0;
}
void body::update(double dt) {
  v_x += dt * f_x / mass;
  v_y += dt * f_y / mass;
  x += dt * v_x;
  y += dt * v_y;
}
