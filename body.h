//
// Created by Dandoh on 1/11/18.
//

#ifndef BTL_BODY_H
#define BTL_BODY_H

#include <iostream>
/**
 * Struct represent a body
 */
struct body {
  double mass;
  // position
  double x, y;
  // velocity
  double v_x = 0, v_y = 0;
  // force
  double f_x = 0, f_y = 0;
  double r, g, b; // color to draw
  friend std::ostream &operator<<(std::ostream &os, const body &body1) {
    os << "(mass: " << body1.mass << " x: " << body1.x << " y: " << body1.y << " f_x: " << body1.f_x << " f_y: "
       << body1.f_y << ")";
    return os;
  }

  body() {}
  body(double mass, double x, double y) : mass(mass), x(x), y(y) {}
  body(double mass, double x, double y, double v_x, double v_y, double f_x, double f_y)
      : mass(mass), x(x), y(y), v_x(v_x), v_y(v_y), f_x(f_x), f_y(f_y) {}

  /**
   * Reset all the force affected on this body
   */
  void reset_force();

  /**
   * Update position of this body
   * @param dt delta time
   */
  void update(double dt);
};

#endif //BTL_BODY_H
