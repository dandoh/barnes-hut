//
// Created by Dandoh on 1/11/18.
//

#include "tree_buffer.h"
void tree_buffer::insert(body b) {
  insert(0, b);
}
void tree_buffer::update_force(body &b) {
  update_force(0, b);
}
void tree_buffer::put_body(int i, body b) {
  assert(quads[i].nw != -1 && quads[i].ne != -1 && quads[i].sw != -1 && quads[i].se != -1);
  assert(quads[i].nw != i && quads[i].ne != i && quads[i].sw != i && quads[i].se != i);

  if (b.x <= quads[i].center_x && b.y >= quads[i].center_y) {
    insert(quads[i].nw, b);
  } else if (b.x >= quads[i].center_x && b.y >= quads[i].center_y) {
    insert(quads[i].ne, b);
  } else if (b.x <= quads[i].center_x && b.y <= quads[i].center_y) {
    insert(quads[i].sw, b);
  } else {
    insert(quads[i].se, b);
  }

}
void tree_buffer::insert(int i, body &b) {
  // make sure the body lie within the quadrant
  if (b.x < -length / 2 || b.x > length / 2 || b.y > length / 2 || b.y < -length / 2) return;

  if (quads[i].mass == 0) {
    quads[i].mass = b.mass;
    quads[i].mass_x = b.x, quads[i].mass_y = b.y;
  } else {
    // create 4 children
    if (quads[i].is_external()) {
      quads.push_back(quad(quads[i].center_x - quads[i].length / 4,
                           quads[i].center_y + quads[i].length / 4,
                           quads[i].length / 2));
      quads[i].nw = (int) (quads.size() - 1);

      quads.push_back(quad(quads[i].center_x + quads[i].length / 4,
                           quads[i].center_y + quads[i].length / 4,
                           quads[i].length / 2));
      quads[i].ne = (int) (quads.size() - 1);

      quads.push_back(quad(quads[i].center_x - quads[i].length / 4,
                           quads[i].center_y - quads[i].length / 4,
                           quads[i].length / 2));
      quads[i].sw = (int) (quads.size() - 1);

      quads.push_back(quad(quads[i].center_x + quads[i].length / 4,
                           quads[i].center_y - quads[i].length / 4,
                           quads[i].length / 2));
      quads[i].se = (int) (quads.size() - 1);
      // done with this node, now children
      put_body(i, body(quads[i].mass, quads[i].mass_x, quads[i].mass_y));
      put_body(i, b);
      // update the mass_x, mass_y
      quads[i].update_mass_and_pos(b);
    } else {
      quads[i].update_mass_and_pos(b);
      put_body(i, b);
    }
  }

}
void tree_buffer::update_force(int i, body &b) {
  if (quads[i].mass == 0) return;

  double d = quads[i].distance_to(b);
  if (quads[i].is_external() || (quads[i].length / d) < Theta) {
    double dx = quads[i].mass_x - b.x;
    double dy = quads[i].mass_y - b.y;
    double dist = sqrt(dx * dx + dy * dy);
    double EPS = 3e4;

    if (dist == 0) {
      return;
    }
    double F = (G * b.mass * quads[i].mass) / (dist * dist + EPS * EPS);

    b.f_x += F * dx / dist;
    b.f_y += F * dy / dist;

  } else {
    update_force(quads[i].nw, b);
    update_force(quads[i].ne, b);
    update_force(quads[i].sw, b);
    update_force(quads[i].se, b);
  }

}
