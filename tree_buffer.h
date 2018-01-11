//
// Created by Dandoh on 1/11/18.
//

#ifndef BTL_TREE_BUFFER_H
#define BTL_TREE_BUFFER_H

#include <vector>
#include <cassert>
#include "quad.h"
#include "body.h"

// Newton constant
const double G = 6.67e-11;
// Theta to approximate force
const double Theta = 0.5;

// wrapper around quad_tree so we can have an array to send over MPI
class tree_buffer {
 public:
  int n; // number of body
  double length; // length of outter most quad
  std::vector<quad> quads; // array of buffer

  tree_buffer(int n, double length) : n(n), length(length) {
    quads.push_back(quad(0, 0, length));
  }
  tree_buffer(int n, double length, quad *qs, int quads_length) : n(n), length(length) {
    for (int i = 0; i < quads_length; i++) {
      quads.push_back(qs[i]);
    }
  }

  /**
   * Insert a body into this quad tree
   * @param b
   */
  void insert(body b);

  /**
   * Update a force to a body by this quad tree
   * @param b
   */

  void update_force(body &b);

 private:
  // put the body into appropriate quadrant
  void put_body(int i, body b);

  void insert(int i, body &b);

  void update_force(int i, body &b);
};

#endif //BTL_TREE_BUFFER_H
