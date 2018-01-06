#include<iostream>
#include<mpi.h>
#include<limits.h>
#include<vector>
#include<assert.h>
#include <cmath>
#define DEBUG(x) std::cout << #x << " = " << x << std::endl;
#define DEBUGA(a, b) std::cout << "[ "; for (auto x = (a); x != (b); x++) std::cout << *x << " ";\
  std::cout << "]" << std::endl


const double G = 6.67e-11;
const double Theta = 0.5;

struct body {
  double mass;
  double x, y;
  double v_x = 0, v_y = 0;
  double f_x = 0, f_y = 0;

  body(double mass, double x, double y) : mass(mass), x(x), y(y) {}
  body(double mass, double x, double y, double v_x, double v_y, double f_x, double f_y)
      : mass(mass), x(x), y(y), v_x(v_x), v_y(v_y), f_x(f_x), f_y(f_y) {}
};

struct quad {

  double mass = 0;
  double mass_x = LLONG_MAX;
  double mass_y = LLONG_MAX;

  double center_x;
  double center_y;
  double length;

  int nw = -1, ne = -1, sw = -1, se = -1;

  // use only when it is leaf (i.e with one body)
  double v_x, v_y, f_x, f_y;

  quad() {}
  quad(double center_x, double center_y, double length) : center_x(center_x), center_y(center_y), length(length) {}

  // update mass as well as it's representative mass X and Y
  void update_mass_and_pos(body b) {
    mass_x = (mass_x * mass + b.x * b.mass) / (b.mass + mass);
    mass_y = (mass_y * mass + b.y * b.mass) / (b.mass + mass);
    mass += b.mass;
  }
  bool is_external() {
    return nw == -1 && ne == -1 && sw == -1 && se == -1;
  }

  double distance_to(body &b) {
    double dx = mass_x - b.x;
    double dy = mass_y - b.y;
    return sqrt(dx * dx + dy * dy);
  }

};

// wrapper around quad_tree so we can have an array to send over MPI
class tree_buffer {
 public:
  int n; // number of body
  double length; // length of outter most quad
  int next;
  quad *quads; // array of buffer

  tree_buffer(int n, double length) : n(n), length(length) {
    quads = new quad[8 * n];
    quads[0] = quad(0, 0, length);
    next = 1;
  }
  ~tree_buffer() {
    delete[] quads;
  }

  void insert(int i, body &b) {
    assert(i < next);

    if (quads[i].mass == 0) {
      quads[i].mass = b.mass;
      quads[i].mass_x = b.x, quads[i].mass_y = b.y;
      quads[i].v_x = b.v_x, quads[i].v_y = b.v_y;
      quads[i].f_x = b.f_x, quads[i].f_y = b.f_y;
    } else {
      // create 4 children
      if (quads[i].is_external()) {
        int nw = next++, ne = next++, sw = next++, se = next++;
        quads[i].nw = nw, quads[i].ne = ne, quads[i].sw = sw, quads[i].se = se;

        quads[nw] = quad(quads[i].center_x - quads[i].length / 4, quads[i].center_y + quads[i].length / 4, length / 2);
        quads[ne] = quad(quads[i].center_x + quads[i].length / 4, quads[i].center_y + quads[i].length / 4, length / 2);
        quads[sw] = quad(quads[i].center_x - quads[i].length / 4, quads[i].center_y - quads[i].length / 4, length / 2);
        quads[nw] = quad(quads[i].center_x + quads[i].length / 4, quads[i].center_y - quads[i].length / 4, length / 2);

        put_body(i, b);
        put_body(i, body(quads[i].mass, quads[i].mass_x, quads[i].mass_y,
                         quads[i].v_x, quads[i].v_y,
                         quads[i].f_x, quads[i].f_y));
        quads[i].update_mass_and_pos(b);
      } else {
        quads[i].update_mass_and_pos(b);
        put_body(i, b);
      }
    }
  }

  void update_force(int i, body &b) {
    if (quads[i].mass == 0) return;

    double d = quads[i].distance_to(b);
    if (quads[i].is_external() || (quads[i].length / d) < Theta) {
      double EPS = 3e4;
      double dx = quads[i].mass_x - b.x;
      double dy = quads[i].mass_y - b.y;
      double dist = sqrt(dx * dx + dy * dy);
      double F = G * b.mass * quads[i].mass / (dist * dist + EPS * EPS);

      b.f_x += F * dx / dist;
      b.f_y += F * dy / dist;

    } else {
      update_force(quads[i].nw, b);
      update_force(quads[i].ne, b);
      update_force(quads[i].sw, b);
      update_force(quads[i].se, b);
    }
  }


  void insert(body b) {
    insert(0, b);
  }

  void update_force(body &b) {
    update_force(0, b);
  }

 private:
  void put_body(int i, body b) {
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

};

int main() {
  tree_buffer tree(10, 16);
  tree.insert(body(12, 1, 2));
  tree.insert(body(6, -3, -5));
  tree.insert(body(6, 1, 1));
  DEBUG(tree.quads[0].mass_x);
}