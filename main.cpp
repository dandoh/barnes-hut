#include<iostream>
#include<mpi.h>
#include<limits.h>
#include<vector>
#include<time.h>
#include<assert.h>
#include<cmath>

#define DEBUG(x) std::cerr << #x << " = " << x << std::endl;
#define DEBUGA(a, b) std::cout << "[ "; for (auto x = (a); x != (b); x++) std::cout << *x << " ";\
  std::cout << "]" << std::endl

const double G = 6.67e-11;
const double Theta = 0.5;

double f_rand(double f_min, double f_max) {
  double f = (double) rand() / RAND_MAX;
  return f_min + f * (f_max - f_min);
}

struct body {
  double mass;
  double x, y;
  double v_x = 0, v_y = 0;
  double f_x = 0, f_y = 0;
  double r, g, b;
  friend std::ostream &operator<<(std::ostream &os, const body &body1) {
    os << "(mass: " << body1.mass << " x: " << body1.x << " y: " << body1.y << " f_x: " << body1.f_x << " f_y: "
       << body1.f_y << ")";
    return os;
  }

  body() {}
  body(double mass, double x, double y) : mass(mass), x(x), y(y) {}
  body(double mass, double x, double y, double v_x, double v_y, double f_x, double f_y)
      : mass(mass), x(x), y(y), v_x(v_x), v_y(v_y), f_x(f_x), f_y(f_y) {}

  void reset_force() {
    f_x = 0;
    f_y = 0;
  }

  void update(double dt) {
    v_x += dt * f_x / mass;
    v_y += dt * f_y / mass;
    x += dt * v_x;
    y += dt * v_y;
  }
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
  std::vector<quad> quads; // array of buffer

  tree_buffer(int n, double length) : n(n), length(length) {
    quads.push_back(quad(0, 0, length));
  }

  void insert(body b) {
    insert(0, b);
  }

  void update_force(body &b) {
    update_force(0, b);
  }

 private:
  void put_body(int i, body b) {
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

  void insert(int i, body &b) {
    if (b.x < -length / 2 || b.x > length / 2 || b.y > length / 2 || b.y < -length / 2) return;

    if (quads[i].mass == 0) {
      quads[i].mass = b.mass;
      quads[i].mass_x = b.x, quads[i].mass_y = b.y;
      quads[i].v_x = b.v_x, quads[i].v_y = b.v_y;

      quads[i].f_x = b.f_x, quads[i].f_y = b.f_y;
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

        put_body(i, body(quads[i].mass, quads[i].mass_x, quads[i].mass_y,
                         quads[i].v_x, quads[i].v_y,
                         quads[i].f_x, quads[i].f_y));
        put_body(i, b);
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

};

int main() {
  freopen("positions.txt", "w", stdout);
  freopen("inputs/galaxy4.txt", "r", stdin);
  srand(time(NULL));

  int initial_size;
  double radius;
  std::cin >> initial_size >> radius;

  double length = 2 * radius;
  int current_size = initial_size;

  body bb[initial_size + 5];

  for (int i = 1; i <= initial_size; i++) {
    double x, y, vx, vy, mass;
    double r, g, b;
    std::cin >> x >> y >> vx >> vy >> mass >> r >> g >> b;
    bb[i] = body(mass, x, y, vx, vy, 0, 0);
    bb[i].r = r / 255; bb[i].g = g / 255; bb[i].b = b / 255;
  }

  std::cout << length << "\n";
  double dt = 0.2;
  for (double t = 0; t < 100; t += dt) {
    int sz = 0;
    for (int i = 1; i <= initial_size; i++) {
      if (bb[i].x <= length / 2 && bb[i].x >= -length / 2
          && bb[i].y <= length / 2 && bb[i].y >= -length / 2) {
        bb[++sz] = bb[i];
      }
    }
    current_size = sz;


    for (int i = 1; i <= current_size; i++) {
      std::cout << bb[i].x << " " << bb[i].y << " "
                << bb[i].r << " " << bb[i].g << " " << bb[i].b
                << "\n";
    }
    std::cout << "\n";

    tree_buffer tree(current_size, length);
    for (int i = 1; i <= current_size; i++) {
      tree.insert(bb[i]);
    }
    for (int i = 1; i <= current_size; i++) {
      bb[i].reset_force();
      tree.update_force(bb[i]);
      bb[i].update(dt);
    }
  }

}