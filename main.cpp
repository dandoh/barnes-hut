#include<iostream>
#include<mpi.h>
#include<limits.h>
#include<vector>
#include<time.h>
#include<assert.h>
#include<cmath>
#include<fstream>
#include<getopt.h>

#define DEBUG(x) std::cerr << #x << " = " << x << std::endl;
#define DEBUGA(a, b) std::cout << "[ "; for (auto x = (a); x != (b); x++) std::cout << *x << " ";\
  std::cout << "]" << std::endl

// Newton constant
const double G = 6.67e-11;
// Theta to approximate force
const double Theta = 0.5;

double f_rand(double f_min, double f_max) {
  double f = (double) rand() / RAND_MAX;
  return f_min + f * (f_max - f_min);
}

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
  void reset_force() {
    f_x = 0;
    f_y = 0;
  }

  /**
   * Update position of this body
   * @param dt delta time
   */
  void update(double dt) {
    v_x += dt * f_x / mass;
    v_y += dt * f_y / mass;
    x += dt * v_x;
    y += dt * v_y;
  }
};

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
  void update_mass_and_pos(body b) {
    mass_x = (mass_x * mass + b.x * b.mass) / (b.mass + mass);
    mass_y = (mass_y * mass + b.y * b.mass) / (b.mass + mass);
    mass += b.mass;
  }

  // not having any children
  bool is_external() {
    return nw == -1 && ne == -1 && sw == -1 && se == -1;
  }

  // distance to a body
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
  tree_buffer(int n, double length, quad *qs, int quads_length) : n(n), length(length) {
    for (int i = 0; i < quads_length; i++) {
      quads.push_back(qs[i]);
    }
  }

  /**
   * Insert a body into this quad tree
   * @param b
   */
  void insert(body b) {
    insert(0, b);
  }

  /**
   * Update a force to a body by this quad tree
   * @param b
   */
  void update_force(body &b) {
    update_force(0, b);
  }

 private:
  // put the body into appropriate quadrant
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

MPI_Datatype mpi_body_type;
MPI_Datatype mpi_quad_type;

std::string
    input_file,
    output_file = "positions.txt",
    mode = "p";
double dt = 0.2;
int num_steps = 250;
double elapsed_time;

void run_serial() {
  int n;
  // length of the field
  double length;
  // radius of the field, length = 2 * radius
  double radius;
  std::ifstream is;
  std::ofstream os;
  is.open(input_file);
  os.open(output_file);

  is >> n >> radius;
  length = 2 * radius;
  os << length << "\n";
  body *bodies = new body[n + 5];

  // read body from file
  for (int i = 0; i < n; i++) {
    double x, y, vx, vy, mass;
    double r, g, b;
    is >> x >> y >> vx >> vy >> mass >> r >> g >> b;
    bodies[i] = body(mass, x, y, vx, vy, 0, 0);
    bodies[i].r = r / 255;
    bodies[i].g = g / 255;
    bodies[i].b = b / 255;
  }
  is.close();

  for (int st = 0; st < num_steps; st++) {

    for (int i = 0; i < n; i++) {
      os << bodies[i].x << " " << bodies[i].y << " "
         << bodies[i].r << " " << bodies[i].g << " " << bodies[i].b
         << "\n";
    }
    os << "\n";

    tree_buffer tree(n, length);
    for (int i = 0; i < n; i++) {
      tree.insert(bodies[i]);
    }
    for (int i = 0; i < n; i++) {
      bodies[i].reset_force();
      tree.update_force(bodies[i]);
      bodies[i].update(dt);
    }
  }

}

void init_body_type() {
  int block_lengths[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype types[10] = {
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE
  };
  MPI_Aint offsets[10] = {
      offsetof(body, mass),
      offsetof(body, x), offsetof(body, y),
      offsetof(body, v_x), offsetof(body, v_y),
      offsetof(body, f_x), offsetof(body, f_y),
      offsetof(body, r), offsetof(body, g), offsetof(body, b)
  };
  MPI_Type_create_struct(10, block_lengths, offsets, types, &mpi_body_type);
  MPI_Type_commit(&mpi_body_type);
}

void init_quad_type() {
  int block_lengths[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype types[10] = {
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT
  };
  MPI_Aint offsets[10] = {
      offsetof(quad, mass),
      offsetof(quad, mass_x), offsetof(quad, mass_y),
      offsetof(quad, center_x), offsetof(quad, center_y),
      offsetof(quad, length),
      offsetof(quad, nw), offsetof(quad, ne), offsetof(quad, sw), offsetof(quad, se)
  };
  MPI_Type_create_struct(10, block_lengths, offsets, types, &mpi_quad_type);
  MPI_Type_commit(&mpi_quad_type);
}

void init_types() {
  init_body_type();
  init_quad_type();
}

void run_parallel() {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  init_types();

  // read file and broadcast number of bodies and field length
  // number of bodies
  int n;
  // length of the field
  double length;
  // radius of the field, length = 2 * radius
  double radius;
  std::ifstream is;
  std::ofstream os;
  if (rank == 0) {
    is.open(input_file);
    os.open(output_file);

    is >> n >> radius;
    length = 2 * radius;
    os << length << "\n";
  }

  // broadcast n and length
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  body *bodies = new body[n + 5];
  // read body from file
  if (rank == 0) {
    for (int i = 0; i < n; i++) {
      double x, y, vx, vy, mass;
      double r, g, b;
      is >> x >> y >> vx >> vy >> mass >> r >> g >> b;
      bodies[i] = body(mass, x, y, vx, vy, 0, 0);
      bodies[i].r = r / 255;
      bodies[i].g = g / 255;
      bodies[i].b = b / 255;
    }
    is.close();
  }

  body *local_bodies = new body[2 * n / size];

  // parameter for scattering the bodies over
  int *send_counts = new int[size];
  int *displs = new int[size];
  int rem = n % size;
  int sum = 0;

  for (int i = 0; i < size; i++) {
    send_counts[i] = n / size;
    if (rem > 0) {
      send_counts[i]++;
      rem--;
    }

    displs[i] = sum;
    sum += send_counts[i];
  }

  for (int st = 0; st < num_steps; st++) {
    // print position of bodies
    if (rank == 0) {
      for (int i = 0; i < n; i++) {
        os << bodies[i].x << " " << bodies[i].y << " "
           << bodies[i].r << " " << bodies[i].g << " " << bodies[i].b
           << "\n";
      }
      os << "\n";
    }
    // build the tree in root process
    int quads_length = 0;
    tree_buffer tree(n, length);

    if (rank == 0) {
      for (int i = 0; i < n; i++) {
        tree.insert(bodies[i]);
      }
      quads_length = (int) tree.quads.size();
    }
    // send info about the array of quads inside the tree
    MPI_Bcast(&quads_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    quad *qs = new quad[quads_length];

    if (rank == 0) {
      std::copy(tree.quads.begin(), tree.quads.end(), qs);
    }

    // send the tree over the network - TODO
    MPI_Bcast(qs, quads_length, mpi_quad_type, 0, MPI_COMM_WORLD);
    tree_buffer local_tree(n, length, qs, quads_length);

    // scatter the body
    MPI_Scatterv(bodies, send_counts, displs, mpi_body_type, local_bodies,
                 send_counts[rank], mpi_body_type, 0, MPI_COMM_WORLD);
    // compute
    for (int i = 0; i < send_counts[rank]; i++) {
      local_bodies[i].reset_force();
      local_tree.update_force(local_bodies[i]);
      local_bodies[i].update(dt);
    }

    // gather body
    MPI_Gatherv(local_bodies, send_counts[rank], mpi_body_type,
                bodies, send_counts, displs, mpi_body_type, 0, MPI_COMM_WORLD);
  }

  delete[] bodies;
}

bool parse_input(int argc, char *argv[]) {
  static const char *opt_string = "i:o:m:d:s:";
  int in;
  try {
    while ((in = getopt(argc, argv, opt_string)) != -1) {
      switch (in) {
        case 'i':input_file = optarg;
          break;
        case 'o':output_file = optarg;
          break;
        case 'm':mode = optarg;
          break;
        case 'd':dt = std::stod(optarg);
          break;
        case 's':num_steps = std::stoi(optarg);
          break;
        default:return false;
      }
    }
  } catch (std::exception const &e) {
    return false;
  }

  // must provide input file
  if (input_file == "") {
    return false;
  }

  return true;
}

void print_statistics() {
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0) {
    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "Mode: " << (mode == "p" ? "Parallel" : "Sequential") << std::endl;
    if (mode == "p") {
      std::cout << "Num core: " << size << std::endl;
    }
    std::cout << "Delta time: " << dt << std::endl;
    std::cout << "Number of steps: " << num_steps << std::endl;
    std::cout << "Elapsed time: " << elapsed_time << "s" << std::endl;
  }
}

int main(int argc,
         char *argv[]) {

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // parse arguments
  if (!parse_input(argc, argv)) {
    if (rank == 0) {
      std::cout << "Error: Invalid command line input" << std::endl;
    }
    return -1;
  }

  double start_time = MPI_Wtime();
  if (mode == "p") {
    run_parallel();
  } else {
    if (rank == 0) {
      run_serial();
    }
  }
  elapsed_time = MPI_Wtime() - start_time;
  print_statistics();

  MPI_Finalize();
}