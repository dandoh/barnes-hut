#include<iostream>
#include<mpi.h>
#include<limits.h>
#include<vector>
#include<time.h>
#include<assert.h>
#include<cmath>
#include<fstream>
#include<getopt.h>
#include<stddef.h>

#include "body.h"
#include "quad.h"
#include "tree_buffer.h"

#define DEBUG(x) std::cerr << #x << " = " << x << std::endl;
#define DEBUGA(a, b) std::cout << "[ "; for (auto x = (a); x != (b); x++) std::cout << *x << " ";\
  std::cout << "]" << std::endl


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