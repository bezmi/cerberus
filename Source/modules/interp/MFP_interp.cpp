#include "MFP_interp.H"
#include "MFP_global.H"
#include "MFP_source.H"

#include <fstream>
#include <sstream>

#include "delaunator.hpp"


//#define DBL_EPSILON 2.2204460492503131e-16
#define EPS 100 * DBL_EPSILON

Interp2D::Interp2D() {}

void Interp2D::init(const std::vector<double> xs, const std::vector<double> ys,
                    const std::vector<double> zs, const bool logy, const bool logx) {
  BL_PROFILE("Interp2D::init");
  std::vector<double> coords;

  vals = zs;

  for (size_t i = 0; i < xs.size(); i++) {
    if (logx) {
      coords.push_back(std::log10(xs[i]));
    } else {
      coords.push_back(xs[i]);
    }
    if (logy) {
      coords.push_back(std::log10(ys[i]));
    } else {
      coords.push_back(ys[i]);
    }
  }

  delaunator::Delaunator d(coords);

  d_triangles = d.triangles;
  d_coords = d.coords;
  for (size_t i = 0; i < d_coords.size(); i++) {
    if (logx and not (i%2)) {
      d_coords[i] = std::pow(10,d_coords[i]);
    } else if (logy and (i%2)) {
      d_coords[i] = std::pow(10,d_coords[i]);
    }
  }
  d_halfedges = d.halfedges;

  is_valid = true;
}

Interp2D::Interp2D(const std::vector<double> xs, const std::vector<double> ys,
                   const std::vector<double> zs, const bool logy, const bool logx) {
  init(xs, ys, zs, logy, logx);
}

Interp2D::Interp2D(const std::string ifile, bool logy, bool logx) {
  BL_PROFILE("Interp2D::Interp2D (read file)");
  std::vector<double> coords;
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> zs;

  std::ifstream input_file;
  input_file.open(ifile);

  std::string line;
  std::string x_st, y_st, z_st;

  while (std::getline(input_file, line, '\n')) {
    std::istringstream line_stream(line);
    getline(line_stream, x_st, ',');
    xs.push_back(stod(x_st));

    getline(line_stream, y_st, ',');
    ys.push_back(stod(y_st));

    getline(line_stream, z_st, ',');
    zs.push_back(stod(z_st));
  }
  x_min = xs[0];
  x_max = xs[0];
  y_min = ys[0];
  y_max = ys[0];

  for (size_t i = 0; i < xs.size(); i++) {
    if (xs[i] < x_min) {
      x_min = xs[i];
    }
    if (ys[i] < y_min) {
      y_min = ys[i];
    }
    if (xs[i] > x_max) {
      x_max = xs[i];
    }
    if (ys[i] > y_max) {
      y_max = ys[i];
    }
  }

  input_file.close();

  init(xs, ys, zs, logy, logx);
}

std::array<size_t, 3> Interp2D::edges_of_triangle(size_t t) const {
  return {3 * t, 3 * t + 1, 3 * t + 2};
}

size_t Interp2D::triangle_of_edge(size_t e) const { return floor(e / 3); }

std::vector<size_t> Interp2D::triangles_adjacent_to_triangle(size_t t) const {
  std::vector<size_t> adjacent_triangles;
  for (auto e : edges_of_triangle(t)) {
    std::size_t opposite = d_halfedges[e];
    if (opposite >= 0) {
      adjacent_triangles.push_back(triangle_of_edge(e));
    }
  }
  return adjacent_triangles;
}


double Interp2D::distance_along_z(double x1, double y1, double z1, double x2,
                                  double y2, double z2, double x3, double y3,
                                  double z3, double x, double y) const {
  BL_PROFILE("Interp2D::distance_along_z");
  double detT = ((y2 - y3) * (x1 - x3)) + ((x3 - x2) * (y1 - y3));
  double z = std::numeric_limits<double>::infinity();

  if (std::abs(detT) >= 1e-16) {
    double lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
    double lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
    double lambda3 = 1 - lambda2 - lambda1;

    if (lambda1 >= 0 - EPS && lambda1 <= 1 + EPS && lambda2 >= 0 - EPS &&
        lambda2 <= 1 + EPS && lambda3 >= 0 - EPS && lambda3 <= 1 + EPS) {
      z = (lambda1 * z1) + (lambda2 * z2) + (lambda3 * z3);
    }
  }
  return z;
}

double Interp2D::point_in_triel(size_t triel, double x, double y) const {
  BL_PROFILE("Interp2D::point_in_triel");
  double d;
  int i, j, ztest = 0;
  double xi, yi;
  double xj, yj;

  for (int i = 2, j = i - 1; i > -1; i--, j--) {
    if (j == -1) j = 3;

    auto ijpoint = d_triangles[triel + j];
    auto iipoint = d_triangles[triel + i];

    xj = d_coords[2 * ijpoint];
    yj = d_coords[2 * ijpoint + 1];

    xi = d_coords[2 * iipoint];
    yi = d_coords[2 * iipoint + 1];

    d = (xj - xi) * (y - yi) - (yj - yi) * (x - xi);
    if (d < 0.0) return -1;
    if (d == 0.0) ztest = 1;
  }
  return (ztest ? 0 : 1);
}

double Interp2D::interpolate(double x, double y)  {
  BL_PROFILE("Interp2D::interpolate");
  // compute the linear barycentric interpolation for x and y
  // Finds the simplex containing x and y via Lawson's oriented walk

  // TODO: change this to the Barycentric Walk algorithm, which should be
  // able to combine the walk and interpolation steps and save on a call to the
  // distance_along_z function.

  // use nearest co-ordinate for out of bounds
  if (x < x_min) {
    x = x_min;
  } else if (x > x_max) {
    x = x_max;
  }

  if (y < y_min) {
    y = y_min;
  } else if (y > y_max) {
    y = y_max;
  }

  size_t triel = last_triel;
  bool found = false;

  while (!found) {
    found = true;

    std::array<std::size_t, 3> edges = edges_of_triangle(triel);

    // loop over edges of the triangle
    for (int e = 0; e < 3; e++) {

      // d_triangles[e] is the index of the first point for edge e
      std::size_t l = d_triangles[edges[e]];
      double lx = d_coords[2 * l];
      double ly = d_coords[2 * l + 1];

      int e_end = (e + 1) % 3;
      std::size_t r = d_triangles[edges[e_end]];
      double rx = d_coords[2 * r];
      double ry = d_coords[2 * r + 1];

      // if this value is > 0, then move to the next triangle
      double or2d = (rx - lx) * (y - ly) - (x - lx) * (ry - ly);
      if (or2d >= 0 + EPS) {
        triel = triangle_of_edge(d_halfedges[edges[e]]);
        found = false;
        break;
      }
    }
  }

  // save the last triangle so we can use it as our starting point next
  last_triel = triel;

  auto ai = d_triangles[3*triel];
  auto bi = d_triangles[3*triel + 1];
  auto ci = d_triangles[3*triel + 2];

  double ax = d_coords[2 * ai];
  double ay = d_coords[2 * ai + 1];
  double az = vals[ai];

  double bx = d_coords[2 * bi];
  double by = d_coords[2 * bi + 1];
  double bz = vals[bi];

  double cx = d_coords[2 * ci];
  double cy = d_coords[2 * ci + 1];
  double cz = vals[ci];

  double z1 = distance_along_z(ax, ay, az, bx, by, bz, cx, cy, cz, x, y);
  return z1;
}
