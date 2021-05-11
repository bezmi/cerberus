#ifndef MFP_INTERP_H
#define MFP_INTERP_H

#include <iostream>
#include <vector>
#include <array>

class Interp2D {
 public:
  Interp2D();
  Interp2D(const std::vector<double> xs, const std::vector<double> ys,
           const std::vector<double> zs, const bool logy = true, const bool logx = true);
  Interp2D(const std::string ifile, bool logy = true, bool logx = true);
  double distance_along_z(double x1, double y1, double z1, double x2, double y2,
                          double z2, double x3, double y3, double z3, double x,
                          double y) const;

  double interpolate(double x, double y);
  size_t triangle_of_edge(size_t e) const;
  std::vector<size_t> triangles_adjacent_to_triangle(size_t t) const;
  std::array<size_t, 3> edges_of_triangle(size_t t) const;
  double point_in_triel(size_t triel, double x, double y) const;
  double x_min = -1;
  double x_max = -1;
  double y_min = -1;
  double y_max = -1;
  std::vector<double> vals;
  std::vector<std::size_t> d_triangles;
  std::vector<double> d_coords;
  std::vector<std::size_t> d_halfedges;
  bool is_valid = false;

  std::size_t last_triel = 0;

 private:
  void init(const std::vector<double> xs, const std::vector<double> ys,
            const std::vector<double> zs, const bool logy, const bool logx);
};

#endif
