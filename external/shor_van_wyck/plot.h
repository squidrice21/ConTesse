#ifndef PLOT
#define PLOT

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>

void plot_polygon(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXi &H,
                  const Eigen::MatrixXd &poly, int start_idx = -1,
                  int end_idx = -1);

void plot_mesh(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd &V,
               const Eigen::MatrixXi &F, const std::vector<int> &id,
               const Eigen::VectorXi &L, // hight faces
               bool show_boundary = false);

void plot_polygons(igl::opengl::glfw::Viewer &viewer,
                   std::vector<Eigen::Vector3d> const &poly1,
                   std::vector<Eigen::Vector3d> const &poly2, int f2_idx);

enum class OTYPE { MESH, POLYGON };

class Object {

public:
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  OTYPE type;
  OTYPE get_type() { return type; }
  int get_v_size() { return V.rows(); }
  int get_f_size() { return F.rows(); }
  Eigen::MatrixXd get_vertices() { return V; }
  Object(Eigen::MatrixXd V0, Eigen::MatrixXi F0, OTYPE t0)
      : V(V0), F(F0), type(t0){};
  ~Object(){};
};

void plots(std::vector<Object> &data);

#endif