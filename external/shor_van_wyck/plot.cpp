#include "plot.h"
#include <igl/boundary_loop.h>
#include <unordered_set>
#include <vector>

void plot_polygon(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXi &H,
                  const Eigen::MatrixXd &poly, int start_idx, int end_idx) {
  viewer.data().clear();
  viewer.core().align_camera_center(poly);

  std::unordered_set<int> highlight_indices;
  if (start_idx >= 0 && end_idx >= 0) {
    for (int i = start_idx; i != end_idx; i = (i + 1) % poly.rows()) {
      highlight_indices.emplace(i);
    }
  }

  for (int i = 0; i < poly.rows(); i++) {
    int i_1 = (i + 1) % poly.rows();
    viewer.data().add_label(poly.row(i), std::to_string(i));

    if (!highlight_indices.count(i))
      viewer.data().add_edges(poly.row(i), poly.row(i_1),
                              Eigen::RowVector3d(1, 0, 0));
    else
      viewer.data().add_edges(poly.row(i), poly.row(i_1),
                              Eigen::RowVector3d(0, 1, 0));
  }
  for (int i = 0; i < H.rows(); i++) {
    if (H(i))
      viewer.data().add_points(poly.row(i), Eigen::RowVector3d(0, 0, 0));
  }
}

void plot_mesh(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd &V,
               const Eigen::MatrixXi &F, const std::vector<int> &id,
               const Eigen::VectorXi &L, // hight faces
               bool show_boundary) {
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.core().align_camera_center(V, F);
  Eigen::MatrixXd C(F.rows(), 3);
  C.setConstant(1);
  if (show_boundary) {
    Eigen::VectorXi bd;
    igl::boundary_loop(F, bd);
    for (int i = 0; i < bd.rows(); i++) {
      int i_1 = (i + 1) % bd.rows();
      viewer.data().add_edges(V.row(bd(i)), V.row(bd(i_1)),
                              Eigen::RowVector3d(0, 1, 0));
    }
  }
  for (int i : id) {
    viewer.data().add_points(V.row(i), Eigen::RowVector3d(1, 0, 0));
  }
  for (int i = 0; i < L.rows(); i++) {
    if (L(i) != 0)
      C.row(i) << 1, 0, 0;
  }
  if (L.sum() > 0)
    viewer.data().set_colors(C);
}

void plot_polygons(igl::opengl::glfw::Viewer &viewer,
                   std::vector<Eigen::Vector3d> const &poly1,
                   std::vector<Eigen::Vector3d> const &poly2, int f2_idx) {
  viewer.data().clear();
  Eigen::MatrixXd poly2_mat;
  poly2_mat.resize(3, 3);
  poly2_mat.row(0) = poly2[0];
  poly2_mat.row(1) = poly2[1];
  poly2_mat.row(2) = poly2[2];

  viewer.core().align_camera_center(poly2_mat);
  for (size_t i = 0; i < poly1.size(); i++) {
    size_t i_1 = (i + 1) % poly1.size();
    viewer.data().add_edges(poly1[i].transpose(), poly1[i_1].transpose(),
                            Eigen::RowVector3d(1, 0, 0));
  }
  for (size_t i = 0; i < poly2.size(); i++) {
    size_t i_1 = (i + 1) % poly2.size();
    viewer.data().add_label(poly2[i], std::to_string(f2_idx));
    viewer.data().add_edges(poly2[i].transpose(), poly2[i_1].transpose(),
                            Eigen::RowVector3d(1, 0, 0));
  }
}