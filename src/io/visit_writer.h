// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include <cstdio>
#include <limits>

#include "common.h"

/*****************************************************************************
 *
 * Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
 * Produced at the Lawrence Livermore National Laboratory
 * LLNL-CODE-400142
 * All rights reserved.
 *
 * This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
 * full copyright notice is contained in the file COPYRIGHT located at the root
 * of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
 *
 * Redistribution  and  use  in  source  and  binary  forms,  with  or  without
 * modification, are permitted provided that the following conditions are met:
 *
 *  - Redistributions of  source code must  retain the above  copyright notice,
 *    this list of conditions and the disclaimer below.
 *  - Redistributions in binary form must reproduce the above copyright notice,
 *    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
 * ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
 * LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
 * DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
 * CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
 * LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
 * OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 *****************************************************************************/

enum {
 VISIT_VERTEX= 1,
 VISIT_LINE =3,
 VISIT_TRIANGLE= 5,
 VISIT_QUAD= 9,
 VISIT_TETRA= 10,
 VISIT_HEXAHEDRON= 12,
 VISIT_WEDGE= 13,
 VISIT_PYRAMID= 14
};

enum {
 VISIT_CELL_VAR = 0,
 VISIT_POINT_VAR = 1
};


class Visit_writer
{
public:
  /* ************************************************************************* //
  //                              visit_writer.h                               //
  // ************************************************************************* */

  /*
  // This file contains function prototypes for writing out point meshes,
  // unstructured meshes, rectilinear meshes, regular meshes, and
  // structured/curvilinear meshes into files that can later be read by VisIt.
  //
  // Each routine assumes that the data being written is three-dimensional.
  // If the data is two-dimensional, you must still write out the data
  // as three-dimensional (ie pad arrays so that they are the correct size, etc).
  // However: the VisIt reader will determine that the data is truly two-
  // dimensional and visualize it as a two-dimensional dataset.
  //
  // All writers have an ASCII vs Binary decision.  The tradeoffs are the
  // standard ones: ASCII is human readable, but slow.  The
  // binary is much faster, but not human readable.  Note: the binary format
  // is portable, since it converts all data to be big-endian (this was a
  // design decision for the format the visit_writer writes to -- the VTK
  // format).
  //
  // If you have multiple grids, you can write out one file for each grid.
  // There are potential pitfalls in doing this, where extra geometry and
  // interpolation problems appear along grid boundaries.  For additional
  // help with this issue, e-mail visit-users@ornl.gov
  */


  /* ****************************************************************************
  //  Function: write_point_mesh
  //
  //  Purpose:
  //      Writes out a point mesh.
  //
  //  Arguments:
  //      filename   The name of the file to write.  If the extension ".vtk" is
  //                 not present, it will be added.
  //      useBinary  '0' to write ASCII, !0 to write binary
  //      npts       The number of points in the mesh.
  //      pts        The spatial locations of the points.  This array should
  //                 be size 3*npts.  The points should be encoded as:
  //                 <x1, y1, z1, x2, y2, z2, ..., xn, yn, zn>
  //      nvars      The number of variables.
  //      vardim     The dimension of each variable.  The size of vardim should
  //                 be nvars.  If var i is a scalar, then vardim[i] = 1.
  //                 If var i is a vector, then vardim[i] = 3.
  //      vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  //
  //  Programmer: Hank Childs
  //  Creation:   September 2, 2004
  //
  // ***************************************************************************/

  void write_point_mesh(
    const char * filename,
    int useBinary,
    int npts,
    float * pts,
    int nvars,
    int * vardim,
    const char * const * varnames,
    float ** vars
  );



  /* ****************************************************************************
  //  Function: write_unstructured_mesh
  //
  //  Purpose:
  //      Writes out a unstructured mesh.
  //
  //
  //  Arguments:
  //      filename   The name of the file to write.  If the extension ".vtk" is
  //                 not present, it will be added.
  //      useBinary  '0' to write ASCII, !0 to write binary
  //      npts       The number of points in the mesh.
  //      pts        The spatial locations of the points.  This array should
  //                 be size 3*npts.  The points should be encoded as:
  //                 <x1, y1, z1, x2, y2, z2, ..., xn, yn, zn>
  //      ncells     The number of cells.
  //      celltypes  The type of each cell.
  //      conn       The connectivity array.
  //      nvars      The number of variables.
  //      vardim     The dimension of each variable.  The size of vardim should
  //                 be nvars.  If var i is a scalar, then vardim[i] = 1.
  //                 If var i is a vector, then vardim[i] = 3.
  //      centering  The centering of each variable.  The size of centering
  //                 should be nvars.  If centering[i] == 0, then the variable
  //                 is cell-based.  If centering[i] != 0, then the variable
  //                 is point-based.
  //      vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  //
  //  Example:
  //      You have two triangles.  The first has points (0,0,0), (0,1,0), and
  //      (1,1,0).  The second has points (0,0,0), (1,1,0), and (1,0,0).
  //
  //      There are four unique points.
  //
  //      float pts[12] = { 0,0,0, 0,1,0, 1,1,0, 1,0,0 };
  //
  //      It is important the points list contain only unique points,
  //      because VisIt is not able to correctly determine the connectivity of a
  //      dataset when points are duplicated.
  //
  //      There are two triangles.
  //      int ncells = 2;
  //
  //      The cells are both triangles.
  //      int celltypes[2] = { VISIT_TRIANGLE, VISIT_TRIANGLE };
  //
  //      The connectivity contains indices into the points list.  The indexing
  //      assumes that each point has size 3 (x,y,z).
  //
  //      int conn[6] = { 0, 1, 2, 0, 2, 3 };
  //
  //  Hint:
  //      When writing an unstructured mesh, it is easy to get the orientation
  //      of a cell backwards.  VisIt typically does okay with this, but it
  //      can cause problems.  To test if this is happening, bring up VisIt on
  //      your newly outputted dataset and make a Pseudocolor plot of
  //      "mesh_quality/volume" for 3D datasets or "mesh_quality/area" for 2D
  //      datasets.  If the cells are inside-out, the volumes or areas will be
  //      negative.
  //
  //
  //  Programmer: Hank Childs
  //  Creation:   September 2, 2004
  //
  // ***************************************************************************/

  void write_unstructured_mesh(
    const char * filename,
    int useBinary,
    int npts,
    float * pts,
    int ncells,
    int * celltypes,
    int * conn,
    int nvars,
    int * vardim,
    int * centering,
    const char * const * varnames,
    float ** vars
  );



  /* ****************************************************************************
  //  Function: write_regular_mesh
  //
  //  Purpose:
  //      Writes out a regular mesh.  A regular mesh is one where the data lies
  //      along regular intervals.  "Brick of bytes/floats",
  //      "Block of bytes/floats", and MRI data all are examples of data that
  //      lie on regular meshes.
  //
  //
  //  Arguments:
  //      filename   The name of the file to write.  If the extension ".vtk" is
  //                 not present, it will be added.
  //      useBinary  '0' to write ASCII, !0 to write binary
  //      dims       An array of size 3 = { nX, nY, nZ }, where nX is the
  //                 number of points in the X-dimension, etc.
  //      nvars      The number of variables.
  //      vardim     The dimension of each variable.  The size of vardim should
  //                 be nvars.  If var i is a scalar, then vardim[i] = 1.
  //                 If var i is a vector, then vardim[i] = 3.
  //      centering  The centering of each variable.  The size of centering
  //                 should be nvars.  If centering[i] == 0, then the variable
  //                 is cell-based.  If centering[i] != 0, then the variable
  //                 is point-based.
  //      vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  //
  //
  //  Programmer: Hank Childs
  //  Creation:   September 2, 2004
  //
  // ***************************************************************************/

  void write_regular_mesh(
    const char * filename,
    int useBinary,
    int * dims,
    int nvars,
    int * vardim,
    int * centering,
    const char * const * varnames,
    float ** vars
  );



  /* ****************************************************************************
  //  Function: write_rectilinear_mesh
  //
  //  Purpose:
  //      Writes out a rectilinear mesh.
  //
  //
  //  Arguments:
  //      filename   The name of the file to write.  If the extension ".vtk" is
  //                 not present, it will be added.
  //      useBinary  '0' to write ASCII, !0 to write binary
  //      dims       An array of size 3 = { nX, nY, nZ }, where nX is the
  //                 number of points in the X-dimension, etc.
  //      x          An array of size dims[0] that contains the x-coordinates.
  //      y          An array of size dims[1] that contains the x-coordinates.
  //      z          An array of size dims[2] that contains the x-coordinates.
  //      nvars      The number of variables.
  //      vardim     The dimension of each variable.  The size of vardim should
  //                 be nvars.  If var i is a scalar, then vardim[i] = 1.
  //                 If var i is a vector, then vardim[i] = 3.
  //      centering  The centering of each variable.  The size of centering
  //                 should be nvars.  If centering[i] == 0, then the variable
  //                 is cell-based.  If centering[i] != 0, then the variable
  //                 is point-based.
  //      vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  //
  //
  //  Example:
  //      You have a rectilinear mesh with x = { 0, 1, 2}, y = { 1, 1.5, 2, 3 },
  //      and z = { 2.5, 3.5 }.
  //
  //      Then dims = { 3, 4, 2 }.
  //
  //  Programmer: Hank Childs
  //  Creation:   September 2, 2004
  //
  // ***************************************************************************/

  void write_rectilinear_mesh(
    const char * filename,
    int useBinary,
    int * dims,
    float * x,
    float * y,
    float * z,
    int nvars,
    int * vardim,
    int * centering,
    const char * const * varnames,
    float ** vars
  );



  /* ****************************************************************************
  //  Function: write_curvilinear_mesh
  //
  //  Purpose:
  //      Writes out a curvilinear mesh.
  //
  //
  //  Arguments:
  //      filename   The name of the file to write.  If the extension ".vtk" is
  //                 not present, it will be added.
  //      useBinary  '0' to write ASCII, !0 to write binary
  //      dims       An array of size 3 = { nI, nJ, nK }, where nI is the
  //                 number of points in the logical I dimension, etc.
  //      pts        An array of size nI*nJ*nK*3.  The array should be layed
  //                 out as (pt(i=0,j=0,k=0), pt(i=1,j=0,k=0), ...
  //                 pt(i=nI-1,j=0,k=0), pt(i=0,j=1,k=0), ...).
  //      nvars      The number of variables.
  //      vardim     The dimension of each variable.  The size of vardim should
  //                 be nvars.  If var i is a scalar, then vardim[i] = 1.
  //                 If var i is a vector, then vardim[i] = 3.
  //      centering  The centering of each variable.  The size of centering
  //                 should be nvars.  If centering[i] == 0, then the variable
  //                 is cell-based.  If centering[i] != 0, then the variable
  //                 is point-based.
  //      vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  //
  //
  //  Programmer: Hank Childs
  //  Creation:   September 2, 2004
  //
  // ***************************************************************************/

  void write_curvilinear_mesh(
    const char * filename,
    int useBinary,
    int * dims,
    float * pts,
    int nvars,
    int * vardim,
    int * centering,
    const char * const * varnames,
    float ** vars
  );


private:
  FILE * fp = NULL;
  int useBinary = 0;
  int numInColumn = 0;

  void end_line(void);
  void open_file(const char * filename);
  void close_file(void);
  void force_big_endian(unsigned char * bytes);
  void write_string(const char * str);
  void new_section(void);
  void write_int(int val);
  void write_float(float val);
  void write_header(void);
  void write_variables(
    int nvars,
    int * vardim,
    int * centering,
    const char * const * varname,
    float ** vars,
    int npts,
    int ncells
  );


};

// Helper to get better auto complete.
struct Visit_args {
  // Outoput file name
  const char * filename = nullptr;
  // write in ascii or binary
  int useBinary = 1;
  // Number of vertices
  int npts = std::numeric_limits<int>::max();
  // Vertices
  float * pts = nullptr;
  // Number of cells (triangles, lines, quads are all cells
  int ncells = std::numeric_limits<int>::max();
  // An array that for each cell contains its type
  // The ones we need are  VISIT_LINE, VISIT_TRIANGLE, VISIT_QUAD
  int * celltypes = nullptr;
  // Connectivity: a flat array that for each cell contains its point indices
  int * conn = nullptr;
  // Number of variables you want to write (variables can be visualized)
  int nvars = 0;
  // Dimension of each variable (can be 1 or 3)
  int * vardim = nullptr;
  // For each variable, does it belong to vertices or cells
  //centering[i] == VISIT_CELL_VAR-> per-cell variable
  //centering[i] == VISIT_POINT_VAR-> per-vertex variable
  int * centering = nullptr;
  // A name for each variable
  const char * const * varnames = nullptr;
  // Values of each variable
  float ** vars = nullptr;

  void write() {
    contess_assert(filename);
    contess_assert(useBinary != std::numeric_limits<int>::max());
    contess_assert(npts != std::numeric_limits<int>::max());
    contess_assert(pts);
    contess_assert(ncells != std::numeric_limits<int>::max());
    contess_assert(celltypes);
    contess_assert(conn);
    if( nvars > 0 ) {
      contess_assert(vardim);
      contess_assert(centering);
      contess_assert(varnames);
      contess_assert(vars);
    }
    Visit_writer().write_unstructured_mesh(
      filename, 
      useBinary, 
      npts, 
      pts, 
      ncells, 
      celltypes,
      conn, 
      nvars, 
      vardim, 
      centering, 
      varnames,
      vars
    );
  }
};