// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

//
//   Copyright 2013 Pixar
//
//   Licensed under the Apache License, Version 2.0 (the "Apache License")
//   with the following modification; you may not use this file except in
//   compliance with the Apache License and the following modification to it:
//   Section 6. Trademarks. is deleted and replaced with:
//
//   6. Trademarks. This License does not grant permission to use the trade
//      names, trademarks, service marks, or product names of the Licensor
//      and its affiliates, except as required to comply with Section 4(c) of
//      the License and to reproduce the content of the NOTICE file.
//
//   You may obtain a copy of the Apache License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the Apache License with the above modification is
//   distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//   KIND, either express or implied. See the Apache License for the specific
//   language governing permissions and limitations under the Apache License.
//

#ifndef SHAPE_UTILS_H
#define SHAPE_UTILS_H

#include "common.h"
#include <map>
#include <optional>
#include <string>
#include <vector>

//------------------------------------------------------------------------------
enum Scheme { kBilinear = 0, kCatmark, kLoop };

//------------------------------------------------------------------------------

struct Shape {
  // full(er) spec here: http://paulbourke.net/dataformats/mtl/
  struct material {

    material();

    std::string name;

    real_t ka[3],  // ambient
        kd[3],     // diffuse
        ks[3],     // specular
        ns,        // specular exponent
        ni,        // optical density (1.0=no refraction, glass=1.5)
        sharpness, // reflection sharpness
        tf[3],     // transmission filter
        d;         // dissolve factor (1.0 = opaque)

    int illum;
  };

  struct tag {

    static std::optional<tag> parseTag(char const *stream);

    std::string genTag() const;

    std::string name;
    std::vector<int> intargs;
    std::vector<real_t> floatargs;
    std::vector<std::string> stringargs;
  };

  static Shape *parseObj(char const *Shapestr, Scheme schme,
                         bool isLeftHanded = false, int axis = 1,
                         bool parsemtl = false);

  void parseMtllib(char const *stream);

  std::string genShape(char const *name) const;

  std::string genObj() const;

  std::string genRIB() const;

  int GetNumVertices() const { return (int)verts.size() / 3; }

  int GetNumFaces() const { return (int)nvertsPerFace.size(); }

  bool HasUV() const { return !(uvs.empty() || faceuvs.empty()); }

  int GetFVarWidth() const { return HasUV() ? 2 : 0; }

  std::vector<real_t> verts;
  std::vector<real_t> uvs;
  std::vector<real_t> normals;
  std::vector<int> nvertsPerFace;
  std::vector<uint32_t> faceverts;
  std::vector<uint32_t> faceuvs;
  std::vector<uint32_t> facenormals;
  std::vector<tag> tags;
  Scheme scheme = kCatmark;
  bool isLeftHanded = false;

  char FindMaterial(char const *name) {
    for (int i = 0; i < (int)mtls.size(); ++i) {
      if (mtls[i].name == name) {
        return i;
      }
    }
    return -1;
  }

  std::string mtllib;
  std::vector<unsigned short> mtlbind;
  std::vector<material> mtls;
};

//------------------------------------------------------------------------------

#endif /* SHAPE_UTILS_H */
