# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

#!/usr/bin/env python3

# Example call:
# python3 src/py/process_svg.py torus_quad_v2_qi_chained.svg --out highlight_vertices.svg
# python3 process_svg.py data/torus_quad_v2_qi_chained.svg --out highlight_vertices.svg

import sys
import os
import re
import subprocess
import operator
import argparse
from svgpathtools import (Line, Path, Arc, svg2paths,
                          wsvg, disvg, misctools)
import xml.etree.ElementTree as ET
import random
from dataclasses import dataclass
from enum import Enum
from typing import List, Set

# Topological simplification parameters

arcLength_threshold = 10
arcLength_pre_threshold = 2

# SVG Rendering style options

random_color_strokes = True
show_points = True
stroke_width = 2

#################### Data structures/types #############


class VertexType(Enum):
    REGULAR = 'c'
    ENDPOINT = 'e'
    JUNCTION = 'j'


VertexID = int
CurveID = int


@dataclass
class Vertex:
    ID: VertexID
    pos: complex
    vtype: VertexType
    z_value: float
    refs: Set[CurveID]   # paths that connect to this vertex


@dataclass
class Curve:
    ID: CurveID
    vs: List[VertexID]


#################### PARSING AN SVG ####################

def get_svg_size(svg_file):
    viewbox = None
    tree = ET.parse(svg_file)
    root = tree.getroot()
    if "viewBox" in root.attrib:
        if ',' in root.attrib['viewBox']:
            viewbox = [float(v)
                       for v in root.attrib['viewBox'].split(',')[0:4]]
            width = viewbox[2] - viewbox[0]
            height = viewbox[3] - viewbox[1]
        else:
            viewbox = [float(v)
                       for v in root.attrib['viewBox'].split(' ')[0:4]]
            width = viewbox[2] - viewbox[0]
            height = viewbox[3] - viewbox[1]
    elif "width" in root.attrib and "height" in root.attrib:
        width = float(root.attrib['width'].strip('px'))
        height = float(root.attrib['height'].strip('px'))
    else:
        print("Error:\tparsing svg failed")
    if not viewbox:
        viewbox = [0, 0, width, height]
    viewbox = tuple([int(coord) for coord in viewbox])
    return viewbox


def build_indices(paths, attributes):
    # convert input datastructures into curve network datastructures

    vertex_index = {}
    path_index = {}

    for i in range(len(paths)):
        if not isinstance(paths[i], Path):
            continue

        vs = []
        z_list = None
        if 'z' in attributes[i]:
            z_list = [float(z) for z in attributes[i]['z'].split(' ')]

        for j in range(len(paths[i])+1):

            if j < len(paths[i]):  # all but the last point on the curve
                pos = paths[i][j][0]
            else:
                pos = paths[i][j-1][1]  # last point on the curve

            type_str = attributes[i]['type_v{}'.format(j)]

            vertID = int(attributes[i][f'index{j}'])

            correspID = None
            if type_str[0] == 'j':
                correspID = int(type_str[2:])

            # check if this vertex already has an ID, and update index
            if correspID not in vertex_index:
                z_value = -1
                if z_list:
                    z_value = z_list[j]
                vertex_index[vertID] = Vertex(
                    vertID, pos, type_str[0], z_value, {i})
            else:
                vertex_index[correspID].refs.add(i)
                vertID = correspID

            vs.append(vertID)

        path_index[i] = Curve(i, vs)

    return vertex_index, path_index

# for debugging


def check_datastructures(vertex_index, path_index):
    # check that each vertex has the correct list of refs
    vind = {v: set() for v in vertex_index}
    for p in path_index:
        for v in path_index[p].vs:
            try:
                vind[v].add(p)
            except Exception as e:
                print('Warning: Lookup failed.')

    for v in vertex_index:
        # if vind[v] != vertex_index[v].refs:
        #     print(vind[v])
        #     print(vertex_index[v].refs)
        # assert vind[v] == vertex_index[v].refs
        pass

    return True


def parse_svg(filename):
    viewbox = get_svg_size(filename)
    p, attr = svg2paths(filename)
    attr[0] = {key: v for key, v in attr[0].items(
    ) if 'type_v' in key or 'index' in key or 'z' in key}

    return p, attr, viewbox

#################### SIMPLIFICATION #############################


def arcLength(path, vertex_index, start=0, end=None):
    vs = path.vs

    if end == None:
        end = len(vs)

    length = 0
    for i in range(start, end-1):
        v0 = vertex_index[vs[i]].pos
        v1 = vertex_index[vs[i+1]].pos

        length = length + Line(v0, v1).length()

    return length


def types(path, vertex_index):  # return a list with all the types
    try:
        return list(map(lambda x: vertex_index[x].vtype, path.vs))
    except Exception as e:
        print('Warning: Lookup failed.')
        return None


def delete_path(pID, vertex_index, path_index):

    assert check_datastructures(vertex_index, path_index)

    path = path_index[pID]

    # update all vertices
    for vnum in range(len(path.vs)):
        try:
            v = vertex_index[path.vs[vnum]]
        except Exception as e:
            print('Warning: Lookup failed.')
            continue
        if pID in v.refs:
            v.refs.remove(pID)
        if len(v.refs) == 0:
            del vertex_index[path.vs[vnum]]
        else:
            # this vertex must have been a junction
            v.vtype == 'c'
    del path_index[pID]

    assert check_datastructures(vertex_index, path_index)


def nextJunctionFromStart(path, vtypes):
    N = len(vtypes)
    for i in range(N-1):
        if vtypes[i] == 'j':
            return i

    return None


def nextJunctionFromEnd(path, vtypes):
    N = len(vtypes)
    for i in range(N-1, 0, -1):
        if vtypes[i] == 'j':
            return i

    return None


def delete_path_to_vertex(pID, vertex_index, path_index, jind, delete_to_vertex):
    # delete path up to a junction vertex, and merge the paths at that point
    # delete_to_vertex: True if deleting 0...jind, False if deleting jind+1...0

    path = path_index[pID]   # this path
    vID = path.vs[jind]      # ID of the junction verterx
    vj = vertex_index[vID]   # junction vertex

    assert vj.vtype == 'j'
    assert check_datastructures(vertex_index, path_index)
#    print(path.vs)

    if delete_to_vertex:
        drange = range(jind)
    else:
        drange = range(jind+1, len(path.vs))

    # update all  vertices before the junction vertex
    for vnum in drange:
        try:
            v = vertex_index[path.vs[vnum]]
            v.refs.remove(pID)
            if len(v.refs) == 0:
                del vertex_index[path.vs[vnum]]
            else:
                # this vertex must have been a junction
                v.vtype = 'c'
        except Exception as e:
            print('Warning: Lookup failed.')

    # remove vertices from path
    if delete_to_vertex:
        path.vs = path.vs[jind:]
    else:
        path.vs = path.vs[:jind+1]

    assert check_datastructures(vertex_index, path_index)

    if len(vj.refs) == 1:
        # the path connected to itself, and now it's now a loop

        vj.vtype = 'c'
        vertex_index[path.vs[-1]].vtype = 'c'
#        path.vs = path.vs[1:]

#        print(vID)
#        print(vj.refs)
#        print(path.vs)

        assert check_datastructures(vertex_index, path_index)
    else:
        # the vertex connected two paths; merge them
        pIDj = (vj.refs-{pID}).pop()   # the ID of the other path
        pathj = path_index[pIDj]       # the other path

        # check if this was the start or end point of the connecting path
        if vID == pathj.vs[0]:
            vs_j = pathj.vs[::-1]
        else:
            vs_j = pathj.vs

        if delete_to_vertex:
            # this path ends with the same vertex the other begins with
            assert vs_j[-1] == path.vs[0]
            # merge the vertex lists (without duplicating the junction vertex)
            path.vs = vs_j[:-1] + path.vs
        else:
            assert vs_j[-1] == path.vs[-1]
            # merge the vertex lists (without duplicating the junction vertex)
            path.vs = vs_j[:-1] + path.vs[::-1]

        # update vertices
        for v in set(vs_j):
            if pIDj in vertex_index[v].refs:
                vertex_index[v].refs.remove(pIDj)
            vertex_index[v].refs.add(pID)

        # delete the old path
        del path_index[pIDj]

        assert check_datastructures(vertex_index, path_index)


def simplify(vertex_index, path_index):

    #    vertex_index = input_vertex_index
    #    path_index = input_path_index

    changed = True
    while changed:
        changed = False

        pids = list(path_index)
        # iterate over all subchains, see which ones can be deleted
        for pID in pids:
            if pID not in path_index:  # the stroke could have been deleted by another iteration
                continue

            path = path_index[pID]

            if not types(path, vertex_index):
                return

            # check if this is a simple path with nothing happening between the endpoints
            # (this includes cases (b,c,d) in Benard et al [2014])
            vtypes = types(path, vertex_index)
            if arcLength(path, vertex_index) < arcLength_pre_threshold and 'e' not in vtypes[1:-1] and 'j' not in vtypes[1:-1]:
                delete_path(pID, vertex_index, path_index)
                changed = True
            else:
                # check for short curves that terminate in junctions (case (a) in BHK2014)
                vtypes = types(path, vertex_index)
                if vtypes[0] == 'e':
                    jnum = nextJunctionFromStart(path, vtypes)

                    if jnum != None and arcLength(path, vertex_index, 0, jnum) < arcLength_pre_threshold:
                        delete_path_to_vertex(
                            pID, vertex_index, path_index, jnum, True)
                        changed = True

                vtypes = types(path, vertex_index)
                N = len(vtypes)
                if vtypes[N-1] == 'e':
                    jnum = nextJunctionFromEnd(path, vtypes)

                    if jnum != None and arcLength(path, vertex_index, jnum, N) < arcLength_pre_threshold:
                        delete_path_to_vertex(
                            pID, vertex_index, path_index, jnum, False)
                        changed = True

    changed = True
    while changed:
        changed = False

        pids = list(path_index)
        # iterate over all subchains, see which ones can be deleted
        for pID in pids:
            if pID not in path_index:  # the stroke could have been deleted by another iteration
                continue

            path = path_index[pID]

            if not types(path, vertex_index):
                return

            # check if this is a simple path with nothing happening between the endpoints
            # (this includes cases (b,c,d) in Benard et al [2014])
            vtypes = types(path, vertex_index)
            if arcLength(path, vertex_index) < arcLength_threshold and 'e' not in vtypes[1:-1] and 'j' not in vtypes[1:-1]:
                delete_path(pID, vertex_index, path_index)
                changed = True
            else:
                # check for short curves that terminate in junctions (case (a) in BHK2014)
                vtypes = types(path, vertex_index)
                if vtypes[0] == 'e':
                    jnum = nextJunctionFromStart(path, vtypes)

                    if jnum != None and arcLength(path, vertex_index, 0, jnum) < arcLength_threshold:
                        delete_path_to_vertex(
                            pID, vertex_index, path_index, jnum, True)
                        changed = True

                vtypes = types(path, vertex_index)
                N = len(vtypes)
                if vtypes[N-1] == 'e':
                    jnum = nextJunctionFromEnd(path, vtypes)

                    if jnum != None and arcLength(path, vertex_index, jnum, N) < arcLength_threshold:
                        delete_path_to_vertex(
                            pID, vertex_index, path_index, jnum, False)
                        changed = True

#################### SAVE A NEW SVG TO DISK ####################


def generate_svg(path_index, vertex_index, viewbox, out_filename, scale=1):
    out_paths = []
    out_attributes = []
    random.seed(0)

    found_ends = []
    found_junctions = []

    for pnum in path_index:
        if random_color_strokes:
            color = [random.randint(128, 255), random.randint(
                0, 128), random.randint(0, 128)]
            random.shuffle(color)
        else:
            color = [0, 0, 0]

        style_str = 'fill:none;stroke-width:%d;stroke-opacity:1;stroke:#%02X%02X%02X' % (
            stroke_width, color[0], color[1], color[2])

        path = path_index[pnum]

        lines = []

        out_path = Path()
        for i in range(len(path.vs)-1):
            try:
                v0 = vertex_index[path.vs[i]].pos * scale
                v1 = vertex_index[path.vs[i+1]].pos * scale
                out_path.append(Line(v0, v1))
            except Exception as e:
                print('Warning: Lookup failed.')

        z_values = []
        for i in range(len(path.vs)):
            try:
                if vertex_index[path.vs[i]].z_value >= 0:
                    z_values.append(str(vertex_index[path.vs[i]].z_value))
            except Exception as e:
                print('Warning: Lookup failed.')

        if show_points:
            try:
                v0 = vertex_index[path.vs[0]]
                if v0.vtype == 'e':
                    found_ends.append(v0.pos)
                elif v0.vtype == 'j':
                    found_junctions.append(v0.pos)

                vN = vertex_index[path.vs[-1]]
                if vN.vtype == 'e':
                    found_ends.append(vN.pos)
                elif vN.vtype == 'j':
                    found_junctions.append(vN.pos)
            except Exception as e:
                print('Warning: Lookup failed.')

        out_paths.append(out_path)
        attr = {'style': style_str}
        if z_values:
            attr['z'] = ' '.join(z_values)
        out_attributes.append(attr)

    node_colors = ['orange']*len(found_ends) + ['green']*len(found_junctions)

    viewbox = [int(v * scale) for v in viewbox]
    print(viewbox)

    if random_color_strokes and show_points:
        wsvg(out_paths, node_radii=[5]*(len(found_ends) + len(found_junctions)),
             nodes=found_ends+found_junctions, node_colors=node_colors,
             attributes=out_attributes, dimensions=(
             int(viewbox[2]-viewbox[0]), int(viewbox[3]-viewbox[1])),
             viewbox=f'{viewbox[0]} {viewbox[1]} {viewbox[2]} {viewbox[3]}', filename=out_filename)
    else:
        with open(out_filename, 'w') as f:
            header = f'<?xml version="1.0" ?>\n<svg viewBox="{viewbox[0]} {viewbox[1]} {viewbox[2]} {viewbox[3]}" version="1.1" xmlns="http://www.w3.org/2000/svg">\n'
            paths = ''
            for i in range(len(out_paths)):
                p = f'<g fill="none" stroke="rgb(0,0,0)" stroke-width="{stroke_width}">'
                p += '<path d="'
                for j in range(len(out_paths[i])):
                    if j == 0:
                        p += f'M{out_paths[i][j][0].real} {out_paths[i][j][0].imag}'
                    elif j + 1 == len(out_paths[i]):
                        p += f' L{out_paths[i][j][0].real} {out_paths[i][j][0].imag}'
                        p += f' L{out_paths[i][j][1].real} {out_paths[i][j][1].imag}'
                    else:
                        p += f' L{out_paths[i][j][0].real} {out_paths[i][j][0].imag}'
                p += '"'
                z_str = out_attributes[i]['z']
                p += f' z="{z_str}" />'
                p += '</g>\n'
                paths += p
            f.write(header + paths + '</svg>\n')


def generate_svg_scale(paths, viewbox, out_filename, scale=1):
    out_paths = []
    out_attributes = []
    random.seed(0)

    found_ends = []
    found_junctions = []

    for p in paths:
        if random_color_strokes:
            color = [random.randint(128, 255), random.randint(
                0, 128), random.randint(0, 128)]
            random.shuffle(color)
        else:
            color = [0, 0, 0]

        style_str = 'fill:none;stroke-width:%d;stroke-opacity:1;stroke:#%02X%02X%02X' % (
            stroke_width, color[0], color[1], color[2])

        lines = []

        out_path = Path()
        for i in range(len(p)):
            v0 = p[i].start * scale
            v1 = p[i].end * scale
            out_path.append(Line(v0, v1))

        out_paths.append(out_path)
        attr = {'style': style_str}
        out_attributes.append(attr)

    node_colors = ['orange']*len(found_ends) + ['green']*len(found_junctions)

    viewbox = [int(v * scale) for v in viewbox]
    print(viewbox)

    with open(out_filename, 'w') as f:
        header = f'<?xml version="1.0" ?>\n<svg viewBox="{viewbox[0]} {viewbox[1]} {viewbox[2]} {viewbox[3]}" version="1.1" xmlns="http://www.w3.org/2000/svg">\n'
        paths = ''
        for i in range(len(out_paths)):
            p = f'<g fill="none" stroke="rgb(0,0,0)" stroke-width="{stroke_width}">'
            p += '<path d="'
            for j in range(len(out_paths[i])):
                if j == 0:
                    p += f'M{out_paths[i][j][0].real} {out_paths[i][j][0].imag}'
                elif j + 1 == len(out_paths[i]):
                    p += f' L{out_paths[i][j][0].real} {out_paths[i][j][0].imag}'
                    p += f' L{out_paths[i][j][1].real} {out_paths[i][j][1].imag}'
                else:
                    p += f' L{out_paths[i][j][0].real} {out_paths[i][j][0].imag}'
            p += '" />'
            p += '</g>\n'
            paths += p
        f.write(header + paths + '</svg>\n')

#################### COMMAND-LINE PROCESSING ####################


def parse_args():
    global random_color_strokes, show_points

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('svg', type=str, nargs='+',
                        help='Input svg file.')
    parser.add_argument('--out', help='Output svg file.')
    parser.add_argument(
        '--plain', help='Whether to render the result with a plain style.', action='store_true')
    parser.add_argument(
        '--scale', help='A scale factor for the output.', type=float)
    parser.add_argument(
        '--skip', help='Skip simplification and only scale the input.', action='store_true')

    args = parser.parse_args()

    if args.plain:
        random_color_strokes = False
        show_points = False

    return args.svg, args.out, args.scale, args.skip


if __name__ == '__main__':
    # parse input arguments
    path_filenames, out_filename, scale, skip = parse_args()
    if not scale:
        scale = 1

    # load SVG from file
    paths, attributes, viewbox = parse_svg(
        path_filenames[0])  # only process the first file

    if not skip:
        # process the path data into a curve network data structure
        vertex_index, path_index = build_indices(paths, attributes)

        # apply topological simplification to the data structure
        simplify(vertex_index, path_index)

        # output the curve network as an SVG file
        generate_svg(path_index, vertex_index, viewbox, out_filename, scale)
    else:
        # Simply scale the input
        generate_svg_scale(paths, viewbox, out_filename, scale)
