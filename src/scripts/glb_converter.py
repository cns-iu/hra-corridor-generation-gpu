# Single mesh converter. Please refer to multi-LOD for the multi-mesh converter
import argparse
import os

import numpy as np
import pygltflib
from pygltflib import GLTF2


def off_to_primitives(off_file_path):
    with open(off_file_path, 'r') as f:
        lines = f.read().splitlines()
        meta_data = lines[1].split(' ')
        n_points = int(meta_data[0])
        n_triangles = int(meta_data[1])

        i = 2
        # skip empty lines
        while i < len(lines) and not lines[i]:
            i += 1

        points = []
        triangles = []

        while n_points:
            point_3 = lines[i].split()
            p1, p2, p3 = float(point_3[0]), float(point_3[1]), float(point_3[2])
            points.append([p1, p2, p3])
            n_points -= 1
            i += 1

        while n_triangles:
            triangle_3 = lines[i].split()
            t1, t2, t3 = int(triangle_3[1]), int(triangle_3[2]), int(triangle_3[3])
            triangles.append([t1, t2, t3])
            n_triangles -= 1
            i += 1

        points = np.array(points, dtype=np.float32)
        triangles = np.array(triangles, dtype=np.uint32)

        return points, triangles


def converter(points, triangles):
    triangles_binary_blob = triangles.flatten().tobytes()
    points_binary_blob = points.tobytes()
    gltf = pygltflib.GLTF2(
        scene=0,
        scenes=[pygltflib.Scene(nodes=[0])],
        nodes=[pygltflib.Node(mesh=0)],
        meshes=[
            pygltflib.Mesh(
                primitives=[
                    pygltflib.Primitive(
                        attributes=pygltflib.Attributes(POSITION=1), indices=0
                    )
                ],
                name='corridor',
            )
        ],
        accessors=[
            pygltflib.Accessor(
                bufferView=0,
                componentType=pygltflib.UNSIGNED_INT,
                count=triangles.size,
                type=pygltflib.SCALAR,
                max=[int(triangles.max())],
                min=[int(triangles.min())],
                name='accessorIndices',
            ),
            pygltflib.Accessor(
                bufferView=1,
                componentType=pygltflib.FLOAT,
                count=len(points),
                type=pygltflib.VEC3,
                max=points.max(axis=0).tolist(),
                min=points.min(axis=0).tolist(),
                name='accessorPositions',
            ),
        ],
        bufferViews=[
            pygltflib.BufferView(
                buffer=0,
                byteLength=len(triangles_binary_blob),
                target=pygltflib.ELEMENT_ARRAY_BUFFER,
            ),
            pygltflib.BufferView(
                buffer=0,
                byteOffset=len(triangles_binary_blob),
                byteLength=len(points_binary_blob),
                target=pygltflib.ARRAY_BUFFER,
            ),
        ],
        buffers=[
            pygltflib.Buffer(
                byteLength=len(triangles_binary_blob) + len(points_binary_blob)
            )
        ],
    )
    gltf.set_binary_blob(triangles_binary_blob + points_binary_blob)

    return gltf

def generate_glb(input_file_off, output_file_glb):
    
    points, triangles = off_to_primitives(input_file_off)
    if len(points) == 0:
        continue
    gltf = converter(points, triangles)
    gltf.save(output_file_glb)

def generate_glbs(input_dir, output_dir):

    files = os.listdir(input_dir)
    for f in files:
        if f.endswith('.off'):
            # step 1
            off_file_path = os.path.join(input_dir, f)
            points, triangles = off_to_primitives(off_file_path)
            if len(points) == 0:
                continue
            # step 2
            gltf = converter(points, triangles)
            corridor_name = f[:-4]
            output_file_path = os.path.join(output_dir, corridor_name + '.glb')
            # output_file_path_2 = os.path.join(output_dir, corridor_name + '.gltf')

            gltf.save(output_file_path)
            # gltf.save(output_file_path_2)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file_off')
    parser.add_argument('output_file_glb')
    args = parser.parse_args()
    
    input_file_off = args.input_file_off
    output_file_glb = args.output_file_glb

    generate_glb(input_file_off, output_file_glb)

