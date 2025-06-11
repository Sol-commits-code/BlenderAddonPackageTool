import bpy
import mathutils
import math
import bmesh
import numpy as np
from collections import deque

from ..config import __addon_name__
from ..preference.AddonPreferences import ExampleAddonPreferences
from ..panels.AddonPanels import mainPanel
from typing import Tuple, List

#Principle class that executes the test script
class test_script(bpy.types.Operator):

    """Executes Test Script"""
    bl_idname = "object.test_script"
    bl_label = "Execute Test Script"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        test_script.main()
        return {'FINISHED'}
    
    def main():
        #main function calls all functions needed to perform drawing - currently used for testing

        # Stand in varibles - Temporary - Make input arguments eventually. 
        faceID = 0
        numIntersections = 1    

        #Perform Drawing for a given face 
        test_script.selectFace(faceID)
        quad = test_script.getverts()
        root = test_script.build_dyadic_binary_tree(max_nodes=5)
        interceptPoints = test_script.print_tree_breadth_first(root)
        test_script.drawLines(quad, interceptPoints)
       

    def selectFace(faceID):
        #Face select test space - delete when done:
        obj = bpy.context.active_object
        if obj and obj.type == 'MESH':
            test_script.select_face_by_index(obj, faceID)  # Change to the desired index
        else:
            print("No mesh object selected.")   

    
    # Print Vertex coordinates of selected face
    def getverts():
        
        obj = bpy.context.active_object
        if obj is None or obj.type != 'MESH':
            print("No mesh object selected.")
        else:
            # Switch to edit mode and get the BMesh
            bpy.ops.object.mode_set(mode='EDIT')
            bm = bmesh.from_edit_mesh(obj.data)

            # Update BMesh to reflect current selection
            bm.faces.ensure_lookup_table()
            vert_coords = []

            for face in bm.faces:
                if face.select:
                    for vert in face.verts:
                        co_world = obj.matrix_world @ vert.co
                        vert_coords.append(co_world) #creates a list of the x,y,z vectors of each vert of the face
                        
                        
            
            # Disect vertex coordinates into their respective vertices for easier data handeling 
            vert1 = vert_coords[0]
            vert2 = vert_coords[1]
            vert3 = vert_coords[2]
            vert4 = vert_coords[3]


            # Create a matrix containing the vert coordinates of the active quad  
            quad = [
                [vert1[0], vert1[1], vert1[2]],
                [vert2[0], vert2[1], vert2[2]],
                [vert3[0], vert3[1], vert3[2]],
                [vert4[0], vert4[1], vert4[2]]
            ]
            return quad     
        
    def drawLines(quad, intersections):

        for intersectionPoints in intersections:
            # Given a position along the edge defined as a percentage between 0 & 1 -> Return the x,y,z coordinates of this position
            line1Point1 = test_script.point_along_line(quad[0], quad[3], intersectionPoints) 
            line1Point2 = test_script.point_along_line(quad[1], quad[2], intersectionPoints)
            test_script.create_edge_between_points(line1Point1, line1Point2)
            


    ## LINE INTERSECT FUNCTIONS     
    def point_along_line(p1, p2, t):
        """
        Compute a point along the line from p1 to p2 at normalized position t.

        Parameters:
            p1 (tuple of float): (x1, y1, z1), the start point (t = 0).
            p2 (tuple of float): (x2, y2, z2), the end point   (t = 1).
            t (float): normalized position along the line, 0 <= t <= 1.

        Returns:
            tuple of float: (x, y, z) coordinate of the interpolated point.

        Raises:
            ValueError: if t is not in the [0.0, 1.0] range.
        """
        x1, y1, z1 = p1
        x2, y2, z2 = p2

        # Linear interpolation: P = p1 + t * (p2 - p1)
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        z = z1 + t * (z2 - z1)

        return (x, y, z)
    
    def addsphere(coords):  
        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.02, location=(coords[0], coords[1], coords[2]))
        return {'FINISHED'}
    
    def select_face_by_index(obj, face_index):
        # Ensure we're in object mode before switching to edit
        if bpy.context.object.mode != 'EDIT':
            bpy.ops.object.mode_set(mode='EDIT')

        # Access the mesh via BMesh
        mesh = bmesh.from_edit_mesh(obj.data)

        # Deselect all faces
        for f in mesh.faces:
            f.select = False

        # Select the target face if index is valid
        if 0 <= face_index < len(mesh.faces):
            mesh.faces.ensure_lookup_table()
            mesh.faces[face_index].select = True
            print(f"Face {face_index} selected.")
        else:
            print(f"Face index {face_index} out of range.")

        # Update mesh to show selection
        bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
   
    def create_edge_between_points(p1, p2, object_name="EdgeObject"):
        """
        Creates a new mesh object with two vertices connected by an edge.

        Args:
            p1, p2: Tuples of (x, y, z) coordinates.
            object_name: Name of the new object created.
        """
        # Create a new mesh and object
        mesh = bpy.data.meshes.new(name=f"{object_name}_Mesh")
        obj = bpy.data.objects.new(object_name, mesh)
        bpy.context.collection.objects.link(obj)

        # Create a BMesh and add the geometry
        bm = bmesh.new()

        v1 = bm.verts.new(p1)
        v2 = bm.verts.new(p2)
        bm.edges.new((v1, v2))

        bm.to_mesh(mesh)
        bm.free()
    
    class Node:
        def __init__(self, value: float, depth: int):
            self.value = value
            self.depth = depth
            self.left = None
            self.right = None

        def __repr__(self):
            return f"Node(value={self.value:.4f}, depth={self.depth})"

    def build_dyadic_binary_tree(max_nodes: int) -> Node:
        """
        Builds a complete binary tree where:
        - Each child is offset from its parent by (1 / 2^depth)
        - Root starts at 0.5 (center of [0,1])
        - All node values stay within [0, 1]

        Args:
            max_nodes (int): Total number of nodes to create.

        Returns:
            Node: Root of the tree.
        """
        if max_nodes <= 0:
            return None

        root_value = 0.5
        root = test_script.Node(value=root_value, depth=1)
        queue = deque([root])
        node_count = 1

        while queue and node_count < max_nodes:
            parent = queue.popleft()
            next_depth = parent.depth + 1
            offset = 1 / (2 ** next_depth)

            # Create left child
            if node_count < max_nodes:
                left_value = parent.value - offset
                left = test_script.Node(value=left_value, depth=next_depth)
                parent.left = left
                queue.append(left)
                node_count += 1

            # Create right child
            if node_count < max_nodes:
                right_value = parent.value + offset
                right = test_script.Node(value=right_value, depth=next_depth)
                parent.right = right
                queue.append(right)
                node_count += 1

        return root

    def print_tree_breadth_first(root):
        """Print node values in level-order."""
        if not root:
            print("Empty tree.")
            return
        
        interceptPoints = []

        queue = deque([root])
        while queue:
            node = queue.popleft()
            interceptPoints.append(node.value)
            if node.left:
                queue.append(node.left)
            if node.right:
                queue.append(node.right)
        return interceptPoints


def register():
    # Register the castRays class
    bpy.utils.register_class(test_script)

def unregister():
    # Register the castRays class
    bpy.utils.register_class(test_script)

if __name__ == "__main__":
    register()

