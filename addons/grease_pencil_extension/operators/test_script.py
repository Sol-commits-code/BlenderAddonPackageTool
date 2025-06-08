import bpy
import mathutils
import math
import bmesh
import numpy as np

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
        test_script.selectFace()
        test_script.getverts()
        return {'FINISHED'}
    
    def selectFace():
        #Face select test space - delete when done:
        obj = bpy.context.active_object
        if obj and obj.type == 'MESH':
            test_script.select_face_by_index(obj, face_index=2  )  # Change to the desired index
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

            print(f"Selected face vertex coordinates (object: {obj.name}):")

            for face in bm.faces:
                if face.select:
                    print(f"\nFace index: {face.index}")
                    for vert in face.verts:
                        co_world = obj.matrix_world @ vert.co
                        vert_coords.append(co_world) #creates a list of the x,y,z vectors of each vert of the face
                        print(f"Vertex: ({co_world.x:.4f}, {co_world.y:.4f}, {co_world.z:.4f})")
                        
            
            # Disect vertex coordinates into their respective vertices for easier data handeling 
            vert1 = vert_coords[0]
            vert2 = vert_coords[1]
            vert3 = vert_coords[2]
            vert4 = vert_coords[3]


            # Create a matrix containing the vert coordinates of the active quad  
            quad = [
                (vert1[0], vert1[1], vert1[2]),
                (vert2[0], vert2[1], vert2[2]),
                (vert3[0], vert3[1], vert3[2]),
                (vert4[0], vert4[1], vert4[2])
            ]

            positionAlongEdge = 0.5 

            pointLine1 = test_script.point_along_line((vert1[0], vert1[1], vert1[2]), (vert4[0], vert4[1], vert4[2]), positionAlongEdge) 
            pointLine2 = test_script.point_along_line((vert2[0], vert2[1], vert2[2]), (vert3[0], vert3[1], vert3[2]), positionAlongEdge)
            #pointLine3 = test_script.point_along_line((vert1[0], vert1[1], vert1[2]), (vert2[0], vert2[1], vert2[2]), positionAlongEdge) 
            #pointLine4 = test_script.point_along_line((vert4[0], vert4[1], vert4[2]), (vert3[0], vert3[1], vert3[2]), positionAlongEdge)
       
            test_script.addsphere(pointLine1)
            test_script.addsphere(pointLine2)
            #test_script.addsphere(pointLine3)
            #test_script.addsphere(pointLine4)  
            


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
        if not (0.0 <= t <= 1.0):
            raise ValueError(f"t must be between 0 and 1 inclusive, got {t!r}")

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


def register():
    # Register the castRays class
    bpy.utils.register_class(test_script)

def unregister():
    # Register the castRays class
    bpy.utils.register_class(test_script)

if __name__ == "__main__":
    register()

