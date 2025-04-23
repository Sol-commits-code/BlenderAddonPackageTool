import bpy
import mathutils
import math
import bmesh

from ..config import __addon_name__
from ..preference.AddonPreferences import ExampleAddonPreferences
from ..panels.AddonPanels import mainPanel

#Principle class that executes the test script
class test_script(bpy.types.Operator):
    """Executes Test Script"""
    bl_idname = "object.test_script"
    bl_label = "Execute Test Script"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        test_script.getverts()
        return {'FINISHED'}
        
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

            print(f"Selected face vertex coordinates (object: {obj.name}):")

            for face in bm.faces:
                if face.select:
                    print(f"\nFace index: {face.index}")
                    for vert in face.verts:
                        co_world = obj.matrix_world @ vert.co
                        print(f"Vertex: ({co_world.x:.4f}, {co_world.y:.4f}, {co_world.z:.4f})")
    
def register():
    # Register the castRays class
    bpy.utils.register_class(test_script)

def unregister():
    # Register the castRays class
    bpy.utils.register_class(test_script)

if __name__ == "__main__":
    register()

