import bpy

from ..config import __addon_name__
from ..preference.AddonPreferences import ExampleAddonPreferences


class OBJECT_OT_add_sphere(bpy.types.Operator):
    """Add a UV Sphere"""
    bl_idname = "object.add_sphere"
    bl_label = "Add Sphere"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, location=(0, 0, 0))
        return {'FINISHED'}

