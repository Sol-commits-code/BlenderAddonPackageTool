import bpy

from ..config import __addon_name__
from ..operators.AddonOperators import OBJECT_OT_add_sphere
from ....common.i18n.i18n import i18n
from ....common.types.framework import reg_order


class VIEW3D_PT_add_sphere_panel(bpy.types.Panel):
    """Creates a Panel in the Object properties window"""
    bl_label = "Add Sphere Panel"
    bl_idname = "VIEW3D_PT_add_sphere"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Add Sphere'

    def draw(self, context):
        layout = self.layout
        layout.operator("object.add_sphere", text="Add Sphere")

