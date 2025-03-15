import bpy

from ..config import __addon_name__
from ..operators.AddonOperators import castRays
from ....common.i18n.i18n import i18n
from ....common.types.framework import reg_order


class mainPanel(bpy.types.Panel):
    """Creates a Panel in the Object properties window"""
    bl_label = "Create Drawing"
    bl_idname = "VIEW3D_PT_add_sphere"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Grease Pencil Extension'

    def draw(self, context):
        layout = self.layout
        layout.operator("object.cast_rays", text="Cast Rays")
            
