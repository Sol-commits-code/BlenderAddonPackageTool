import bpy

from ..config import __addon_name__
#from ..operators.AddonOperators import castRays
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

        layout.label(text="Ray Count:")
        layout.prop(context.scene, "rays_x", text="Number of Rays: X")
        layout.prop(context.scene, "rays_y", text="Number of Rays: Y")
        

        layout.label(text="Execute Drawing:")
        layout.operator("object.cast_rays", text="Cast Rays")
        layout.operator("object.test_script", text="Execute Test Script")
            
    def register():
      
        bpy.types.Scene.rays_x = bpy.props.IntProperty(
            name="Number of Rays: X", default=10, min=1, description="Number of cast rays in x direction relative to camera"
        )
        bpy.types.Scene.rays_y = bpy.props.IntProperty(
            name="Number of Rays: Y", default=10, min=1, description="Number of cast rays in y direction relative to camera"
        )

    def unregister():
    
        del bpy.types.Scene.rays_x
        del bpy.types.Scene.rays_y


