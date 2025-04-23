"""
Filename: AddonOperators.py
Description: This script will perform a ray casting in the blender workspace. The rays will emit 
From the princple camera in the scene. The number of rays is set by the ray_x and ray_y user inputs 
as defined in the mainPanel class in AddonPanels. It its current state, the function spawns a red 
shpere at the point of each intersection  

Author: Sol Sinclair
Date Created: 2025-03-15
Last Modified: 2025-03-15
Version: 1.0

License: [Specify license, e.g., MIT, GPL, or 'All rights reserved']
Dependencies:
    - Blender 3.x (if used in Blender)
    - NumPy (if applicable)
    - Any other necessary libraries

Notes:
    - TODO return a list containing intersection coordinates for each object 
"""

import bpy
import mathutils
import math
import bmesh

from ..config import __addon_name__
from ..preference.AddonPreferences import ExampleAddonPreferences
from ..panels.AddonPanels import mainPanel


#Principle class that perform the ray casting function 
class castRays(bpy.types.Operator):
    """Cast Rays to get Intersection Point"""
    bl_idname = "object.cast_rays"
    bl_label = "Cast Rays"
    bl_options = {'REGISTER', 'UNDO'}

    # following def executes the casting operation
    def execute(self, context):
        scene = bpy.context.scene
        camera = bpy.context.scene.camera
        
        # Check if scene has an active camera if true, perform casting
        if camera:
            castRays.cast_ray_from_camera(scene, camera)
        else:
            print("No active camera found!")
        return {'FINISHED'}

    # Adds spawned shepres into a blender collection for easy deletion from blender environment
    def get_or_create_collection(name="Raycast_Dots"):
        """Gets or creates a collection to store raycast dots."""
        if name in bpy.data.collections:
            return bpy.data.collections[name]
        else:
            new_collection = bpy.data.collections.new(name)
            bpy.context.scene.collection.children.link(new_collection)
            return new_collection

    # Meat and Potatos function 
    def cast_ray_from_camera(scene, camera):
        """Casts rays from the camera within its FOV and marks intersection points with a small sphere."""

        # Get the dependency graph
        depsgraph = bpy.context.evaluated_depsgraph_get()

        # Get camera transformation matrix
        cam_matrix = camera.matrix_world

        # Set ray origin at the camera's location
        ray_origin = cam_matrix.translation

        # Get camera FOV (horizontal and vertical)
        sensor_width = camera.data.sensor_width
        sensor_height = camera.data.sensor_height
        aspect_ratio = sensor_width / sensor_height
        fov_x = camera.data.angle  # Horizontal FOV in radians
        fov_y = fov_x / aspect_ratio  # Vertical FOV

        # Get or create the collection for storing dots
        dot_collection = castRays.get_or_create_collection()

        # Get user input for number of cast rays as defined in AddonPanels
        num_rays_x = bpy.context.scene.rays_x
        num_rays_y = bpy.context.scene.rays_y

        for i in range(num_rays_x):
            for j in range(num_rays_y):
                # Convert pixel space (normalized) to camera view space
                u = (i / (num_rays_x - 1)) * 2 - 1  # Ranges from -1 to 1 (left to right)
                v = (j / (num_rays_y - 1)) * 2 - 1  # Ranges from -1 to 1 (bottom to top)
                
                # Map to FOV
                angle_x = u * (fov_x / 2)  # Scale to horizontal FOV
                angle_y = v * (fov_y / 2)  # Scale to vertical FOV
                
                # Convert angles to world direction
                view_vector = mathutils.Vector((
                    math.tan(angle_x),
                    math.tan(angle_y),
                    -1  # Looking forward in the -Z direction
                )).normalized()

                # Transform view vector to world space
                view_vector = cam_matrix.to_3x3() @ view_vector
                
                # Perform ray casting using the dependency graph
                hit, location, normal, index, object, matrix = scene.ray_cast(depsgraph, ray_origin, view_vector)

                if hit:
                    print(f"Ray ({i}, {j}) hit {object.name} at {location}")
                    castRays.draw_dot_at_point(location, dot_collection)

    # Spawns the Shperes 
    def draw_dot_at_point(location, collection, radius=0.1):
        """Creates a small sphere at the given location and adds it to the specified collection."""
        
        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
        dot = bpy.context.object
        dot.name = "Raycast_Dot"
        dot.data.materials.append(castRays.get_or_create_material("RedMaterial"))

        # Add the dot to the specified collection
        collection.objects.link(dot)

        # Remove it from the default scene collection
        if dot.name in bpy.context.scene.collection.objects:
            bpy.context.scene.collection.objects.unlink(dot)

    # Colors the shperes red 
    def get_or_create_material(name):
        """Creates or retrieves a material with the given name."""
        mat = bpy.data.materials.get(name)
        if mat is None:
            mat = bpy.data.materials.new(name=name)
            mat.diffuse_color = (1, 0, 0, 1)  # Red color
        return mat
    
def register():
    # Register the castRays class
    bpy.utils.register_class(castRays)

def unregister():
    # Register the castRays class
    bpy.utils.register_class(castRays)

if __name__ == "__main__":
    register()

