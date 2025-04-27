import bpy
import mathutils
import math
import bmesh
import numpy as np

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
            vert_coords = []

            print(f"Selected face vertex coordinates (object: {obj.name}):")

            for face in bm.faces:
                if face.select:
                    print(f"\nFace index: {face.index}")
                    for vert in face.verts:
                        co_world = obj.matrix_world @ vert.co
                        vert_coords.append(co_world) #creates a list of the x,y,z vectors of each vert of the face
                        print(f"Vertex: ({co_world.x:.4f}, {co_world.y:.4f}, {co_world.z:.4f})")
                        
            
            print("now testing")
            #edge1 = [vert_coords[0], vert_coords[1]]
            edge2 = [vert_coords[2], vert_coords[3]]
            edge3 = [vert_coords[0], vert_coords[3]]
            edge4 = [vert_coords[1], vert_coords[2]]
            #print(edge1)

            vert1 = vert_coords[0]
            vert2 = vert_coords[1]
            vert3 = vert_coords[2]
            vert4 = vert_coords[3]

            # --- Example Usage ---

            quad = [
                (vert1[0], vert1[1], vert1[2]),
                (vert2[0], vert2[1], vert2[2]),
                (vert3[0], vert3[1], vert3[2]),
                (vert4[0], vert4[1], vert4[2])
            ]

            P = (0.5, 0.2, 2)  # interior point
            m1, m2 = 0, 100
            hits = test_script.intersect_lines_with_quad(quad, P, m1, m2)
            print(hits)
            line1 = (hits['line1'])
            line2 = (hits['line2'])
            line1_Vert1 = line1[0]
            line1_Vert2 = line1[1]
            line2_Vert1 = line2[0]
            line2_Vert2 = line2[1]
        
            test_script.addsphere(0.5, 0.2, 2)
            test_script.addsphere(line1_Vert1[0], line1_Vert1[1], line1_Vert1[2])
            test_script.addsphere(line1_Vert2[0], line1_Vert2[1], line1_Vert2[2])
            test_script.addsphere(line2_Vert1[0], line2_Vert1[1], line2_Vert1[2])
            test_script.addsphere(line2_Vert2[0], line2_Vert2[1], line2_Vert2[2])

            #test_script.addsphere(line1[0], line1[1], line1[2])
            #test_script.addsphere(line1[0], line1[1], line1[2])

            #test_script.addsphere(line2[0], line2[1], line2[2])
            #test_script.addsphere(line2[0], line2[1], line2[2])

        


            #test_script.addsphere(0.8, 0.4, 0.0)

    
    ## LINE INTERSECT FUNCTIONS     
    

    def intersect_lines_with_quad(vertices, point, m1, m2):
        """
        vertices : list of four (x,y,z) tuples, ordered around the convex quad
        point    : (x,y,z) tuple lying on the same plane
        m1, m2   : floats, the slopes of two lines in the quad's local frame
        
        Returns a dict:
        {
            'line1': [P1_back, P1_fwd],
            'line2': [P2_back, P2_fwd],
        }
        where each Pi is an (x,y,z) tuple of intersection points.
        """

        # --- basic 3D vector ops ---
        def sub(a,b):   return (a[0]-b[0],  a[1]-b[1],  a[2]-b[2])
        def add(a,b):   return (a[0]+b[0],  a[1]+b[1],  a[2]+b[2])
        def dot(a,b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        def cross(a,b): 
            return (a[1]*b[2] - a[2]*b[1],
                    a[2]*b[0] - a[0]*b[2],
                    a[0]*b[1] - a[1]*b[0])
        def mul(a,s):   return (a[0]*s, a[1]*s, a[2]*s)
        def norm(a):
            mag = math.sqrt(dot(a,a))
            return (a[0]/mag, a[1]/mag, a[2]/mag)
        
        # 1) Build an orthonormal (u,v) basis in the quad's plane:
        v1, v2, v3, _ = vertices
        u = norm(sub(v2, v1))
        n = cross(sub(v2, v1), sub(v3, v1))
        v = norm(cross(n, u))
        O = point  # origin for local coords

        # 2) Project quad verts into 2D local frame (relative to O):
        verts2d = []
        for V in vertices:
            d = sub(V, O)
            verts2d.append((dot(d, u), dot(d, v)))

        # 3) Build edge-list in that 2D frame:
        edges = []
        for i in range(4):
            a = verts2d[i]
            b = verts2d[(i+1) % 4]
            edges.append((a, b))
        
        # 2D cross-product (scalar)
        def cross2(a, b):
            return a[0]*b[1] - a[1]*b[0]

        # 4) Intersection solver for a given slope m:
        def find_hits(m):
            d = (1.0, m)   # direction of our line in 2D
            hits = []
            for a, b in edges:
                e = (b[0] - a[0], b[1] - a[1])
                det = cross2(d, e)
                if abs(det) < 1e-8:
                    continue  # parallel (no unique intersection)
                # solve a + s e = t d  =>  s = cross2(a, d)/det,  t = cross2(a, e)/det
                s = cross2(a, d) / det
                t = cross2(a, e) / det
                if 0.0 <= s <= 1.0:
                    # compute intersection in the 2D frame
                    x2, y2 = t*d[0], t*d[1]
                    # lift back to 3D:  P = O + x2*u + y2*v
                    P3 = add(O, add(mul(u, x2), mul(v, y2)))
                    hits.append((t, P3))
            # sort by t so we get the “backward” and “forward” hits
            hits.sort(key=lambda it: it[0])
            if len(hits) >= 2:
                return [hits[0][1], hits[-1][1]]
            else:
                return []

        return {
            'line1': find_hits(m1),
            'line2': find_hits(m2),
        }
        
    def addsphere(x_sphere, y_sphere, z_sphere):  
        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.05, location=(x_sphere, y_sphere, z_sphere))
        return {'FINISHED'}

def register():
    # Register the castRays class
    bpy.utils.register_class(test_script)

def unregister():
    # Register the castRays class
    bpy.utils.register_class(test_script)

if __name__ == "__main__":
    register()

