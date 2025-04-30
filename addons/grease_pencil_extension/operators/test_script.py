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

            grad1 = test_script.gradientcalc((vert1[0], vert1[1]), (vert2[0], vert2[1]))
            grad2 = test_script.gradientcalc((vert3[0], vert3[1]), (vert4[0], vert4[1]))    
            print(grad1)
            print(grad2)

            # --- Example Usage ---

            quad = [
                (vert1[0], vert1[1], vert1[2]),
                (vert2[0], vert2[1], vert2[2]),
                (vert3[0], vert3[1], vert3[2]),
                (vert4[0], vert4[1], vert4[2])
            ]
            desiredGrad = 0.2
            dummyGrad = test_script.convertGrad(grad1, desiredGrad)
            print(dummyGrad)
            P = (0, 0, 0)  # interior point3
            m1, m2 = dummyGrad, 10
            hits = test_script.intersect_lines_with_quad(quad, P, m1, m2)
            print(hits)
            line1 = (hits['line1'])
            line2 = (hits['line2']) 
            line1_Vert1 = line1[0]
            line1_Vert2 = line1[1]
            line2_Vert1 = line2[0]
            line2_Vert2 = line2[1]
            sphere1 = (line1_Vert1[0], line1_Vert1[1], line1_Vert1[2])
            sphere2 = (line1_Vert2[0], line1_Vert2[1], line1_Vert2[2])
            #sphere3 = (line2_Vert1[0], line2_Vert1[1], line2_Vert1[2])
            #sphere4 = (line2_Vert2[0], line2_Vert2[1], line2_Vert2[2])
        
            test_script.addsphere(P)
            test_script.addsphere(sphere1)
            test_script.addsphere(sphere2)
            #test_script.addsphere(sphere3)
            #test_script.addsphere(sphere4)
    
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
        
    def addsphere(coords):  
        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.02, location=(coords[0], coords[1], coords[2]))
        return {'FINISHED'}
    
    def convertGrad(TrueGrad, desiredGrad):  
        # Helper function to calculate the dummy grad to dope the transformed x-y axis 
        # Converts Desired grad and actual grad to degrees sums them and returns the result as gradient
        GradConvert = -1*TrueGrad
        TrueGradRad = math.atan(GradConvert)
        TrueGradDegrees = math.degrees(TrueGradRad)
        desiredGradRad = math.atan(desiredGrad)
        desiredGradDegrees = math.degrees(desiredGradRad)
        summedGradDegrees = TrueGradDegrees+desiredGradDegrees
        finalGrad = math.tan(math.radians(summedGradDegrees))
        return finalGrad

    def gradientcalc(p1, p2):
        """
        Calculate the gradient (slope) of the line through two 2D points.

        Parameters:
            p1 (tuple of float): (x1, y1)
            p2 (tuple of float): (x2, y2)

        Returns:
            float: (y2 - y1) / (x2 - x1)

        Raises:
            ValueError: if x2 == x1 (vertical line has undefined slope).
        """
        x1, y1 = p1
        x2, y2 = p2

        dx = x2 - x1
        if dx == 0:
            raise ValueError(f"Vertical line through x = {x1!r}; slope is undefined.")

        return (y2 - y1) / dx
    
    def intersect_lines_3d(p1, d1, p2, d2, tol=1e-6):
        """
        Given two lines in 3D:
        L1: X = p1 + s * d1
        L2: X = p2 + t * d2
        where p1, p2 are points (x,y,z) and d1, d2 are direction vectors (dx,dy,dz),
        determine if they intersect. If so, return the intersection point as (x,y,z).
        Otherwise return None.

        Parameters:
            p1 (tuple of float): a point on line1
            d1 (tuple of float): direction vector of line1
            p2 (tuple of float): a point on line2
            d2 (tuple of float): direction vector of line2
            tol (float): tolerance for floating-point comparisons

        Returns:
            tuple of float or None: intersection point, or None if no intersection.
        """

        # vector ops
        def sub(a,b):   return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
        def add(a,b):   return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
        def mul(a, s):  return (a[0]*s,   a[1]*s,   a[2]*s)
        def dot(a,b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        def cross(a,b):
            return (a[1]*b[2] - a[2]*b[1],
                    a[2]*b[0] - a[0]*b[2],
                    a[0]*b[1] - a[1]*b[0])
        def norm(a):    return math.sqrt(dot(a,a))

        # 1) Check if directions are parallel
        cr = cross(d1, d2)
        if norm(cr) < tol:
            # parallel (or colinear)
            # check if (p2 - p1) is also parallel to d1 → colinear
            if norm(cross(sub(p2, p1), d1)) < tol:
                raise ValueError("Lines are colinear (infinitely many intersections).")
            else:
                return None  # strictly parallel, no intersection

        # 2) To solve p1 + s d1 = p2 + t d2, pick the two coords where cross(d1, d2) has largest magnitude
        #    so the projection yields a nonsingular 2×2 system.
        abs_cr = list(map(abs, cr))
        # index of the largest component in cross
        k = abs_cr.index(max(abs_cr))
        # use the other two axes
        idx = [0,1,2]
        idx.remove(k)  # now idx = [i, j]

        i, j = idx

        # Build and solve the 2×2 system:
        #    p1[i] + s*d1[i] = p2[i] + t*d2[i]
        #    p1[j] + s*d1[j] = p2[j] + t*d2[j]
        A = (( d1[i], -d2[i] ),
            ( d1[j], -d2[j] ))
        B = ( p2[i] - p1[i],
            p2[j] - p1[j] )

        detA = A[0][0]*A[1][1] - A[0][1]*A[1][0]
        if abs(detA) < tol:
            return None  # degenerate in this projection plane

        # Cramer’s rule
        s = ( B[0]*A[1][1] - B[1]*A[0][1] ) / detA
        t = ( A[0][0]*B[1] - A[1][0]*B[0] ) / detA

        # 3) Compute the two candidate points and check they coincide
        X1 = add(p1, mul(d1, s))
        X2 = add(p2, mul(d2, t))
        if norm(sub(X1, X2)) > tol:
            return None  # they miss (skew lines)

        # 4) Success!
        return X1

def register():
    # Register the castRays class
    bpy.utils.register_class(test_script)

def unregister():
    # Register the castRays class
    bpy.utils.register_class(test_script)

if __name__ == "__main__":
    register()

