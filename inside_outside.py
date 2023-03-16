import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from scipy.cluster.hierarchy import fclusterdata

def inside_outside(cell_points):
    """
    Define a flat log-prior for use with the PTsampler

    Parameters:
    -----------
    cell_points : ndarray of floats, shape (npoints, 3)
        An array whose rows define the locations in 3D of npoints cells.

    Returns:
    --------
    outside_cells : ndarray of bools
        An array whose nth element is True if the cell whose location is given 
        by the nth row of cell_points is classified as on the outside of the set 
        of cells, and is False otherwise.

    cluster_data : ndarray of floats, shape (npoints, 2)
        An array whose nth row comprises two elements: (i) the quantity given by 
        log10(m + 0.01), where m denotes the minimum distance from the cell 
        whose location is given by the nth row of cell_points to the surface of 
        the set of cells; (ii) the variance in distances from this cell to the 
        surface of the set of cells.
    """
    # Construct the convex hull of the Delaunay triangulation of the set of cells
    tri = Delaunay(cell_points)
    hull = ConvexHull(tri.points)
    hull_points = hull.points
    hull_simplices = hull.simplices

    # Get the unit outward normal to each face of the convex hull
    unit_outward_normals = hull.equations[:,:-1]

    num_faces = len(hull_simplices)
    num_cells = len(cell_points)
    dist_from_cells_to_faces = np.inf * np.ones((num_cells, num_faces))

    u = hull_points[hull_simplices[:,1]] - hull_points[hull_simplices[:,0]]
    v = hull_points[hull_simplices[:,2]] - hull_points[hull_simplices[:,0]]
    n = np.cross(u, v)
    n = n / np.sqrt(n[0]*n[0] + n[1]*n[1])
    n_cross_u = np.cross(n, u)
    v_cross_n = np.cross(v, n)

    # Get the distance from each cell to each face's coplanar surface
    for cell_idx in range(num_cells):
        
        cell_point = cell_points[cell_idx]
        w = cell_point - hull_points[hull_simplices[:,0]]
        
        gamma = np.einsum('ij,ij->i', n_cross_u, w)
        beta = np.einsum('ij,ij->i', v_cross_n, w)
        alpha = 1 - gamma - beta

        for face_idx in range(num_faces):

            # Determine if the proj(X0) is in the triangular domain via barycentric coordinates
            # See https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
            if ((0 <= alpha[face_idx] <= 1) and (0 <= beta[face_idx] <= 1) and (0 <= gamma[face_idx] <= 1)):
                dist_from_cells_to_faces[cell_idx, face_idx] = np.dot(unit_outward_normals[face_idx], hull_points[hull_simplices[face_idx][0]] - cell_points[cell_idx])
            
            else:
                # Projection of cell onto plane containing this face
                P = alpha[face_idx]*hull_points[hull_simplices[face_idx,0]] + beta[face_idx]*hull_points[hull_simplices[face_idx,1]] + gamma[face_idx]*hull_points[hull_simplices[face_idx,2]]

                # For each edge of this face
                for vertex in range(3):

                    A = hull_points[hull_simplices[face_idx][vertex]]
                    B = hull_points[hull_simplices[face_idx][(vertex+1)%3]]
                    
                    # Let AB be a line segment in 3D, and P be the projection of 
                    # a point X onto a plane containing this line segment. If 
                    # angle PAB is > 90 degrees, then A is the closest point on 
                    # AB to P, and the distance from X to AB is |XA|. Else, if 
                    # angle PBA is > 90 degrees, then B is the closest point on 
                    # AB to P, and the distance from X to AB is |XB|. Else, the 
                    # closest point on AB to P is C = A - BA*(PA.AB)/|AB|^2, and 
                    # the distance from X to AB is |XC|.
                    if np.dot(P - A, B - A) < 0:
                        dist_from_cell_to_edge = np.linalg.norm(cell_point - A) 
                    elif np.dot(P - B, A - B) < 0:
                        dist_from_cell_to_edge = np.linalg.norm(cell_point - B)
                    else:
                        C = A - (A - B)*np.dot(A - P, B - A) / (np.dot(B - A, B - A))
                        dist_from_cell_to_edge = np.linalg.norm(cell_point - C)
                    
                    dist_from_cells_to_faces[cell_idx, face_idx] = min(dist_from_cells_to_faces[cell_idx, face_idx], dist_from_cell_to_edge)

    # Calculate (and transform) the minimum distance from each cell to the convex hull
    min_distances = np.min(dist_from_cells_to_faces, axis=1)
    transformed_min_distances = np.log10(min_distances + 0.01)

    # Calculate the variance in distances from each cell to the convex hull
    var_distances = np.var(dist_from_cells_to_faces, axis=1)

    # Perform clustering, having scaled (transformed) mins and variances to lie in [0, 1]
    cluster_data = np.transpose([transformed_min_distances, var_distances])
    cluster_data = cluster_data - np.min(cluster_data, axis=0)
    cluster_data = cluster_data / np.max(cluster_data, axis=0)

    cluster_indices = fclusterdata(cluster_data, t=2, criterion='maxclust', method='ward')

    # Find the mean (transformed) minimum-distance among cells in each cluster
    mean_min_cluster_1 = np.mean(cluster_data[np.where(cluster_indices == 1)][:,0])
    mean_min_cluster_2 = np.mean(cluster_data[np.where(cluster_indices == 2)][:,0])
    outside_cluster = 1 if (mean_min_cluster_1 < mean_min_cluster_2) else 2

    outside_cells = np.where(cluster_indices == outside_cluster)

    return outside_cells, cluster_data
