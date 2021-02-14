from scipy.ndimage import laplace
import numpy as np


def leaky_hill(s, K, lam, min, max):
    # get rid of the very small negative values
    s[s < 0] = 0

    h = (max - min) * s ** lam / (K ** lam + s ** lam) + min
    return h

def leaky_inverse_hill(s, K, lam, min, max):
    # get rid of the very small negative values
    s[s < 0] = 0

    h = (max - min) * K ** lam / (K ** lam + s ** lam) + min
    return h

def hill(s, K, lam):
    s[s < 0] = 0
    h = s**lam / (K**lam + s**lam)
    return(h)


def ficks(s, w):
    return(laplace(s) / np.power(w, 2))

def get_vertex_coordinates(vertex_numbers, n_rows, n_cols):
    '''
    use to get grid coordinates of vertices

    args:
        vertex_numbers: the numbers of the vertices you want coordinates for 0 <= vertex_number < n_rows * n_cols
        n_rows, n_cols: number of rows and columns in the finite difference simulation, a total of n-rows*n_cols vertices

    returns:
        vertex_coordinates: the coordinates on the finite difference grid of the supplied vertex number: [[r0, c0]; [r1,c1]; ... [rn,cn]]
            these use matrix indexing, in the format (row, col) starting from the top left of the grid
    '''

    vertex_coordinates = np.hstack((vertex_numbers // n_rows, vertex_numbers % n_cols))

    return vertex_coordinates

def get_vertex_positions(vertex_numbers, n_rows, n_cols, w):
    '''
    use to get the positions (in mm) of vertices on the real grid

    args:
        vertex_numbers: the numbers of the vertices you want coordinates for 0 <= vertex_number < n_rows * n_cols
        n_rows, n_cols: number of rows and columns in the finite difference simulation, a total of n-rows*n_cols vertices
        w: the distance between finite difference vertices
    returns:
        vertex_positions: the positions on the finite difference grid of the supplied vertex number (in mm from the top left of the grid):
            [[r0, c0]; [r1,c1]; ... [rn,cn]]
    '''


    vertex_coordinates = get_vertex_coordinates(vertex_numbers, n_rows, n_cols)

    vertex_positions = vertex_coordinates * w

    return vertex_positions


def assign_vertices(vertex_positions, node_positions, node_radius):
    '''
    assigns vertices to be part of nodes in node_positions with radius: node radius.


    args:
        vertex_positions: the positions of the vertices to be tested
        node_positions, node_radius: positions and radius of the nodes we want vertices for
    returns:
        vertex_numbers: the numbers of the vertices that are within on of the nodes
        indicators: vector with an index for each vertex indicating whether it is inside a node (value = 1) or outside all nodes (value = 0)

     NOTE: this assigns position based on real life position, not the grid coordinates i.e the distance in mm
    '''


    indicators = np.zeros(len(vertex_positions))

    if node_positions == []:
        return [], indicators

    if node_positions[0] is not None:
        node_positions = np.array(node_positions)
        differences = vertex_positions - node_positions[:, None]

        vertex_numbers = np.where(np.linalg.norm(differences, axis=2) < node_radius)[1].reshape(-1, 1)

        indicators[vertex_numbers] = 1

    indicators = np.array(indicators, dtype=np.int32)



    return vertex_numbers, indicators

# this is the only one you really need to use
def get_node_coordinates(node_positions, node_radius, n_rows, n_cols, w):
    '''
       gets the coordinates of the vertices inside the nodes with position node_positions with radius: node radius.

       args:
           vertex_positions: the positions of the vertices to be tested
           node_positions, node_radius: positions and radius of the nodes we want vertices for
           n_rows, n_cols: the number of rows and cols on the finite difference grid
       returns:
           coordinates: the coordinates of the vertices that are within on of the nodes

        NOTE: this assigns position based on real life position, not the grid coordinates i.e the distance in mm
       '''


    # use the individual functions if repeating these two lines for each node type is too slow
    all_vertex_numbers = np.arange(n_rows * n_cols).reshape(-1, 1)  # reshpae to colum vector
    all_vertex_positions = get_vertex_positions(all_vertex_numbers, n_rows, n_cols, w)

    vertex_numbers, vertex_indicators = assign_vertices(all_vertex_positions, node_positions, node_radius)
    coordinates = get_vertex_coordinates(vertex_numbers, n_rows, n_cols)

    return coordinates
