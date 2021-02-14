import numpy as np
import scipy
from scipy.ndimage import laplace
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from matplotlib import cm, animation, rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML
import time
import os
from matplotlib.backends.backend_pdf import PdfPages



def ficks(s, w):
    return(laplace(s) / np.power(w, 2))


def hill(s, K, lam):
    s[s < 0] = 0
    h = s**lam / (K**lam + s**lam)
    return(h)

# Results from Clemens' data use these parameters in your model with the hill_AHL function for AHL production
kd = 17.8
n = 1.75
min =  4230
max = 54096
def hill_AHL(conc, n, kd, min, max):
    # get rid of the very small negatvie values
    conc[conc<0] = 0


    h = (min + (max-min)*(conc**n/(kd**n + conc**n)))/max

    return h

def multi_plots(sim, title=""):
    f, ax = plt.subplots(3, 3, sharex=True, sharey=False, figsize=(15, 15))

    f.suptitle(title, fontsize=40)
    im1 = ax[0, 0].imshow(sim[3], interpolation="none", cmap=cm.viridis, vmin=0, vmax=1)
    ax[0, 0].set_title("Sender")
    divider = make_axes_locatable(ax[0, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im1, cax=cax, shrink=0.8)

    im2 = ax[0, 1].imshow(sim[5], interpolation="none", cmap=cm.viridis, vmin=0, vmax=1)
    ax[0, 1].set_title("Receiver")
    divider = make_axes_locatable(ax[0, 1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im2, cax=cax, shrink=0.8)

    im3 = ax[0, 2].imshow(sim[1], interpolation="none", cmap=cm.viridis, vmin=0, vmax=5)
    ax[0, 2].set_title("Arabinose")
    divider = make_axes_locatable(ax[0, 2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im3, cax=cax, shrink=0.8)

    im4 = ax[1, 0].imshow(sim[3], interpolation="none", cmap=cm.viridis, vmin=0)
    ax[1, 0].set_title("LuxI")
    divider = make_axes_locatable(ax[1, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im4, cax=cax, shrink=0.8)

    im5 = ax[1, 1].imshow(sim[4], interpolation="none", cmap=cm.viridis, vmin=0)
    ax[1, 1].set_title("C6")
    divider = make_axes_locatable(ax[1, 1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im5, cax=cax, shrink=0.8)

    im6 = ax[1, 2].imshow(sim[6], interpolation="none", cmap=cm.viridis, vmin=0)
    ax[1, 2].set_title("GFP")
    divider = make_axes_locatable(ax[1, 2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im6, cax=cax, shrink=0.8)

    im7 = ax[2, 0].imshow(sim[2], interpolation="none", cmap=cm.viridis, vmin=0, vmax=100)
    ax[2, 0].set_title("Nutrients")
    divider = make_axes_locatable(ax[2, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im7, cax=cax, shrink=0.8)

    ax[2, 1].axis('off')
    ax[2, 2].axis('off')

    return(f)

def plot_nn(sim, title=""):
    f, ax = plt.subplots(3, 3, sharex=True, sharey=False, figsize=(15, 15))

    f.suptitle(title, fontsize=40)
    im1 = ax[0, 0].imshow(sim[3], interpolation="none", cmap=cm.viridis, vmin=0, vmax=1)
    ax[0, 0].set_title("Input Layer")
    divider = make_axes_locatable(ax[0, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im1, cax=cax, shrink=0.8)

    im2 = ax[0, 1].imshow(sim[5], interpolation="none", cmap=cm.viridis, vmin=0, vmax=1)
    ax[0, 1].set_title("Hidden Layer")
    divider = make_axes_locatable(ax[0, 1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im2, cax=cax, shrink=0.8)

    im3 = ax[0, 2].imshow(sim[9], interpolation="none", cmap=cm.viridis, vmin=0, vmax=1)
    ax[0, 2].set_title("Output Layer")
    divider = make_axes_locatable(ax[0, 2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im3, cax=cax, shrink=0.8)

    im5 = ax[1, 0].imshow(sim[4], interpolation="none", cmap=cm.viridis, vmin=0 )
    ax[1, 0].set_title("AHL 1")
    divider = make_axes_locatable(ax[1, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im5, cax=cax, shrink=0.8)

    im5 = ax[1, 1].imshow(sim[8], interpolation="none", cmap=cm.viridis, vmin=0)
    ax[1, 1].set_title("AHL 2")
    divider = make_axes_locatable(ax[1, 1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im5, cax=cax, shrink=0.8)

    im6 = ax[1, 2].imshow(sim[6], interpolation="none", cmap=cm.viridis, vmin=0)
    ax[1, 2].set_title("GFP")
    divider = make_axes_locatable(ax[1, 2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im6, cax=cax, shrink=0.8)

    im7 = ax[2, 0].imshow(sim[2], interpolation="none", cmap=cm.viridis, vmin=0, vmax=100)
    ax[2, 0].set_title("Nutrients")
    divider = make_axes_locatable(ax[2, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im7, cax=cax, shrink=0.8)

    im7 = ax[2, 1].imshow(sim[1], interpolation="none", cmap=cm.viridis, vmin=0)
    ax[2, 1].set_title("Arabinose")
    divider = make_axes_locatable(ax[2, 1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im7, cax=cax, shrink=0.8)

    ax[2, 1].axis('off')
    ax[2, 2].axis('off')

    return(f)


def plots(sim, names):

    n_plots = sim.shape[0]
    x = int(np.ceil(n_plots / 3))

    f, ax = plt.subplots(x, 3, sharex=True, sharey=False, figsize=(15, 15))

    for i, val in enumerate(ax.flatten()):

        if i < n_plots:
            im = ax.flatten()[i].imshow(sim[i], cmap=cm.viridis, vmin=0)
            ax.flatten()[i].set_title(names[i])
            divider = make_axes_locatable(ax.flatten()[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            f.colorbar(im, cax=cax, shrink=0.8, label='')
        else:
            ax.flatten()[i].axis("off")



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


