from shapely.geometry import Polygon, LineString
import pandas as pd
import numpy as np
from shapely.ops import split
import math
from scipy.integrate import quad

def arrange_circular_bars(diameter, cover, radial_num_bars):
    rebar_coords = []
    if radial_num_bars > 0 and cover > 0:
        radius = diameter / 2
        theta = 2 * math.pi / radial_num_bars
        bar_radius = radius - cover
        x = [bar_radius * math.sin(theta * i) for i in range(radial_num_bars)]
        y = [bar_radius * math.cos(theta * i) for i in range(radial_num_bars)]
    rebar_coords = np.array([x,y]).T
    return rebar_coords

def arrange_rectangular_bars(b, h, cover, link_dia, phi, n_x, n_y):
    rebar_coords = []
    cover = float(cover)
    b = float(b)
    h = float(h)
    link_dia = float(link_dia)
    phi = float(phi)
    n_x = int(n_x)
    n_y = int(n_y)
    spacing_x = float((b - 2*cover - 2*link_dia - phi) / (n_x - 1))
    spacing_y = float((h - 2*cover - 2*link_dia - phi) / (n_y - 1))
        
    #create rebar in top + bottom layer
    for i in range(n_x):
        x1 = cover + link_dia + phi/2 + i*spacing_x
        y1 = cover + link_dia + phi/2
        rebar_coords.append([x1, y1])
        x2 = x1
        y2 = h - cover - link_dia - phi/2
        rebar_coords.append([x2,y2])
    # create side bars
    for j in range(n_y-2):
        x1 = cover + link_dia + phi/2
        y1 = cover + link_dia + phi/2 + spacing_y * (j+1)
        rebar_coords.append([x1, y1])
        x2 = b - cover - link_dia - phi/2
        y2 = y1
        rebar_coords.append([x2,y2])
    return np.array(rebar_coords)

def clip_polygon_at_y(polygon, y_value):
    # Create a horizontal line at y_value
    line = LineString([(-1e10, y_value), (1e10, y_value)])
    # Split the polygon by the line
    result = split(polygon, line)
    #Return multiple geoemtry objects
    return result.geoms

# Function to clip polygon at a certain x value
def clip_polygon_at_x(polygon, x_value):
    # Create a horizontal line at x_value
    line = LineString([(x_value, -1e10), (x_value, 1e10)])
    # Split the polygon by the line
    result = split(polygon, line)
    #Return multiple geoemtry objects
    return result.geoms

def integrate_part_of_circle(neutral_axis, diameter, lambd):
    r = diameter / 2 
    distance_from_centre = diameter / 2 - lambd * neutral_axis
    if distance_from_centre < diameter/2:
        # Calculate the area using the integral
        result, error = quad(lambda x, r: np.sqrt(r**2 - distance_from_centre**2), distance_from_centre, r, args=(r,))
        area = 2 * result
        # Calculate the angle for the centroid calculation
        theta = 2 * np.arccos(distance_from_centre/ r)
        area = r**2/2 *(theta - np.sin(theta))
        # Calculate the centroid y-coordinate (only y because x is 0 / not applicable)
        centroid = diameter/2 - (4 * r * (np.sin(theta / 2) ** 3)) / (3 * (theta - np.sin(theta)))
        lever_arm = neutral_axis - centroid
    else:
        area = math.pi*diameter**2/4
        lever_arm = neutral_axis - diameter/2
        centroid_y = diameter/2
    return area, centroid, lever_arm

def compute_area_and_centroid(polygon):
    # Get the coordinates of the polygon
    x, y = polygon.exterior.coords.xy
    # Convert to numpy arrays for easy calculations
    x = np.array(x)
    y = np.array(y)
    # Area Calculation
    A = 0.5 * np.sum(y[:-1] * x[1:] - y[1:] * x[:-1])
    # Centroid Calculation
    Cx = np.sum((x[:-1] + x[1:]) * (y[:-1] * x[1:] - y[1:] * x[:-1])) / (6*A)
    Cy = np.sum((y[:-1] + y[1:]) * (y[:-1] * x[1:] - y[1:] * x[:-1])) / (6*A)
    return A, (Cx, Cy)

def compute_second_moment_area(polygon: Polygon):
    c_x, c_y = polygon.centroid.x, polygon.centroid.y
    coords = list(polygon.exterior.coords)
    Ix = 0.0
    Iy = 0.0
    Ixy = 0.0

    for i in range(len(coords) - 1):
        x0, y0 = coords[i]
        x1, y1 = coords[i + 1]

        Ix += 1/12* (y0 - y1)*(x1+x0-2*c_x)*(x1**2+x0**2+2*c_x*(c_x-x0-x1))        
        Iy += 1/12* (x0 - x1)*(y1+y0-2*c_y)*(y1**2+y0**2+2*c_y*(c_y-y0-y1))
    return Ix, Iy

def rotate_points_3d(points, angle, help=False):
    if help:
            print('''
            Rotates a set of 2D points around the Z axis by a given
                  Inputs:   points - array of points to be rotated
                            angle - angle in degrees to rotate points (rotates anticlockwise by default)
                  Outputs:  rotated points
                  ''')
    # Rotate a set of points around the z-axis by a given angle
    theta = np.radians(angle)
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
    return np.dot(points, rotation_matrix.T)

