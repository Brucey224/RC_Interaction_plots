

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