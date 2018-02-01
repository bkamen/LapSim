from bs4 import BeautifulSoup
import pyproj
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate

p = pyproj.Proj(proj='utm', ellps='WGS84')  # use kwargs

kml_file = ['LeMans_left.kml', 'LeMas_right.kml']


def process_coordinate_read(coords):
    lat = []
    lon = []
    for i in range(1,len(coords)-1):
        #print(coords[i].split(','))
        lon.append(coords[i].split(',')[0])
        lat.append(coords[i].split(',')[1])
    return lat, lon


def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial

    window_size = np.abs(np.int(window_size))
    order = np.abs(np.int(order))

    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def process_kml_file(file):
    with open(file, 'r') as f:
        s = BeautifulSoup(f, 'xml')

        coords = s.find_all('coordinates')
        coords = str(coords).split()
        #print(coords)

        lat, lon = process_coordinate_read(coords)

    x = np.array([])
    y = np.array([])

    for i, j in zip(lon, lat):
        x_, y_ = p(i, j)
        x = np.append(x, x_)
        y = np.append(y, y_)

    dist = scipy.integrate.cumtrapz(np.abs(np.sqrt(np.gradient(x)**2+np.gradient(y)**2)))
    dist = np.append(np.array([0]), dist)
    dist_new = np.arange(0, dist[-1])
    x = np.interp(dist_new, dist, x)
    y = np.interp(dist_new, dist, y)

    x = savitzky_golay(x, window_size=40, order=4)
    y = savitzky_golay(y, window_size=40, order=4)

    #plt.plot(x, y)
    #print(dist[-1])

    return x, y


def find_centerline(x1, y1, x2, y2):
    x_ct = np.array([])
    y_ct = np.array([])

    x_left = np.array([])
    x_right = np.array([])
    y_left = np.array([])
    y_right = np.array([])

    width = np.array([])
    for i, j in zip(x1, y1):
        dist = np.sqrt((i-x2)**2+(j-y2)**2)
        idx_nb = np.argmin(dist) # neighbour index
        x_ct = np.append(x_ct, i-0.5*(i-x2[idx_nb]))
        y_ct = np.append(y_ct, j-0.5*(j-y2[idx_nb]))
        x_left = np.append(x_left, i)
        y_left = np.append(y_left, j)
        x_right = np.append(x_right, x2[idx_nb])
        y_right = np.append(y_right, y2[idx_nb])
        width = np.append(width, np.sqrt((i-x2[idx_nb])**2+(j-y2[idx_nb])**2))


    dist = scipy.integrate.cumtrapz(np.abs(np.sqrt(np.gradient(x_ct) ** 2 + np.gradient(y_ct) ** 2)))
    dist = np.append(np.array([0]), dist)
    dist_new = np.arange(0, dist[-1])
    x_ct = np.interp(dist_new, dist, x_ct)
    y_ct = np.interp(dist_new, dist, y_ct)
    x_left = np.interp(dist_new, dist, x_left)
    y_left = np.interp(dist_new, dist, y_left)
    x_right = np.interp(dist_new, dist, x_right)
    y_right = np.interp(dist_new, dist, y_right)
    return x_ct, y_ct, x_left, y_left, x_right, y_right, dist_new, width


def get_curvature(x, y):
    dx = np.gradient(x)
    dy = np.gradient(y)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)

    curv = (dx * ddy - dy * ddx) / (dx**2 + dy**2)**1.5

    return curv


def cost_function(x,y):
    curv = get_curvature(x, y)
    rms_curv = np.sum(curv**2)

    dist = np.append(np.array([0]), scipy.integrate.cumtrapz(np.abs(np.sqrt(np.gradient(x) ** 2 + np.gradient(y) ** 2))))

    cost = rms_curv*100 + dist[-1]

    print(rms_curv*5e4)
    print(dist[-1])

    return cost


x1, y1 = process_kml_file(kml_file[0])
x2, y2 = process_kml_file(kml_file[1])
xct, yct, x_left, y_left, x_right, y_right, dist, width = find_centerline(x1, y1, x2, y2)

print(width)

curv = get_curvature(xct, yct)

# lower boundary for x
x_ub = np.fmax(x_left, x_right)
y_ub = np.fmax(y_left, y_right)

x_lb = np.fmin(x_left, x_right)
y_lb = np.fmin(y_left, y_right)


plt.plot(x_left, y_left)
plt.plot(x_right, y_right)
plt.plot(xct, yct)
#plt.plot(x_ub, y_ub, 'o')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#plt.plot(dist, curv)
#plt.show()

