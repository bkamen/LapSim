from bs4 import BeautifulSoup

kml_file = 'LeMans_left.kml'


def process_coordinate_read(coords):
    lat = []
    lon = []
    for i in range(1,len(coords)):
        lat.append(coords[i].split(',')[0])
        lon = coords[i].split(',')[1]


with open(kml_file, 'r') as f:
    s = BeautifulSoup(f, 'xml')

    coords = s.find_all('coordinates')
    coords = str(coords).split()
    print(coords)

