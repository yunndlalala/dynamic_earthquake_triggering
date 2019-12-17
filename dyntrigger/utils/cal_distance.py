from math import radians, cos, sin, asin, sqrt


def haversine(lon1, lat1, lon2, lat2):  # degree
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Translate degree to radian.
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula.
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of the Earth.
    return c * r, c * r / 111


if __name__ == "__main__":
    print(haversine(103, 26.5, -76.36, 1.93))
    print('finish')
