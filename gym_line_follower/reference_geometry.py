from shapely.geometry import Point, MultiPoint
from shapely.geometry.polygon import Polygon
from shapely.affinity import rotate, translate
import numpy as np


class ReferenceGeometry:

    def __init__(self, geometry, origin=(0., 0.), origin_angle=0.):
        self.origin = origin
        self.origin_angle = origin_angle

        self.position = self.origin
        self.rotation = self.origin_angle

        self.geometry = geometry
        self.plottable = self._get_plottable(geometry)

    def move(self, new_position, new_rotation):
        """
        Move to new position and orientation.
        :param new_position:
        :param new_rotation:
        :return: None
        """
        x_off = new_position[0] - self.position[0]
        y_off = new_position[1] - self.position[1]
        ang_off = new_rotation - self.rotation

        self.geometry = rotate(self.geometry, ang_off, origin=self.position, use_radians=True)
        self.geometry = translate(self.geometry, x_off, y_off)

        self.position = new_position
        self.rotation = new_rotation
        self.plottable = self._get_plottable(self.geometry)

    def _get_plottable(self, geometry):
        """
        Return tuple of x, y coordinates that can be directly plotted with pyplot.
        :param geometry: shapely.geometry instance
        :return: tuple of x, y
        """
        raise NotImplementedError


class CameraWindow(ReferenceGeometry):
    """
    Polygon of arbitrary shape that can be used for determining visibility of points.
    """

    def __init__(self, window_points, *args, **kwargs):
        """
        Init geometry.
        :param window_points: points of the polygon
        """
        self.window_points = window_points

        geometry = Polygon(window_points)
        super(CameraWindow, self).__init__(geometry, *args, **kwargs)

    def _get_plottable(self, geometry):
        return geometry.exterior.xy

    def visible_points(self, points, return_coords=True):
        if not isinstance(points, MultiPoint):
            points = np.array(points)
            multipoint = MultiPoint(points)
        else:
            multipoint = points

        visible = multipoint.intersection(self.geometry)

        if not return_coords:
            return visible

        coords = []
        if visible.is_empty:
            return np.zeros((0, 2))  # empty array
        elif isinstance(visible, Point):
            coords.append(visible.coords[0])
        elif hasattr(visible, 'geoms'):
            for geom in visible.geoms:
                if isinstance(geom, Point):
                    coords.append(geom.coords[0])
        else:
            # Unknown geometry type
            return np.zeros((0, 2))

        return np.array(coords)


    def convert_to_local(self, geom):
        """
        Convert geometry to local coordinate system.
        :param geom: geometry to be converted, shapely.geometry
        :return: geometry in local c.s.
        """
        x_off = self.origin[0] - self.position[0]
        y_off = self.origin[1] - self.position[1]
        ang_off = self.origin_angle - self.rotation
        out = rotate(geom, ang_off, origin=self.position, use_radians=True)
        out = translate(out, x_off, y_off)
        return out

    def convert_points_to_local(self, points):
        geom = self.convert_to_local(MultiPoint(points))
        if geom.is_empty:
            return np.zeros((0, 2))
        elif isinstance(geom, Point):
            return np.array([geom.coords[0]])
        elif hasattr(geom, 'geoms'):
            coords = [pt.coords[0] for pt in geom.geoms if isinstance(pt, Point)]
            return np.array(coords)
        else:
            return np.zeros((0, 2))


    def get_local_window(self):
        return self.__class__(self.window_points, self.origin, self.origin_angle)


class ReferencePoint(ReferenceGeometry):

    def __init__(self, xy_shift):
        point = Point(xy_shift)
        super(ReferencePoint, self).__init__(point)

    def _get_plottable(self, geometry):
        return geometry.x, geometry.y

    def get_xy(self):
        return self.geometry.x, self.geometry.y


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    win_points = [(0.280, 0.135), (0.280, -0.135), (0.112, -0.060), (0.112, 0.060)]
    points = np.array([(0.2, 0.2), (0.2, 0.05), (0.2, -0.05), (0.2, -0.2)])

    window = CameraWindow(win_points)
    visible_pts = window.visible_points(points)

    plt.axes().set_aspect("equal")
    plt.plot(*window.plottable, ".--")
    for p in points:
        plt.plot(*p, "g*")
    for p in visible_pts:
        plt.plot(*p, "r*")
    plt.show()
