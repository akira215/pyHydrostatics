import math

class Point2D:  
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return "Point2D(%s,%s)"%(self.x, self.y)

    def __repr__(self):
        return str(self)
    
    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        return Point2D(x, y)
    
    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        return Point2D(x, y)
    
    def __mul__ (self, other:float):
        return Point2D(self.x * other, self.y * other)

    def dot_product(self, other):
        return self.x * other.x + self.y * other.y

    def distanceFromOrigin(self):
        origin = Point2D(0, 0)
        dist = self.distanceFromPoint(self, origin)
        return dist
  
    def distanceTo(self, point):
        return  math.sqrt((self.x - point.x)**2 + (self.y - point.y)**2)
    
