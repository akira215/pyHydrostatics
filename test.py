from point2D import Point2D
from polygon import Polygon,Polygons

polygoneInit = [Point2D(0,0),Point2D(0.5,-1),Point2D(1.5,-0.2),Point2D(2,-0.5),Point2D(2,0),Point2D(1.5,1),Point2D(0.3,0),Point2D(0.5,1)]

#polygoneClockwise
polygoneTest = [Point2D(0.0,0.0),Point2D(0.5,1.0),Point2D(0.3,0.0),Point2D(1.5,1.0),Point2D(2.0,0.0),Point2D(2.0,-0.5),Point2D(1.5,-0.2),Point2D(0.5,-1.0)]

#shifted
polygoneTest = [Point2D(2,-0.5),Point2D(2,0),Point2D(1.5,1),Point2D(0.3,0),Point2D(0.5,1),Point2D(0,0),Point2D(0.5,-1),Point2D(1.5,-0.2)]

#carre
#polygoneTest = [Point2D(0,0.0),Point2D(2,0),Point2D(2,2),Point2D(0,2)]

#carre clockwise
#polygoneTest = [Point2D(0.0,0.0),Point2D(2,0),Point2D(2,-2.0),Point2D(0,-2.0)]

#Intial
#polygoneTest = [Point2D(0,0),Point2D(0.0,-1),Point2D(1.5,-0.2),Point2D(2,-0.5),Point2D(2,0),Point2D(1.5,1),Point2D(0.3,0.0),Point2D(0.5,1)]

#initial with only one cut
#polygoneTest = [Point2D(0,0),Point2D(0.5,-1),Point2D(1.5,-0.2),Point2D(2,-0.5),Point2D(2,0),Point2D(1.5,0),Point2D(0.3,0.2),Point2D(0.5,1)]

l = [Point2D(-0.5,0.6),Point2D(2.7,-0.75)]
#l = [Point2D(-0.5,0),Point2D(2.7,0)]
p = Polygon(polygoneTest)
print ('list',p.triangulate_polygon())
print ('area',p.getArea())
print ('CoG',p.getCog())
print ('poly',p)
#print('split',p.splitPolygon(l))
polytotal = Polygons(p.splitPolygon(l))
print ('areaTotal',polytotal.getArea())
print ('CoG total',polytotal.getCog())
print ('cutted edges',p.edges)
print('********Polyfromside**************')
pSide = Polygons(p.getPolygonsFromSide(l,1))
print('getfromside',pSide)
print ('areaTotal',pSide.getArea())
print ('CoG total',pSide.getCog())
print ('cutted edges',pSide.edges)
print('***************************')
print ('poly',p)
print ('poly safe',p)

