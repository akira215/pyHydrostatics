import math
from point2D import Point2D

def sign(x:float)->int:
    '''return 0 if 0.0,-1 for negative values, +1 for positive ones'''
    if math.isclose( x, 0.0,abs_tol=1e-10):
        return 0
    if x > 0:
        return 1
    return -1

def nextVertex(n:int,i:int,di:int)->int:
    ''' return the circular next vertex '''
    return (i+di)%n

def MtoSegment(P0:Point2D,P1:Point2D,M:Point2D)->float:
    ''' Compute Lz from Point M relative to segment P0,P1
    This allow to know from which side from the segment is M and 
    compare distances of several points to this segment'''
    return (P1.x-P0.x)*(M.y-P0.y)-(P1.y-P0.y)*(M.x-P0.x)

def IsInTriangle(triangle:list,M:Point2D):
    '''Return True if M is striclty inside the triangle
    Triangle points shall be given counterclockwise'''
    P0 = triangle[0]
    P1 = triangle[1]
    P2 = triangle[2]
    return MtoSegment(P0,P1,M) > 0 and MtoSegment(P1,P2,M) > 0 and MtoSegment(P2,P0,M) > 0

def farthestVertex(polygon:list,P0:Point2D,P1:Point2D,P2:Point2D,indices)->int:
    '''Return the indice of the polygon vertex from triangle P0,P1,P2 
    which is the farthest from segment P1,P2'''
    n = len(polygon)
    distance = 0.0
    j = None
    for i in range(n):
        if not(i in indices):
            M = polygon[i]
            if IsInTriangle([P0,P1,P2],M):
                d = abs(MtoSegment(P1,P2,M))
                if d > distance:
                    distance = d
                    j = i
    return j

def mostLeftVertex(polygon:list)->int:
    '''Return the indice of the vertex the most on the left of the polygon.
    If more than one vertice have the same abscisse, any of them could be return'''
    n = len(polygon)
    x = polygon[0].x
    j = 0
    for i in range(1,n):
        if polygon[i].x < x:
            x = polygon[i].x
            j = i
    return j

def newPolygon(polygon:list,start:int,end:int)->list:
    ''' Generate a polygone from indice start to end, 
    considering the cyclic condition'''
    n = len(polygon)
    p = []
    i = start
    while i!=end:
        p.append(polygon[i])
        i = nextVertex(n,i,1)
    p.append(polygon[end])
    return p

def computeInter(edge:list,line:list)->Point2D:
    '''Compute Intersection with a segment'''
    ptA = edge[0]
    ptB = edge[1]
    ptC = line[0]
    ptD = line[1]

    vAB = ptB - ptA
    vWL = ptD - ptC
    div = vAB.x * vWL.y - vAB.y * vWL.x
    m = 0.0
    k = 0.0
    if (div != 0):
        m = (vAB.x * ptA.y - vAB.x * ptC.y - vAB.y * ptA.x + vAB.y * ptC.x ) / div
        k = (vWL.x * ptA.y - vWL.x * ptC.y - vWL.y * ptA.x + vWL.y * ptC.x ) / div
        if (0.0<=m) and (m<1) and (0<=k) and (k<1):
            #there is an intersection
            ptInt = ptA +  vAB * k
            return ptInt
    #No intersection
    return None


class Polygon:
    ''' Object to compute areas, cog and splitting, 
    by triangulating the polygon. This one will be counterclockwie
    oriented as soon as computation will be required'''
    def __init__(self, vertices:list):
        self.vertices = vertices
        self.trianglesList = []
        self.area = None
        self.perimeter = None
        self.cog = None
        self.splitPoly = None
        self.intersections = None 
        self.edges = None #list of segment of the cutting line
        self.isCounterClockwise = False
    
    def __str__(self):
        strPoly = 'Polygon['
        first = True
        for p in self.vertices:
            if first:
                first = False
            else:
               strPoly +=', '     
            strPoly += str(p)
        strPoly +=']' 
        return strPoly

    def __repr__(self):
        return str(self)

    def __triangulate_polygon_recursive(self,polygon:list,trianglesList:list)->list:
        '''Recursive function to divide polygone in 2 parts, starting by the most left vertex
        call recursively new polygones until getting triangles, and fill the triangle list'''
        n = len(polygon)
        
        j0 = mostLeftVertex(polygon)
        j1 = nextVertex(n,j0,1)
        j2 = nextVertex(n,j0,-1)
        P0 = polygon[j0]
        P1 = polygon[j1]
        P2 = polygon[j2]
        j = farthestVertex(polygon,P0,P1,P2,[j0,j1,j2])
        if j==None:
            trianglesList.append([P0,P1,P2])
            polygon_1=newPolygon(polygon,j1,j2)
            if len(polygon_1)==3:
                trianglesList.append(polygon_1)
            elif len(polygon_1)==2:
                #the polygon was a triangle not more to do
                return trianglesList
            else:
                self.__triangulate_polygon_recursive(polygon_1,trianglesList)
        else:
            polygon_1 = newPolygon(polygon,j0,j)
            polygon_2 = newPolygon(polygon,j,j0)
            if len(polygon_1)==3:
                trianglesList.append(polygon_1)
            else:
                self.__triangulate_polygon_recursive(polygon_1,trianglesList)
            if len(polygon_2)==3:
                trianglesList.append(polygon_2)
            else:
                self.__triangulate_polygon_recursive(polygon_2,trianglesList)
        return trianglesList

    def __triangleArea(self, tr:list)->float:
        if len(tr)<3:
            return 0
        a = 0.5 * ( tr[0].x * (tr[1].y - tr[2].y) + tr[1].x * (tr[2].y - tr[0].y) + tr[2].x * (tr[0].y - tr[1].y))
        
        return a


    def triangulate_polygon(self)->list:
        '''Search triangles and return the triangle list'''
        if len(self.trianglesList) == 0:
            self.trianglesList = self.__triangulate_polygon_recursive(self.vertices,self.trianglesList)
        return self.trianglesList
    
    def __compute(self):
        '''Compute the area of the triangulated polygone'''

        if self.area == None:
            self.safeOrder()
            area = 0.0
            Mx = 0.0
            My = 0.0
            self.triangulate_polygon()
            for triangle in self.trianglesList:
                xcog_tr = 0.0
                ycog_tr = 0.0
                area_tr = 0.0
                for pt in triangle:
                    xcog_tr += pt.x
                    ycog_tr += pt.y
                
                area_tr = self.__triangleArea(triangle)
                Mx += xcog_tr / 3.0 * area_tr
                My += ycog_tr / 3.0 * area_tr

                area += area_tr

            self.area = abs(area)
            self.cog = Point2D(Mx/area,My/area)

    '''******************* Ordering ******************************************'''

    def __checkPolyOrientation(self)->int:
        '''Return +1 for counterclockwise -1 for clockwise'''
        
        '''find the lowest point at left side to check the angle'''
        shortList = [0]
        n = len(self.vertices)
        x = self.vertices[0].x
        
        index = 0
        for i in range(1,n):
            if math.isclose(self.vertices[i].x,x):
                shortList.append(i)
            elif self.vertices[i].x < x:
                x = self.vertices[i].x
                shortList =[]
                index = i
        
        if len(shortList):
            y = self.vertices[shortList[0]].y
            for i in shortList:
                if self.vertices[i].y < y:
                    y = self.vertices[i].y
                    index = i
        
        A = self.vertices[nextVertex(n,index,-1)]
        B = self.vertices[index]
        C = self.vertices[nextVertex(n,index,+1)]

        return sign((B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y))
    
    def safeOrder(self):
        '''Arrange the list counterclockwise for safe triangulation operations'''
        if self.isCounterClockwise:
            return
        
        if self.__checkPolyOrientation() == -1:
            self.vertices.reverse()
        
        self.isCounterClockwise

    
    '''******************* Intersection s******************************************'''
   
    def __SplitEdges(self,line:list):
        self.splitPoly =[]
        self.intersections = []
        n = len(self.vertices)
        o = 0
        for i in range(n):
            startPt = self.vertices[i]
            endPt = self.vertices[nextVertex(n,i,1)]
            #(-) right side (+) left side (0) from the line
            edgeStartSide = sign(MtoSegment(line[0],line[1],startPt))
            edgeEndSide = sign(MtoSegment(line[0],line[1],endPt))

            self.splitPoly.append({'v':startPt, 'side':edgeStartSide, 'dist':0.0, 'next':0,'prev':0,'x':0})

            if edgeStartSide == 0:
                #Vertex is on the line
                self.intersections.append(i+o)
            else:
                if edgeStartSide * edgeEndSide < 0.0:
                    #edge crossing the split line
                    interPt = computeInter([startPt,endPt],line)
                    #TODO: check is is None ?
                    o += 1

                    self.splitPoly.append({'v':interPt, 'side':0, 'dist':0.0, 'next':0,'prev':0,'x':0})
                    self.intersections.append(i+o)
        
        #Connect the double linked list
        n = len(self.splitPoly)
        for i in range(n):
            self.splitPoly[i]['next'] = nextVertex(n,i,+1)
            self.splitPoly[i]['prev'] = nextVertex(n,i,-1)
                    
        return 
    
    def __SortIntersections(self,line:list):
        '''Sort the intersection list using distance to origine and fill the dist field'''

        def distanceFromFirst(i):
            '''return a signed distance from the start point of the line'''
            pt1 = line[0]
            pt2 = line[1]
            ptM = self.splitPoly[i]['v']
            return (ptM.x - pt1.x)*(pt2.x-pt1.x) + (ptM.y-pt1.y)*(pt2.y-pt1.y)
        
        self.intersections.sort(key=distanceFromFirst)

        for i in  self.intersections:
            self.splitPoly[i]['dist'] = distanceFromFirst(i)

        return

    def __createEdges(self,srcPt,destPt):
        '''create new corresponding edge'''
        self.splitPoly.append(dict(self.splitPoly[srcPt]))
        a=len(self.splitPoly) - 1
        self.splitPoly.append(dict(self.splitPoly[destPt]))
        b=len(self.splitPoly) - 1

        self.splitPoly[a]['next'] = destPt
        self.splitPoly[a]['prev'] = self.splitPoly[srcPt]['prev']
        self.splitPoly[b]['next'] = srcPt
        self.splitPoly[b]['prev'] = self.splitPoly[destPt]['prev']

        self.splitPoly[self.splitPoly[srcPt]['prev']]['next'] = a
        self.splitPoly[srcPt]['prev'] = b

        self.splitPoly[self.splitPoly[destPt]['prev']]['next'] = b
        self.splitPoly[destPt]['prev'] = a
            
    def __SplitPolygon(self):
        '''Split polygon creating new edges'''
        self.edges = []
        useSrc = None
        l = len(self.intersections)
        i = 0
        while i < l:

            #find source point
            srcPt = useSrc
            useSrc = None

            while (i < l) and (srcPt is None):
                curIndex = self.intersections[i]
                #curSide = splitPoly[curIndex]['side'] #Shall be 0
                prevSide = self.splitPoly[self.splitPoly[curIndex]['prev']]['side']
                nextSide = self.splitPoly[self.splitPoly[curIndex]['next']]['side']

                curDist = self.splitPoly[curIndex]['dist'] 
                prevDist = self.splitPoly[self.splitPoly[curIndex]['prev']]['dist']
                nextDist = self.splitPoly[self.splitPoly[curIndex]['next']]['dist']

                
                if ((prevSide == +1 and nextSide == -1) or 
                    (prevSide == +1 and nextSide == 0 and nextDist < curDist)  or 
                    (prevSide == 0 and nextSide == -1 and prevDist < curDist)):
                    srcPt = self.intersections[i]
                else:
                    i+=1
                #i+=1
            
            #find destination point
            destPt = None
            while (i < l) and (destPt is None):
                curIndex = self.intersections[i]
                #curSide = splitPoly[curIndex]['side'] #Shall be 0
                prevSide = self.splitPoly[self.splitPoly[curIndex]['prev']]['side']
                nextSide = self.splitPoly[self.splitPoly[curIndex]['next']]['side']

                if ((prevSide == -1 and nextSide == +1) or 
                    (prevSide == 0 and nextSide == +1)  or 
                    (prevSide == -1 and nextSide == 0 ) or
                    (prevSide == -1 and nextSide == -1 ) or
                    (prevSide == +1 and nextSide == +1 )):
                    destPt = self.intersections[i]
                else:
                    i+=1
                            
            if not((srcPt is None) or (destPt is None)):
                #Bridge source and destination
                self.edges.append([self.splitPoly[srcPt]['v'],self.splitPoly[destPt]['v']])
                self.__createEdges(srcPt,destPt)


                #Is it a config in which a vertex needs to be reused as source vertex ?
                if self.splitPoly[self.splitPoly[self.splitPoly[srcPt]['prev']]['prev']]['side'] == +1 :
                    useSrc = self.splitPoly[srcPt]['prev']
                elif self.splitPoly[self.splitPoly[destPt]['next']]['side'] == -1 :
                    useSrc = destPt
            i+=1
    
        return

    def __CollectPolys(self):
        resPolys =[]

        for i,e in enumerate(self.splitPoly):
            if not e['x']:
                poly = []
                curEdge = e
                finished = False
                while not finished:
                    curEdge['x'] = 1
                    poly.append(curEdge['v'])
                    if curEdge['next'] == i:
                        finished = True
                    curEdge = self.splitPoly[curEdge['next']]
                
                resPolys.append(poly)
                
        return resPolys

    def splitPolygon(self,line:list):
        self.safeOrder()
        self.__SplitEdges(line)
        self.__SortIntersections(line)
        self.__SplitPolygon()
        return {'pList':self.__CollectPolys(),'edges':self.edges}

    def getPolygonsFromSide(self,line:list,side:int):
        '''return a list of poly cutted with the line on the selected side (-1:right or +1:left) of the line'''
        polyList = self.splitPolygon(line)['pList']
        polyFromSide = []

        #check if the polygon is on right side
        for p in polyList:
            i = 0
            l = len(p)
            finished = False
            while (i < l) and not finished:
                s = sign(MtoSegment(line[0],line[1],p[i]))
                if s !=0: #is side == 0 the vertex is on the line
                    finished = True
                    if s == side:
                        polyFromSide.append(p)
                else: # s == 0 the vertex is on the line
                    i+=1

        return {'pList':polyFromSide,'edges':self.edges}


    def getArea(self)->float:
        '''Compute the area of the triangulated polygone'''
        self.__compute()
        return self.area

    def getCog(self)->Point2D:
        '''Compute the area of the triangulated polygone'''
        self.__compute()
        return self.cog
    
    def getEdges(self)->list:
        '''Compute the area of the triangulated polygone'''
        self.__compute()
        return self.edges


class Polygons:
    '''Class to deal with multiples polygons'''
    def __init__(self, cuttedPoly:dict = None, pList:list = None):
            self.p_list = []
            self.edges = []
            self.area = None
            self.cog = None
            if not (pList is None):
                for p in pList:
                    self.p_list.append(Polygon(p))
            elif not (cuttedPoly is None):
                for p in cuttedPoly['pList']:
                    self.p_list.append(Polygon(p))
                self.edges = list(cuttedPoly['edges'])
           

    def __str__(self):
        strPoly = 'Polygons['
        first = True
        for p in self.p_list:
            if first:
                first = False
            else:
               strPoly +=', '     
            strPoly += str(p)
        strPoly +=']' 
        return strPoly

    def __repr__(self):
        return str(self)

    def __compute(self):
        if self.area == None:
            self.area = 0.0
            Mx=0.0
            My=0.0
            for p in self.p_list:
                self.area += p.getArea()
                Mx += p.cog.x * p.getArea()
                My += p.cog.y * p.getArea()
            if self.area:
                self.cog = Point2D(Mx /self.area, My /self.area)
            else:
                self.cog = Point2D(0.0, 0.0)

    
    def getArea(self)->float:
        '''Compute the area of the triangulated polygone'''
        self.__compute()
        return self.area
    
    def getCog(self)->Point2D:
        '''Compute the cog of the triangulated polygone'''
        self.__compute()
        return self.cog
    
    def getEdges(self)->list:
        '''Return the edges from the previous cut'''
        return self.edges
