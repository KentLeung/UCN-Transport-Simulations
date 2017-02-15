from __future__ import division
from visual import *

### Draw imported geometry using VPython objects
def draw_geom(regions):
  # Special region codes are not currently rendered.
  # If you need this, or know of a pretty way to represent it, do let me know.
  pieces = [];
  patches = [];
  endplate = vector(0,0,0);
  direction = vector(0,0,0);
  for i in regions:
    if i[1]==1: #Box-type region.
      #First, determine orientation of box by rotating unit vector.
      direction = rotate(vector(0,0,1),radians(i[8]),vector(0,1,0));
      xaxis = rotate(vector(1,0,0), radians(i[8]), axis = vector(0,1,0));
      direction = direction.rotate(angle = radians(i[9]), axis = xaxis);
      #Next, determine location of box. Two options here:
      #[0] basepoint is -1,-1,-1: object is centered on previous piece
      #[1] basepoint is anything else: This three-vector is the new basepoint.
      if not (i[2]==-1 and i[3]==-1 and i[4]==-1):
        endplate = vector(i[2],i[3],i[4])
      piece = box(pos = endplate+direction*i[5]/2, axis = direction, 
                  length = i[5], width = i[6], height = i[7], 
                  color = color.white, opacity = .2);
      endplate = piece.pos+piece.axis*piece.length/2;
    elif i[1]==2: #Cylinder
      direction = rotate(vector(0,0,1),radians(i[7]),vector(0,1,0));
      xaxis = rotate(vector(1,0,0), radians(i[7]), axis = vector(0,1,0));
      direction = direction.rotate(angle = radians(i[8]), axis = xaxis);
      if not (i[2]==-1.0 and i[3]==-1.0 and i[4]==-1.0):
        endplate = vector(i[2],i[3],i[4])
      piece = cylinder(pos=endplate, axis = direction*i[6], radius = i[5]/2,
                       color = color.magenta,opacity=.4);
      # color-code regions by their Fermi potential
      if i[19]<10:
        piece.color = color.white;
      elif i[19]<200:
        piece.color = color.cyan;
      else:
        piece.color = color.orange;
      endplate = endplate+piece.axis;
    if i[0]==0:
      piece.color = color.white;
      piece.material = materials.emissive;
    pieces.append(piece);
  #print elements, ' defined, and ',len(pieces),' rendered.'
  return pieces
  
### Draw all paths from imported text file.
def draw_simfile(hitplaces, fraction=1):
  total = hitplaces[-1][0]+1
  print 'The total number of neutrons is ', int(total)
  thisoften = round(total/(fraction*total))
  print 'Drawing every ',str(int(thisoften)),' neutrons.'
  counter = -1
  trails = []
  for point in hitplaces:
    if (counter+1)%thisoften==0:
      if point[0]!=counter: #new neutron- start a new trail
        trails.append(curve(color = color.red))
        counter = point[0]
      trails[-1].append(pos=vector(point[1],point[2],point[3]))
    else:
      counter = point[0]
  return trails

def draw_errorpoints():
    file = open('errorpoints.txt','r')
    filetxt = file.read()
    file.close()
            
    errors = filetxt.splitlines()
            
    errorpointlocations = []
    for line in errors:
        errorpointlocations.append(line.split(','))
    
    errorpointlocations = [map(float,x) for x in errorpointlocations]

    errorpoints = []
            
    for point in errorpointlocations:
        errorpoints.append(points(pos=(point[0],point[1],point[2]), size=5, color=color.green))
                           
    return errorpoints

