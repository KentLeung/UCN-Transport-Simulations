import os
### READ REGIONS FROM TEXT FILE... RETURN REGIONS AS STRUCTURE

#############################################################
# Current method reads everybody into list of strings, I think.
# Manipulating immutable strings is stupid expensive...
# Consider rewriting this with a more memory efficient method/structure.
#############################################################
# ARR Jan 2012
def read_regions(regionsfile = 'Regionfile'):
  # READ REGIONS FILE INTO 2-D LIST
  regionsfile = open(regionsfile,'r')
  rawtext = regionsfile.read() #extract text
  
  regionsfile.close()
  
  rawtext = rawtext.replace('*','-1,-1,-1') #allows map() to be called later.
  rawtext = rawtext.replace(' ',',') #Convert into pseudo-csv file
  for i in range(0,7):
    rawtext = rawtext.replace(',,',',') #Make csv nicely delimited
  lineregions = rawtext.splitlines()
  elements = len(lineregions)-2
  stringregions = []
  
  for i in range(1,elements+1):
    stringregions.append(lineregions[i].split(','));

  #Convert read text values to floats
  tempreg = [map(float,x) for x in stringregions]
  regions = []
  for region in tempreg:
    region[0]=int(region[0])
    region[1]=int(region[1])
    for counter in range(-5,0):
      region[counter] = int(region[counter])
    regions.append(region)

  return regions


################################################################################
### READ CONNECTIONS FILE INTO 2-D LIST- return text file
################################################################################
# This is not ideal, but it is presently the best way to handle the data that
# can be given by the connection file. As a matter of fact, this information is
# not presently used by 'render_geometry'. It may be in the future, so we may as
# well be in the habit of reading it now. It's pretty cheap.
################################################################################
# ARR Jan, 2012

def read_connex(connectionfile = 'Connexfile'):
  connectionsfile = open(connectionfile,'r');
  conntext = connectionsfile.read();
  
  connectionsfile.close();
  
  #some (but not all) elements are stored as lists. Preserve this information
  conntext = conntext.replace('\t','-');
  conntext = conntext.replace(' ','-');
  for i in range(0,7):
    conntext = conntext.replace('--','-'); #Format nicely, again.
  lineconn = conntext.splitlines();
  
  seams = len(lineconn)-2;
  
  stringconnex = [];
  
  for i in range(3,seams+2):
    stringconnex.append(lineconn[i].split('-'));
  
  #Connections are still formatted as strings. 
  #This is because some of their elements are lists.
  #Also note: last element in each row is '/'
  return stringconnex
  
def read_simfile(filename):
  #READ HITPLACES INTO MEMORY, DRAW PATHS
  hitsfile = open(filename,'r');
  rawtext = hitsfile.read(); #extract text
  
  hitsfile.close();
  
  linehits = rawtext.splitlines();
  stringhitplaces = [];
  
  for lines in linehits:
    stringhitplaces.append(lines.split(','));
  
  #Convert read text values to floats
  hitplaces = [map(float,x) for x in stringhitplaces]
  return hitplaces
