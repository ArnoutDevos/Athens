#!/usr/bin/env python
##############################################################################
# EVOLIFE  www.dessalles.fr/Evolife                    Jean-Louis Dessalles  #
#            Telecom ParisTech  2014                       www.dessalles.fr  #
##############################################################################


##############################################################################
#  Curves                                                                    #
##############################################################################


""" EVOLIFE: Module Curves:
	Stores data that can be used to plot curves or to store into a file
"""


import sys
if __name__ == '__main__':  sys.path.append('../..')  # for tests

from Evolife.Tools.Tools import transpose, error


##################################################
# Evolife colours                                #
##################################################
EvolifeColours = ['#808080', 'black', 'white', 'blue', 'red', 'yellow', '#A06000', '#0080A0', '#FF80A0', '#94DCDC', 
			'#008000', '#009500', '#00AA00', '#00BF00', '#00D400', '#00E900', '#00FE00', '#64FF64', '#78FF78', '#8CFF8C', '#A0FFA0', '#B4FFB4',
			'#800000', '#950000', '#AA0000', '#BF0000', '#D40000', '#E90000', '#FE0000', '#FF6464', '#FF7878', '#FF8C8C', '#FFA0A0', '#FFB4B4',
			'#000080', '#000095', '#0000AA', '#0000BF', '#0000D4', '#0000E9', '#0000FE', '#6464FF', '#7878FF', '#8C8CFF', '#A0A0FF', '#B4B4FF',
		   ]
EvolifeColourNames = ['grey', 'black', 'white', 'blue', 'red', 'yellow', 'brown', 'blue02', 'pink', 'lightblue', 
			'green', 'green1', 'green2', 'green3', 'green4', 'green5', 'green6', 'green7', 'green8', 'green9', 'green10', 'green11',  # 21
			'red0', 'red1', 'red2', 'red3', 'red4', 'red5', 'red6', 'red7', 'red8', 'red9', 'red10', 'red11', # 33
			'blue0', 'blue1', 'blue2', 'blue3', 'blue4', 'blue5', 'blue6', 'blue7', 'blue8', 'blue9', 'blue10', 'blue11',
			]

# Green colours between 10 and 21
# Red colours between 22 and 33
# Blue colours between 34 and 45
Shades = {'green':(10, 21), 'red':(22,33), 'blue':(34, 45)}

def Shade(x, BaseColour='green', Min=0, Max=1, darkToLight=True, invisible='white'):
	" compute a shade for a given base colour "
	if Min != Max and x >= Min and x <= Max:
		shades = Shades[BaseColour]
		if darkToLight:
			return EvolifeColourNames[shades[0] + int(((x - Min) * (shades[1] - shades[0])) / (Max - Min))]
		else:
			return EvolifeColourNames[shades[1] - int(((Max - x) * (shades[1] - shades[0])) / (Max - Min))]
	return invisible
	
def EvolifeColourID(Colour_designation, default=(4,'red')):
	ID = None
	try:
		if str(Colour_designation).isdigit() and int(Colour_designation) in range(len(EvolifeColours)):	# colour given by number
			ID = int(Colour_designation)
		elif Colour_designation in EvolifeColourNames:	# colour given by name
			ID = EvolifeColourNames.index(Colour_designation)
		elif Colour_designation.startswith('#') and Colour_designation in EvolifeColours:	# colour given by code
			ID = EvolifeColours.index(Colour_designation)
		if ID is not None:	return (ID, EvolifeColours[ID])
	except AttributeError:	pass
	return default


##################################################
# Stroke: drawing element (point or segement)    #
##################################################
class Stroke:
	"""	stores coordinates as:
		(x, y, colour, thickness)
	"""
	def __init__(self, Coordinates):
		DefCoord = 10	# default value
		DefColour = 4	# default value
		DefSize = 3		# default value
		DefaultStroke = (DefCoord, DefCoord, DefColour, DefSize)
		if Coordinates:
			self.Coord = Coordinates[:4] + DefaultStroke[min(len(Coordinates), 4):4] # completing with default values
			(self.x, self.y, self.colour, self.size) = self.Coord
		else:
			self.Coord = None
			
	def point(self):
		return (self.x, self.y)
	
	def scroll(self):
		self.y -= 1
		C1 = list(self.Coord)
		C1[1] -=1
		self.Coord = tuple(C1)
	
	def __add__(self, Other):	# allows to add with None
		if Other.Coord:	return self.Coord + Other.Coord
		else:		return self.Coord

	def __repr__(self):	return '%s, %s, %s, %s' % (str(self.x), str(self.y), str(self.colour), str(self.size))
	
##################################################
# Curve: stores points to display a curve        #
##################################################
class Curve(object):
	""" Holds a complete (continuous) curve in memory
	"""
	def __init__(self, colour, ID, ColName=None, Legend=None):
		self.ID = ID	# typically: number of the curve
		self.colour = colour	# Evolife colour
		self.Name = str(ID)	# default, but will receive a string
		try:	self.ColName = EvolifeColourNames[ID]	# default
		except IndexError:	self.ColName = colour
		if ColName is not None:  
			self.ColName = ColName
			self.Name = ColName	# second default
		self.Legend = Legend if Legend is not None else self.Name
		self.erase()

	def erase(self):
		self.start((0,0))

	def start(self,StartPos):
		self.CurrentPosition = 0   # Current position for reading
		self.positions = [StartPos] # Stores successive points
		self.discontinuities = []

	def name(self, N = ""):
		if N != "":	self.Name = N
		return self.Name

	def legend(self, L=""):
		if L:	self.Legend = L
		return self.Legend
		
	def last(self):
		return self.positions[-1]

	def add(self, Pos, Draw=True):
		if not Draw:
			self.discontinuities.append(self.length())
		self.positions.append(Pos)

	def length(self):
		return len(self.positions)
	
	def X_coord(self):
		" list of x-coordinates "
		return tuple(map(lambda P: P[0], self.positions))
		
	def Y_coord(self):
		" list of y-coordinates "
		return tuple(map(lambda P: P[1], self.positions))

	def Avg(self, start=0):
		" compute average value of Y_coord "
		#ValidValues = [Y for Y in self.Y_coord()[start:] if Y >= 0]
		ValidValues = [P[1] for P in self.positions if P[0] >= start and P[1] >= 0]
		if len(ValidValues):
			return int(round(float(sum(ValidValues)) / len(ValidValues)))
		else:
			return 0

	def __iter__(self):
		# defines the class as an iterator
		return self

	def __next__(self):
		"""2.6-3.x version"""
		return self.next()

	def next(self):
		" Iteratively returns segments of the curve "
		if self.CurrentPosition+1 in self.discontinuities:
			# one segment must be skipped
			self.CurrentPosition += 1
		if self.length() < 2 or self.CurrentPosition >= self.length()-1:
			self.CurrentPosition = 0	# ready for later use
			raise StopIteration
		self.CurrentPosition += 1
		return (self.positions[self.CurrentPosition-1], self.positions[self.CurrentPosition])
			
	def __repr__(self):
		return self.Name

##################################################
# Curves: list of curves                         #
##################################################
class Curves(object):
	""" Stores a list of 'Curves'
	"""
	
	def __init__(self):
		self.Colours = EvolifeColours
		self.Curves = [Curve(Colour, Number, EvolifeColourNames[Number]) for (Number,Colour) in enumerate(EvolifeColours)]


	def start_Curve(self, Curve_id, location):
		""" defines where a curve should start
		"""
		try:
			self.Curves[Curve_id].start(location)
		except IndexError:
			error("Curves: unknown Curve ID")
			

	def Curvenames(self, Names):
		""" records names for Curves
		"""
		Str =  '\nDisplay: \n\t'
		try:
			for Curve_description in Names:
				(Curve_designation, Name, Legend) = Curve_description + (0, '', '')[len(Curve_description):]
				for P in self.Curves:
					CurveId = EvolifeColourID(Curve_designation, default=None)[0]
					if P.ID == CurveId:
						P.name(Name)
						P.legend(Legend if Legend else Name)
						Str += '\n\t%s:\t%s' % (P.ColName, P.legend())
			Str += '\n'
		except IndexError:
			error("Curves: unknown Curve ID")
		return Str
		
	def Legend(self):
		""" returns legend string representing active curves
		"""
		return [(P.ID, P.colour, P.ColName, P.Name, P.Legend) for P in self.Curves if len(P.positions) > 1]
	
	def dump(self, ResultFileName=None, ResultHeader='', DumpStart=0):
		""" saves Curves to a file
		"""
		if ResultFileName == None:
			return
		# DumpStart = points below this x-value are removed from the computation of average values

		# dump: classification of Curves sharing x-coordinates
		X_coordinates = list(set([P.X_coord() for P in self.Curves]))
		X_coordinates.sort(key=lambda x: len(x), reverse=True)
		if len(X_coordinates) <= 2:
			# only one Curve or several Curves sharing x-coordinates
			active_Curves = [P for P in self.Curves if P.X_coord() == X_coordinates[0]]
			# print 'saving Curves %s to %s' %  (active_Curves,FileName)
			Coords = [('Year',) + tuple([P.name() for P in active_Curves])]
			Coords += transpose([X_coordinates[0]] \
								 + [P.Y_coord() for P in active_Curves])
		else:
			active_Curves = self.Curves
			Coords = reduce(lambda x,y: x+y, [P.positions for P in self.Curves
											  if len(P.positions) > 1]) 
			
		File_dump = open(ResultFileName + '.csv', 'w')
		for C in Coords:
			File_dump.write(';'.join([str(x) for x in C]))
			File_dump.write('\n')
		File_dump.close()

		# editing the header
		HeaderLines = ResultHeader.split('\n')
		HeaderLines[0] += 'LastStep;'
		# Writing Curve names sorted by colours at the end of the first line
		HeaderLines[0] += ';'.join([P.name() for P in active_Curves])
		Header = '\n'.join(HeaderLines)

		# storing average values
		Averages = open(ResultFileName + '_res.csv', 'w')
		Averages.write(Header)
		Averages.write(str(active_Curves[0].X_coord()[-1]) +';')	# storing actual max time value
		Averages.write(';'.join([str(P.Avg(DumpStart)) for P in active_Curves]))
		Averages.write('\n')
		Averages.close()
		return
			


if __name__ == "__main__":

	print(__doc__)

__author__ = 'Dessalles'
