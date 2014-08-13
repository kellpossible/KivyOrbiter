from math import *

class Vector3(object):
	@staticmethod
	def ZERO():
		return Vector3(0,0,0)
	
	@staticmethod
	def X():
		return Vector3(1,0,0)
	@staticmethod
	def Y():
		return Vector3(0,1,0)
	@staticmethod
	def Z():
		return Vector3(0,0,1)

	@staticmethod
	def fromVector(other):
		return Vector3(other.x, other.y, other.z)
	
	@staticmethod
	def fromArray(ary):
		return Vector3(ary[0], ary[1], ary[2])
	def __init__(self, x=0.0, y=0.0, z=0.0):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.array = [self.x, self.y, self.z]

	def add(self, other):
		v = Vector3(0, 0, 0)
		v.x = self.x + other.x
		v.y = self.y + other.y
		v.z = self.z + other.z
		return v
		
	def sub(self, other):
		v = Vector3(0, 0, 0)
		v.x = self.x - other.x
		v.y = self.y - other.y
		v.z = self.z - other.z
		return v
	
	def copy(self, other):
		self.x = other.x
		self.y = other.y
		self.z = other.z
		return self
	
	def neg(self):
		return Vector3(-self.x, -self.y, -self.x)

	def mag(self):
		return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
	
	def dot(self, other):
		#returns a scalar
		return self.x * other.x + self.y * other.y + self.z * other.z
	
	
	def fmul(self, other):
		#assuming other is an integer
		return Vector3(self.x * other, self.y * other, self.z * other)
		
	def fdiv(self, other):
		return Vector3(self.x/other, self.y/other, self.z/other)
		
	def norm(self):
		m = self.mag()
		v = self.fdiv(m)
		return v
	
	def unit(self):
		return self.norm()
		
	def cross(self, other):
		v = Vector3()
		v.x = self.y * other.z - self.z * other.y
		v.y = self.z * other.x - self.x * other.z
		v.z = self.x * other.y - self.y * other.x
		return v
	
	def angle(self, other):
		return acos(self.dot(other)/(self.mag()*other.mag()))
	
	def rotate(self, rotation):
		r_vec_ary = [0.0,0.0,0.0]
		axis = 0
		while axis < len(rotation.matrix):
			i=0
			while i < len(rotation.matrix):
				r = rotation.matrix[axis][i]
				s = self.array[i]
				r_vec_ary[axis] +=  r*s
				
				#print("r: {0} s: {1}".format(r, s))
				i+=1
			axis+=1
			
		return Vector3.fromArray(r_vec_ary)
	
	def distance(self, other):
		return abs(other.sub(self).mag())

	def to_tuple(self):
		return (self.x, self.y, self.z)
	

	def __str__(self):
		return "Vec: {0}, {1}, {2}".format(self.x, self.y, self.z)	

	def __repr__(self):
		return str(self)
	
	
	
class Rotation(object):
	@staticmethod
	def aroundVector(u, a):
		"""u is unit vector for axis
			a is angle"""
		cosa = cos(a)
		sina = sin(a)
		matrix = [	[cosa + (u.x**2)*(1-cosa), 		u.x*u.y*(1-cosa) - u.z*sina, 	u.x*u.z*(1-cosa) + u.y*sina],
						[u.y*u.x*(1-cosa) + u.z*sina, 	cosa + (u.y**2)*(1-cosa), 		u.y*u.z*(1-cosa) - u.x*sina],
						[u.z*u.x*(1-cosa) - u.y*sina, 	u.z*u.y*(1-cosa) + u.x*sina,	cosa + (u.z**2)*(1-cosa)]]
		
		
		return Rotation(matrix)
	
	@staticmethod
	def aroundX(a):
		cosa = cos(a)
		sina = sin(a)
		matrix = [	[1, 0, 0],
					[0, cosa, -sina],
					[0, sina, cosa]]
		return Rotation(matrix)
	
	@staticmethod
	def aroundY(a):
		cosa = cos(a)
		sina = sin(a)
		matrix = [	[cosa, 0, sina],
					[0, 1, 0],
					[-sina, 0, cosa]]
		return Rotation(matrix)
	
	@staticmethod
	def aroundZ(a):
		cosa = cos(a)
		sina = sin(a)
		matrix = [	[cosa, -sina, 0],
					[sina, cosa, 0],
					[0, 0, 1]]
		return Rotation(matrix)
	
	def __init__(self, matrix):
		self.matrix = matrix

if __name__ == "__main__":
	v = Vector3(2, 0, 0)
	print(v.norm())
	
