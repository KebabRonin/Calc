import numpy as np, math

EPS = 1e-6

class MatRara:
	def __init__(self, n):
		self.n = n
		self.m = [[] for _ in range(n)]

	def add_elem(self, val, lin, col):
		line = self.m[lin]
		idd = len(line)
		for idx, i in enumerate(line):
			if i[1] == col:
				i[0] += val
				return
			elif i[1] > col and (idd == len(line) or i[1] < line[idd][1]):
				idd = idx
		line.insert(idd, [val, col])

	def __str__(self):
		s = "-*-\n"
		for l in self.m[:5] + ["..."] +  self.m[-5:]:
			s += l.__str__() + '\n'
		s += "-*-"
		return s

def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS

def read_MatRara(file_name):
	mat = None
	with open(file_name, "rt") as f:
		ls = f.readlines()
		n = int(ls[0])
		mat = MatRara(n)
		for l in ls[1:]:
			# print(l)
			l = l.split(',')
			if len(l) == 3:
				val, lin, col = float(l[0]), int(l[1]), int(l[2])
				mat.add_elem(val, lin, col)
	return mat
def read_b(file_name):
	b = None
	with open(file_name, "rt") as f:
		ls = f.readlines()
		n = int(ls[0])
		b = np.zeros((n,))
		for idx in range(n):
			b[idx] = float(ls[idx+1])
	return b

def bonus_sum_MatRara(a, b):
	if len(a) != len(b):
		return

	for i in range(len(a)):
		raise NotImplemented

def diag_nul(a):
	for idx, l in enumerate(a.m):
		for c in l:
			if c[1] == idx and c[0] != 0:
				break
		else:
			print(idx, l)
			return False
	return True

norma = lambda x: np.linalg.norm(x, ord=2)

def gauss_streidel(a, b):
	xc = xp = 0
	kmax = 10_000
	dx = 100
	k = 0
	while not is_zero(dx) and k <= kmax and dx <= 1e8:
		xp = xc
		dx = math.fabs(xp - xc)
		k += 1
	if is_zero(dx):
		return xc
	else:
		return 'divergenta'

a_1 = read_MatRara("tema4/a_1.txt")
b_1 = read_b("tema4/b_1.txt")
print(a_1.n, len(b_1))
print(diag_nul(a_1))