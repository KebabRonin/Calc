import numpy as np, math

EPS = 1e-10
# matrar [0] = val, [1] = col

norma = lambda x: np.linalg.norm(x, ord=math.inf)

class MatRara:
	def __init__(self, n):
		self.n = n
		self.m = [[] for _ in range(n)]

	def add_elem(self, val, lin, col):
		if val == 0:
			return
		line = self.m[lin]
		idd = len(line)
		for idx, i in enumerate(line):
			if i[1] == col:
				i[0] += val
				if i[0] == 0:
					line.remove(i)
					print("Removed 0 elem after sum")
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

	def __eq__(self, a):
		print("Using custom eq")
		if a.n != self.n:
			return False
		for i in range(self.n):
			if len(a.m[i]) != len(self.m[i]):
				print("Equality: line lenghts not equal")
				print(len(a.m[i]), len(self.m[i]))
				print(a.m[i])
				print(self.m[i])
				return False
			for aa, bb in zip(a.m[i], self.m[i]):
				if not (aa[0] == bb[0] and aa[1] == bb[1]):
					return False
		return True

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
	if a.n != b.n:
		return

	c = MatRara(a.n)

	for i in range(a.n):
		amax, bmax = len(a.m[i]), len(b.m[i])
		ida = 0
		idb = 0
		while ida < amax and idb < bmax:
			if a.m[i][ida][1] == b.m[i][idb][1]:
				c.add_elem(a.m[i][ida][0], i, a.m[i][ida][1])
				c.add_elem(b.m[i][idb][0], i, b.m[i][idb][1]) # should add them together
				ida += 1
				idb += 1
			elif a.m[i][ida][1] < b.m[i][idb][1]:
				c.add_elem(a.m[i][ida][0], i, a.m[i][ida][1])
				ida += 1
			else:
				c.add_elem(b.m[i][idb][0], i, b.m[i][idb][1])
				idb += 1
		while ida < amax:
			c.add_elem(a.m[i][ida][0], i, a.m[i][ida][1])
			ida += 1
		while idb < bmax:
			c.add_elem(b.m[i][idb][0], i, b.m[i][idb][1])
			idb += 1
	return c

def get_diag(a):
	diag = np.zeros((a.n,))
	for idx, l in enumerate(a.m):
		for c in l:
			if c[1] == idx:
				diag[idx] = c[0]
				break
		else:
			diag[idx] = 0
	return diag

def diag_nul(a):
	return not np.any(get_diag(a))

norma = lambda x: np.linalg.norm(x, ord=2)

def gauss_streidel_MatRara(a, b):
	xc = np.arange(a.n, dtype=float) + 1
	kmax = 10_000
	dx = 100
	k = 0
	diag = get_diag(a)
	if diag_nul(a):
		return "No solution, diag is null"
	while not is_zero(dx) and k <= kmax and dx <= 1e8:
		dx = 0
		for i in range(a.n):
			column = [(l[0] * xc[l[1]]) for l in a.m[i] if l[1] != i]
			temp = (b[i] - sum(column)) / diag[i]
			dx = max(dx, xc[i] - temp)
			xc[i] = temp
		k += 1
	if is_zero(dx):
		return xc, k
	else:
		return 'divergenta'


for i in range(6):
	a = read_MatRara(f"tema4/a_{i}.txt")
	b = read_b(f"tema4/b_{i}.txt")
	rez = gauss_streidel_MatRara(a, b)
	if isinstance(rez, str):
		print(i, rez)
	else:
		x, steps = rez
		print(i, x, steps, "pasi")
		verif = np.zeros_like(b)
		for i in range(a.n):
			verif[i] = sum([l[0] * x[l[1]] for l in a.m[i]])
		print("norma:", norma(verif - b))

# bonus:
# a = read_MatRara("tema4/a.txt")
# b = read_MatRara("tema4/b.txt")
# aplusb = read_MatRara("tema4/aplusb.txt")
# c = bonus_sum_MatRara(a, b)
# print(c)
# print(f"{c == aplusb}")