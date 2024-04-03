import numpy as np, math

EPS = 1e-10
# matrar [0] = val, [1] = col

norma = lambda x: np.linalg.norm(x, ord=math.inf)

#vector_rar(vector_rar)
class MatRara:
	def __init__(self, n):
		self.n = n
		self.il = [-1 for _ in range(n)] # inceput_linii
		self.vals = []

	def line_end(self, lin):
		l_end = len(self.vals)
		for l in range(lin+1, self.n):
			if self.il[l] != -1:
				l_end = self.il[l]
				break
		return l_end

	def add_elem(self, val, lin, col):
		if val == 0:
			return
		if self.il[lin] == -1:
			self.il[lin] = self.line_end(lin-1 if lin > 0 else 0)
		line = self.vals[self.il[lin]:self.line_end(lin)]

		idd = len(line)
		for idx, i in enumerate(line):
			if i[1] == col:
				i[0] += val
				if i[0] == 0:
					self.vals.pop(self.il[lin] + idx)
					for ll in range(lin+1, self.n):
						if self.il[ll] != -1:
							self.il[ll] -= 1
					print("Removed 0 elem after sum")
				return
			elif i[1] > col and (idd == len(line) or i[1] < line[idd][1]):
				idd = idx
		self.vals.insert(self.il[lin]+idd, [val, col])
		for ll in range(lin+1, self.n):
			if self.il[ll] != -1:
				self.il[ll] += 1

	def __str__(self):
		s = "-*-\n"
		s += f"inceput_linii:{self.il[:5]}..\n"
		s += f"vals:{self.vals[:5]}..\n"
		for l in range(5):
			s += self.il[l].__str__() + self.vals[self.il[l]:self.line_end(l)].__str__() + '\n'
		s += '...\n'
		for l in range(-1, -5, -1):
			s += self.il[l].__str__() + self.vals[self.il[l]:self.line_end(l)].__str__() + '\n'
		s += "-*-"
		return s

	def __eq__(self, a):
		print("Using custom eq")
		if a.n != self.n:
			return False
		for lid in range(self.n):
			if a.il[lid] != self.il[lid]:
				print("[EQ] Lines start differently")
				return False
			else:
				lstart = a.il[lid]
				if lstart == -1:
					continue
				lenda, lends = a.line_end(lid), self.line_end(lid)
				if lenda != lends:
					print("[EQ] Different number of elements on a line")
					return False
				for aaa, bbb in zip(a.vals[lstart:lenda], self.vals[lstart:lends]):
					if not (aaa[0] == bbb[0] and aaa[1] == bbb[1]):
						print(lid, aaa, bbb, self.vals[lstart:lends], a.vals[lstart:lenda], self.il[lid-1:lid+1], a.il[lid-1:lid+1], sep='\n')
						print("[EQ] Different elems on same position")
						return False
		return True

def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS

def read(file_name):
	print(f"Reading {file_name}")
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
	print(f"Read {file_name}")
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

def bonus_sum(a, b):
	if a.n != b.n:
		return

	c = MatRara(a.n)

	for l in range(a.n):
		if a.il[l] != -1:
			for co in a.vals[a.il[l]:a.line_end(l)]:
				c.add_elem(co[0], l, co[1])
	for l in range(b.n):
		if b.il[l] != -1:
			for co in b.vals[b.il[l]:b.line_end(l)]:
				c.add_elem(co[0], l, co[1])
	return c

def get_diag(a):
	diag = np.zeros((a.n,))
	for l in range(a.n):
		for c in a.vals[a.il[l]:a.line_end(l)]:
			if c[1] == l:
				diag[l] = c[0]
				break
	return diag

def diag_nul(a):
	return not np.any(get_diag(a))

norma = lambda x: np.linalg.norm(x, ord=2)

def gauss_streidel(a, b):
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
			column = [(l[0] * xc[l[1]]) for l in a.vals[a.il[i]:a.line_end(i)] if l[1] != i]
			temp = (b[i] - sum(column)) / diag[i]
			dx = max(dx, xc[i] - temp)
			xc[i] = temp
		# print(xc)
		k += 1
	if is_zero(dx):
		return xc, k
	else:
		return 'divergenta'


for i in reversed(range(6)):
	a = read(f"tema4/a_{i}.txt")
	b = read_b(f"tema4/b_{i}.txt")
	rez = gauss_streidel(a, b)
	if isinstance(rez, str):
		print(i, rez)
	else:
		x, steps = rez
		print(i, x, steps, "pasi")
		verif = np.zeros_like(b)
		for l in range(a.n):
			verif[l] = sum([(c[0] * x[c[1]]) for c in a.vals[a.il[l]:a.line_end(l)]])
		print("norma:", norma(verif - b))

# bonus:
a = read("tema4/a.txt")
b = read("tema4/b.txt")
aplusb = read("tema4/aplusb.txt")
c = bonus_sum(a, b)
print(c)
print(f"{c == aplusb}")