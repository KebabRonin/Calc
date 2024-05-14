import numpy as np, math

EPS = 1e-10
# matrar [0] = val, [1] = col

norma = lambda x: np.linalg.norm(x, ord=math.inf)

#vector_rar(vector_rar)
class MatRara:
	def __init__(self, n):
		self.n = n
		self.m = dict()

	def add_elem(self, val, lin, col):
		if val == 0:
			return
		if lin not in self.m:
			self.m[lin] = []
		line = self.m[lin]
		# for idx, l in enumerate(self.m):
		# 	if l[1] == lin:
		# 		line = l[0]
		# 		break
		# 	elif l[1] > lin and (idd == len(self.m) or l[1] < self.m[idd][1]):
		# 		idd = idx
		# else:
		# 	# didnt break
		# 	self.m.insert(idd, [[], lin])
		# 	line = self.m[idd][0]

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
		keis = sorted(list(self.m.keys()))
		for l in keis[:5]:
			s += self.m[l].__str__() + '\n'
		s += '...\n'
		for l in keis[-5:]:
			s += self.m[l].__str__() + '\n'
		s += "-*-"
		return s

	def __eq__(self, a):
		print("Using custom eq")
		if a.n != self.n:
			return False
		if len(a.m) != len(self.m):
			print("Equality: nr non-zero lines not equal")
			return False
		for aa in a.m:
			if aa not in self.m:
				print(aa, 1)
				return False
			else:
				for aaa, bbb in zip(a.m[aa], self.m[aa]):
					if not (aaa[0] == bbb[0] and aaa[1] == bbb[1]):
						print(aa, 2, aaa, bbb)
						return False
		for bb in self.m:
			if bb not in a.m:
				print(bb, 1)
				return False
			else:
				for aaa, bbb in zip(a.m[bb], self.m[bb]):
					if not (aaa[0] == bbb[0] and aaa[1] == bbb[1]):
						print(bb, 1, aaa, bbb)
						return False
		return True

def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS

def read(file_name):
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

def bonus_sum(a, b):
	if a.n != b.n:
		return

	c = MatRara(a.n)

	for l in a.m:
		for co in a.m[l]:
			c.add_elem(co[0], l, co[1])
	for l in b.m:
		for co in b.m[l]:
			c.add_elem(co[0], l, co[1])
	return c

def get_diag(a):
	diag = np.zeros((a.n,))
	for l in range(a.n):
		if l in a.m:
			for c in a.m[l]:
				if c[1] == l:
					diag[l] = c[0]
					break
		else:
			diag[l] = 0
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
			if i not in a.m:
				column = 0
			else:
				column = [(l[0] * xc[l[1]]) for l in a.m[i] if l[1] != i]
			temp = (b[i] - sum(column)) / diag[i]
			dx = max(dx, xc[i] - temp)
			xc[i] = temp
		# print(xc)
		k += 1
	if is_zero(dx):
		return xc, k
	else:
		return 'divergenta'


for i in range(6):
	a = read(f"tema4/a_{i}.txt")
	b = read_b(f"tema4/b_{i}.txt")
	rez = gauss_streidel(a, b)
	if isinstance(rez, str):
		print(i, rez)
	else:
		x, steps = rez
		print(i, x, steps, "pasi")
		verif = np.zeros_like(b)
		for l in a.m:
			verif[l] = sum([(c[0] * x[c[1]]) for c in a.m[l]])
		print("norma:", norma(verif - b))

# bonus:
a = read("tema4/a.txt")
b = read("tema4/b.txt")
aplusb = read("tema4/aplusb.txt")
c = bonus_sum(a, b)
print(c)
print(f"{c == aplusb}")