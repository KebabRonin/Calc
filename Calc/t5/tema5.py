import numpy, math, scipy

# p == n -> Jacobi. Cholesky
# p > n  -> Descompuneri

# np.linalg.cholesky()

PRECISION = 8
norma = lambda x: numpy.linalg.norm(x, ord=1)

EPS = 10 ** -PRECISION


def read_data(file_name: str):
	with open(file_name, "rt") as f:
		data = f.read()

	data = list(map(lambda a: float(a), data.split()))

	p, n = int(data[0]), int(data[1])
	data = data[2:]
	data = data[:n*p]

	if len(data) != n * p:
		raise Exception("Invalid data format. Format: first number is n(matrix dimensions), next n * n are matrix values, next n are s(vector) values")

	A_init = numpy.zeros(shape=(p, n), dtype=float)
	for i, v in enumerate(data):
		A_init[i//n, i%n] = v
	return A_init


def rand_data(n, p):
	A_init = numpy.random.rand(n, p)
	return A_init


def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS

def svd(a):
	print("SVDVals:", scipy.linalg.svdvals(a))
	print("Rank:", numpy.linalg.matrix_rank(a))
	print("Cond:", numpy.linalg.cond(a))
	mp = numpy.linalg.pinv(a)
	ls = (numpy.transpose(a) @ a) ** -1 @ numpy.transpose(a)
	print("Moore-Penrose:", mp)
	print("Least squares:", ls)
	print("Norma:", norma(mp - ls))

def getpq(a):
	el = -1
	pb, qb = None, None
	for p in range(a.shape[0]):
		for q in range(p):
			if math.fabs(a[p, q]) > el:
				el = math.fabs(a[p, q])
				pb = p
				qb = q
	# print(a, a[pb, qb], "\n")
	return pb, qb

def sign(a):
	return 1 if a >= 0 else -1

def getrot(a, p, q):
	alpha = (a[p, p] - a[q, q]) / (2*a[p, q])
	t = -alpha + sign(alpha) * math.sqrt(alpha**2 + 1)
	c = 1 / math.sqrt(1 + t**2)
	s = t / math.sqrt(1 + t**2)
	return c, s, t

def jacobi(a1):
	a = numpy.copy(a1)
	k, kmax = 0, 1_000
	u = numpy.identity(a.shape[0], dtype=float)
	while k <= kmax:
		p, q = getpq(a)
		if is_zero(a[p, q]):
			break
		c, s, t = getrot(a, p, q)
		# update a
		for j in range(a.shape[0]):
			if j != p and j != q:
				a[p, j]= c*a[p, j]+s*a[q, j]
				a[q, j]=a[j, q]=-s*a[j, p]+c*a[q, j]
		for j in range(a.shape[0]):
			if j != p and j != q:
				a[j, p] = a[p, j]
		a[p, p] += t * a[p, q]
		a[q, q] -= t * a[p, q]
		a[p, q] = a[q, p] = 0
		# update u
		for i in range(u.shape[0]):
			uv = u[i, p]
			u[i, p] = c * u[i, p] + s * u[i, q]
			u[i, q] = -s*uv+c*u[i, q]
		k += 1

	print("norma:", norma(a1@u-u@a))
	print("Afinal:\n", a)
	print("U:\n", u)

def sir_ch(a1):
	ainit=a1
	while True:
		l = numpy.linalg.cholesky(ainit)
		print(l)
		ainit = l @ numpy.transpose(l)
		a = numpy.transpose(l) @ l
		# print(a, ainit, a-ainit)
		if norma(a - ainit) < EPS:
			break
	return a


a = read_data("t5.in")
p, n = a.shape
if   p == n:
	jacobi(a)
	try:
		print(sir_ch(a))
	except numpy.linalg.LinAlgError as e:
		print("Cholesky can't be applied:", e)
elif p > n:
	svd(a)
