import numpy, math, scipy

# p == n -> Jacobi. Cholesky
# p > n  -> Descompuneri

# np.linalg.cholesky()

PRECISION = 6
norma = lambda x: numpy.linalg.norm(x, ord=1)

EPS = 10 ** -PRECISION


def read_data(file_name: str):
	with open(file_name, "rt") as f:
		data = f.read()

	data = list(map(lambda a: float(a), data.split()))

	p, n = int(data[0]), int(data[1])
	data = data[2:]

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




a = read_data("t5.in")
p, n = a.shape
if   p == n:
	jacobi(a)
	sir_ch(a)
elif p > n:
	svd(a)
