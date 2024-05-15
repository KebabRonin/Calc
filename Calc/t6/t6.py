import numpy as np, json
import matplotlib.pyplot as plt

#

PRECISION = 8
norma = lambda x: abs(x)

EPS = 10 ** -PRECISION

# pentru un x dat, aprox f(x) folosind:
# * Newton progresiv si Aitken
# * cmm patrate, Horner

def read_data(file_name = None):
	global decl_f
	if file_name is None:
		f = decl_f
		x0, xn, n = float(input('x0=')), float(input('xn=')), int(input('n='))
		if x0 > xn:
			x0, xn = xn, x0
		h = (xn - x0)/(n-1)
		xs = [x0+i*h for i in range(n-1)] + [xn]
		fxs = [f(x) for x in xs]
		return xs, fxs, f

	data = json.load(open(file_name, "rt"))
	if data['method'] == 'given':
		return data['xs'], data['fxs'], lambda x: 0
	elif data['method'] == 'compute':
		exec(data['f'])
		x0, xn, n = data['x0'], data['xn'], data['n']
		if x0 > xn:
			x0, xn = xn, x0
		h = (xn - x0)/(n-1)
		xs = [x0+i*h for i in range(n-1)] + [xn]
		fxs = [f(x) for x in xs]
		return xs, fxs, f
	else:
		raise Exception("Method not good")


def Newton_prog(xs, fxs):
	n = len(xs)
	h = (xs[-1] - xs[0])/(len(xs)-1)

	# Schema Aitken
	ait = [[] for _ in range(n)]
	ait[0] = [fxs[i+1]-fxs[i] for i in range(n-1)]
	for i in range(1, n-1):
		ait[i] = [ait[i-1][j+1]-ait[i-1][j] for j in range(n-1-i)]
	# for i in range(n-1):
	# 	print(ait[i])

	def Ln(x):
		t = (x - xs[0])/h
		ss = [t]
		for k in range(2, n):
			ss.append(ss[-1]*(t-k+1)/k)

		# print(prev, fxs[0])
		l = fxs[0]
		l += sum([ait[k][0]*ss[k] for k in range(n-1)])
		return l
	return Ln

def cmm_patrate(xs, fxs):
	global M
	# b = np.stack([xs]+[np.ones(len(xs)) for _ in range(3)]).T
	ff = fxs
	# return np.linalg.lstsq(b, ff)[0]
	return np.polyfit(xs, fxs, deg=M)


def Horner(polinom, x):
	d = 0
	for i in range(len(polinom)):
		d = polinom[i] + d * x
	return d

def plott(xs, actual, **kwargs):
	for idx, k in enumerate(kwargs):
		plt.plot(xs, [kwargs[k](x) for x in xs], linewidth=5*(len(kwargs)-idx), label=k)
	plt.plot(xs, [actual_f(x) for x in xs], linewidth=1, color='blue', label='actual')
	# plt.plot(xs, fxs, linewidth=1, color='blue', label='actual')

	plt.legend()
	plt.show()

def decl_f(x):
	return x ** 4 - 12 * (x ** 3) + 30 * (x ** 2) + 12

M = 4
xs, fxs, actual_f = read_data()
print("xs, fxs:", xs, fxs)
x_ = xs[0] + (xs[1] - xs[0])/2
print("==== NEWTON")
ll_newt = Newton_prog(xs, fxs)
# print(ll(x_), actual_f(x_))
print(f"{ll_newt(x_)=}\n{norma(ll_newt(x_) - actual_f(x_))=}")

print("==== CMM PATRATE")
poly = cmm_patrate(xs, fxs)
print(poly)
print(f"{Horner(poly, x_)}\n{norma(Horner(poly, x_) - actual_f(x_))=}")
print("Sum norme:", sum([norma(Horner(poly, x) - actual_f(x)) for x in xs]))
ll_cmm = lambda x: Horner(poly, x)

plott(xs, actual=actual_f, newton=ll_newt, patrate=ll_cmm, )