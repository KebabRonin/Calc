
def Horner(polinom, x):
	d = 0
	for i in range(len(polinom)):
		d = polinom[i] + d * x
	return d