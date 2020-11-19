#!/opt/local/bin/python

content = [line.rstrip() for line in open("POSCAR", 'r').readlines()]

scale=float(content[1])

for i in range(2,5):
	vec=[float(val) for val in content[i].split()]
	res=scale*(vec[0]**2+vec[1]**2+vec[2]**2)**(1/2)

	if res<0: print(-res)
	else: print(res)
