#!/opt/local/bin/python

content = [line.rstrip() for line in open("POSCAR", 'r').readlines()]

scale=float(content[1])

a=[float(val) for val in content[2].split()]
b=[float(val) for val in content[3].split()]
c=[float(val) for val in content[4].split()]

vol=scale**3*(a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]))

if vol<0: print(-vol)
else: print(vol)
