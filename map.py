import matplotlib.pyplot as plt

with open('data/map_data.txt') as f:
  lines = f.readlines()

raw = [l.rstrip('\n').split('\t') for l in lines]
x = [float(r[0]) for r in raw]
y = [float(r[1]) for r in raw]
labels = [r[2] for r in raw]

plt.scatter(x,y)
for i, l in enumerate(labels):
  plt.annotate(l, xy=(x[i],y[i]), size=15)

plt.axis('equal')
plt.show()
