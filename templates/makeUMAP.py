import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import umap
import random

random.seed(111)

##unsupervised full
csv = np.genfromtxt("dims.txt", delimiter=",") #read in SVD left singular vectors
fit = umap.UMAP(n_neighbors=50, min_dist = 0.6, random_state=111) #general whole integrated map umaps //  number of neighbors: 150 for CD10+/- data // 50 for other data
#fit = umap.UMAP(n_neighbors = 80, min_dist = 1, random_state=111) #for "lineage" UMAP - num. neighbors: 46 for CD10- mesenchymal data, 89 for PDGFRb mesenchymal data
u = fit.fit_transform(csv)
np.savetxt("xx_umapCoords.csv", u, delimiter = ",") #save the UMAP coords files


