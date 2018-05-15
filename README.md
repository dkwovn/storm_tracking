# storm_tracking
Extratropical storms are related to extreme hazadous weather and long-range air pollutants transport. They usually manifest themseves as "blobs" of high vorticity. Everyday, there are hundreds of such "blobs" floating around the globe but only a small fraction of them are legitimate storms. A computer vision algorithm is applied to find the storms by requiring smooth trajectories of those candidate storms. Such problem is also known as correspondence problem to track multiple moving objects from continuious image series.

~~~~~~~~~~
The computer vision algorithm was proposed by Salari and Sethi (1990) and introduced to the field of atmospheric sciences by Hodges (1995)
~~~~~~~~~~

# Step 1
The first step is to find those objects. Depth-First Search is used to find those candidate storms and then their centroids are calculated as "feature points".

# Step 2
After the feature points are found, hypothetical trajectories are formed from these points. A cost function quantifying the smoothness of the trajectory is associated with each trajectory. Then a greedy exchange algorithm is used to swap feature points from hypothetical trajectories until the sum of all the cost reaches a local minimum. 

