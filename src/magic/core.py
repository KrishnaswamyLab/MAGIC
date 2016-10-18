import numpy as np
import time
import sys
import warnings
from numpy import matlib

from sklearn.neighbors import NearestNeighbors
from scipy import sparse, stats
from numpy import linalg
from scipy.sparse.linalg import norm
import networkx as nx

def wishbone(data, s, k=15, l=15, num_graphs=1, num_waypoints=250, 
	verbose=True, metric='euclidean', voting_scheme='exponential', 
	branch=True, flock_waypoints=2, band_sample=False, partial_order=[],
	search_connected_components=True):

	if verbose:
		print('Building lNN graph...')
	

	# Construct nearest neighbors graph
	start = time.process_time()
	nbrs = NearestNeighbors(n_neighbors=l+1, metric=metric).fit(data)  
	lnn = nbrs.kneighbors_graph(data, mode='distance' ) 
	lnn = np.transpose(lnn)
	print('lNN computed in : %.2f seconds' % (time.process_time()-start))

	#set up return structure
	trajectory = []
	waypoints = []
	branches = []
	bas = []

	# generate klNN graphs and iteratively refine a trajectory in each
	for graph_iter in range(num_graphs):
		if k!=l:
			klnn = _spdists_klnn(lnn, k, verbose)
		else:
			klnn = lnn

		# Make the graph undirected
		klnn = _spdists_undirected(klnn)
		klnn.setdiag(0)
		klnn.eliminate_zeros()

		#run traj. landmarks
		traj, dist, iter_l, paths_l2l = _trajectory_landmarks( klnn, data, [s], num_waypoints, partial_order, verbose, metric, flock_waypoints, band_sample, branch)
		if branch:
			if verbose:
				print ('Determining branch point and branch associations...')
			RNK, bp, diffdists, Y = _splittobranches(traj, traj[0], data, iter_l, dist, paths_l2l)


		# calculate weighed trajectory
		W_full = _weighting_scheme(voting_scheme, dist)

		if branch:
			W = _muteCrossBranchVoting(W_full, RNK, RNK[s], iter_l, Y)
		else:
			W = W_full

		
		# save initial solution - start point's shortest path distances
		t = traj[0, :]
		t = [t, np.sum(np.multiply(traj, W), axis=0)]

		# iteratively realign trajectory (because landmarks moved)
		converged, user_break, realign_iter = False, False, 1
		if verbose:
			print('Running iterations...')

		while converged == False and user_break == False:
			realign_iter = realign_iter + 1
			print('Iteration: %d' % realign_iter)

			np.copyto(traj, dist)
			traj = _realign_trajectory(t, dist, iter_l, traj, 0, len(dist), realign_iter)

			if branch:
				RNK, bp, diffdists, Y = _splittobranches(traj, traj[0],data, iter_l, dist,paths_l2l)
				W = _muteCrossBranchVoting(W_full, RNK, RNK[s], iter_l,Y)
			# calculate weighed trajectory
			t.append(np.sum(np.multiply(traj, W), axis=0))
			
			#check for convergence
			fpoint_corr = stats.pearsonr(np.transpose(t[realign_iter]), np.transpose(t[realign_iter - 1]))[0]
			if verbose:
				print('Correlation with previous iteration:  %.4f' % fpoint_corr)
			converged = fpoint_corr > 0.9999
		
			if (realign_iter % 16) == 0:
				# break after too many alignments - something is wrong
				user_break = True
				print('\nWarning: Force exit after ' + str(realign_iter) + ' iterations')

		print(str(realign_iter-1) + ' realignment iterations')

		# save final trajectory for this graph			
		iter_traj = t[realign_iter][:]
		# Normalize the iter_trajectory
		iter_traj = (iter_traj - iter_traj.min()) / (iter_traj.max() - iter_traj.min())
		trajectory.append(iter_traj)
		waypoints.append(iter_l)
		
		if branch:
			# Recalculate branches post reassignments
			RNK, bp, diffdists, Y = _splittobranches(traj, traj[0], data, iter_l, dist,paths_l2l)
			branches.append(RNK)
			bas.append(Y)

			# Set return values
			branches = branches[0]
			bas = bas[0]
		else:
			branches = None
			bas = None

		trajectory = trajectory[0]
		waypoints = waypoints[0]
		return dict(zip(['Trajectory', 'Waypoints', 'Branches', 'BAS'],
		  [trajectory, waypoints, branches, bas]))









# All the utility functions
#randomly removing l-k edges for each iteration of graph_num
def _spdists_klnn(spdists,k,verbose):
    
    #TODO - avoid converting back to non-sparse
    if sparse.issparse(spdists):
        spdists = spdists.toarray()

    for idx in range(len(spdists)):
        neighs = np.where(spdists[:,idx] != 0)[0]
        l = len(neighs)
        remove_indices = neighs[np.random.choice(l,l-k,replace=False)]
        spdists[remove_indices,idx] = 0
    
    return sparse.csr_matrix(spdists)

#making graph undirected
def _spdists_undirected(spdists):

	with warnings.catch_warnings():
		warnings.simplefilter('ignore')  # Catch the sparse efficiency warning
		[i,j] = np.nonzero(spdists)
		spdists[(j,i)] = spdists[np.nonzero(spdists)]
	return spdists

#calculated weighting matrix
def _weighting_scheme(voting_scheme, dist):
    if voting_scheme == 'uniform':
        W_full = np.ones(dist.shape())
    elif voting_scheme == 'exponential':
        dist_transposed = dist.transpose()
        std = np.empty(len(dist_transposed))
        for i in range(len(dist_transposed)):
            std[i] = np.std(dist_transposed[i])
        sdv = np.mean (std)*3
        W_full = np.exp( -.5 * np.power((dist / sdv), 2))
    elif voting_scheme == 'linear':
        W_full = matlib.repmat(dist.max(), len(dist), 1) - dist

    # The weighing matrix must be a column stochastic operator
    W_full_transposed = W_full.transpose()
    W_full_sum = np.empty(len(W_full_transposed))
    for i in range(len(W_full_transposed)):
        W_full_sum[i] = W_full_transposed[i].sum()

    W_full = np.divide(W_full, matlib.repmat( W_full_sum, len( W_full), 1 ))
    return W_full


def _realign_trajectory(t, dist, l, traj, start_range, end_range, realign_iter):
    for idx in range(start_range, end_range):
        #find position of landmark in previous iteration
        idx_val = t[realign_iter - 1][l[idx]]
        #convert all cells before starting point to the negative
        before_indices = np.where(t[realign_iter - 1] < idx_val)[0]
        for before_idx in before_indices:
            traj[idx][before_idx] = -dist[idx][before_idx]
        #set zero to position of starting point
        traj[idx] = traj[idx] + idx_val

    traj = np.subtract(traj, np.min(traj))
    return traj


#determining initial trajectory
def _trajectory_landmarks(spdists, data, s, waypoints, partial_order, 
    verbose, metric, flock_waypoints, band_sample, branch):
    
    #if given a list of possible starting points, choose one
	if verbose:
		print('Determining waypoints if not specified...')
	start = time.process_time()

	if len(s) > 1:
		s = np.random.choice(s,1,replace=False)

	#if not given landmarks list, decide on random landmarks
	dijkstra = sparse.csgraph.dijkstra(spdists, directed=False, indices=s)
	if isinstance(waypoints, int):
		n_opts = np.arange(0, len(data))
		if band_sample:
			n_opts = []
			window_size = 0.1
			max_dist = dijkstra.max()
			prc = 0.998
			while prc > 0.08:
				band = [i for i in range(len(dijkstra)) if dijkstra[i]>= ((prc-window_size)*max_dist) and dijkstra[i]<=prc*max_dist]
				n_opts = np.append(n_opts, np.random.choice( band, min(len(band), waypoints[0] - 1 - len(partial_order)), replace=False ))
				prc = prc - window_size

		waypoints = np.random.choice(n_opts, waypoints-1-len(partial_order), replace=False)
		waypoints = [int(waypoints[i]) for i in range(len(waypoints))]

		if branch:
			tailk=30
			tailband = np.where(dijkstra>=np.percentile(dijkstra, 85))[0]
			tailk = min(len(tailband), tailk, np.floor(len(waypoints)/2))
			to_replace = np.random.randint(len(waypoints)-1, size=tailk)
			tailband_sample = np.random.choice( tailband, size=tailk, replace=False)
			for i in range(len(to_replace)):
				waypoints[to_replace[i]] = tailband_sample[i]

	# Flock wayoints
	if flock_waypoints:
		nbrs = NearestNeighbors(n_neighbors=20, metric=metric).fit(data) 
	for f in range(flock_waypoints):
		IDX = nbrs.kneighbors([data[i] for i in waypoints], return_distance=False)
		for i in range(len(waypoints)): 
			med_data = np.median(data[IDX[i,:],:],axis=0)
			waypoints[i] = nbrs.kneighbors(med_data.reshape(1, -1), n_neighbors=1, return_distance=False)[0][0]


	if s not in partial_order:
		partial_order = np.append(s,partial_order) #partial_order includes start point

	l = np.append(partial_order,waypoints) #add extra landmarks if user specified
	l = l.astype(int)

	# # calculate all shortest paths
	print('Determining shortest path distances and perspectives....')
	graph = nx.Graph(spdists)

	for i in range(len(l)):
		temp = nx.single_source_dijkstra(graph, l[i])
		paths = temp[1]
		# dist[i, :] = temp[0].values()

		if i == 0:
			dist = [list(temp[0].values())]
			paths_l2l = [[paths[li] for li in l]]
		else:
			dist.append(list(temp[0].values()))
			paths_l2l.append([paths[li] for li in l])

		unreachable = np.where(dist[i]==np.inf)[0]
		if(len(unreachable) > 0):
			dist[i][unreachable] = max(max(dist))
		if(verbose):
			sys.stdout.write('.')
	print("")

	#adjust paths according to partial order by redirecting
	dist = np.array(dist)
	nPartialOrder = len(partial_order)   
	for radius in range(1,nPartialOrder+1): 
		for landmark_row in range(1,nPartialOrder+1):
			if landmark_row + radius <= nPartialOrder:
				a = landmark_row
				b = landmark_row + (radius-1)
				c = landmark_row + radius
				dist[a-1][partial_order[c-1]] = dist[a-1][partial_order[b-1]] + dist[b-1][partial_order[c-1]]
			if landmark_row - radius >= 1:
				a = landmark_row
				b = landmark_row - (radius-1)
				c = landmark_row - radius
				dist[a-1][partial_order[c-1]] = dist[a-1][partial_order[b-1]] + dist[b-1][partial_order[c-1]]

	#align to dist_1
	traj = np.empty_like (dist)
	np.copyto(traj,dist)

	for idx in range(1,len(partial_order)):
		closest_landmark_row = argmin(dist,axis=0) #closest landmark will determine directionality
		traj[idx, closest_landmark_row < idx] = -dist[idx][closest_landmark_row < idx]
		traj[ idx, : ] = traj[ idx, : ] + dist[0][l[idx]]


	if len(l) > len(partial_order):
		traj = _realign_trajectory(dist, dist, l, traj, len(partial_order), len(l), 1)

	print('Time for determining distances and perspectives: %.2f seconds' % (time.process_time()-start))

	return traj, dist, l, paths_l2l



def _splittobranches(trajs, t, data, landmarks, dist, paths_l2l):

	proposed = matlib.repmat(t[landmarks], len(trajs), 1)
	reported = np.array([trajs[i][landmarks] for i in range(len(landmarks))])

	#square matrix of the difference of perspectives landmark to landmark
	diffdists = np.absolute(reported - proposed)
	diffdists = (diffdists.T + diffdists)/2

	# get second eigen vector of diffdists
	EigenVals, EigenVecs = linalg.eig(diffdists)
	sorted_idxs = np.argsort(np.absolute(EigenVals))[::-1]
	
	evec2 = np.multiply(EigenVecs[:, sorted_idxs[1]], -1)
	idx = np.argsort(evec2)
	evec2[np.where(evec2 == 0)[0]] = 0

	#assign last positive 5 and last negative 5
	t_l = t[landmarks]
	b1_suspect = np.where(evec2 < 0)[0] #suspects for branch 1
	b2_suspect = np.where(evec2 > 0)[0] #suspects for branch 2

	b1_sorted_inds = np.argsort(t_l[b1_suspect])[::-1]
	b2_sorted_inds = np.argsort(t_l[b2_suspect])[::-1]

	c = np.ones(len(landmarks))
	c[b1_suspect[b1_sorted_inds[:min(5, len(b1_sorted_inds))]]] = 2
	c[b2_suspect[b2_sorted_inds[:min(5, len(b2_sorted_inds))]]] = 3

	trunk = c[0]

	c_branch = np.setdiff1d(np.unique(c).T, trunk)

	#return with current branches if branch is not found
	if(len(c_branch) == 1):
		print('Branch not found\n')
		I = np.absolute(dist[:len(landmarks)][:]).min()
		RNK = c[I]
		Y = np.zeros((len(RNK)))
		pb = t[landmarks[(np.where(c == c_branch[0])[0][0])]]

		return RNK, pb, diffdists, Y

	brancha = np.where(c == c_branch[0])[0]
	branchb = np.where(c == c_branch[1])[0]
	paths_branch_a = [paths_l2l[i] for i in brancha]
	paths_branch_b = [paths_l2l[i] for i in branchb]
	
	fork_p = []
	for i in range(len(paths_branch_a)):
		paths_branch_a_to_b = [paths_branch_a[i][k] for k in branchb]
		for j in range(len(paths_branch_a_to_b)):
			if paths_branch_a_to_b[j]:
				fork_p.append(t[paths_branch_a_to_b[j]].min())
			else:
				print('no path from l:' + str(brancha[i]) + ' to l:' + str(branchb[j]))

	for i in range(len(paths_branch_b)):
		paths_branch_b_to_a = [paths_branch_b[i][k] for k in brancha]
		for j in range(len(paths_branch_b_to_a)):
			if paths_branch_b_to_a[j]:
				fork_p.append(t[paths_branch_b_to_a[j]].min())
			else:
				print('no path from l:' + str(branchb[i]) + ' to l:' + str(brancha[j]))

	#reassign to clusters based on branch point
	pb = np.percentile(fork_p, 10)

	# #further adjust the branch point until a minimum delta is seen

	# #landmarks before branch point
	# lm_bp = np.where(t[landmarks] < pb)[0]
	
	# #sort according to trajectory
	# lm_order = np.argsort(t[landmarks[lm_bp]])[::-1]

	# #track towards beginning until the difference to the nearest landmark
	# #of the other branch is less than delta
	# for i in range(len(lm_order)):
	# 	iter_lm = lm_bp[lm_order[i]]
	# 	next_lm = np.intersect1d(np.where(np.sign(evec2[lm_bp]) != np.sign(evec2[iter_lm]))[0], np.where(t[landmarks[lm_bp]].T < t[landmarks[iter_lm]])[0])[0]
	# 	if np.absolute(evec2[iter_lm] - evec2[lm_bp[next_lm]]) < 0.05:
	# 		if i > 0:
	# 			pb = t[landmarks[lm_bp[lm_order[i]]]]
	# 		break

	c_new = c
	I = np.argmin(np.absolute(dist[0:landmarks.size, :]), axis=0)
	RNK = c_new[I]
	c_new[np.where(t[landmarks].T <= pb)[0]] = c[0]
	c_new[np.intersect1d(np.where(evec2 < 0)[0], np.where(t[landmarks].T >= pb)[0])] = c_branch[0]
	c_new[np.intersect1d(np.where(evec2 > 0)[0], np.where(t[landmarks].T >= pb)[0])] = c_branch[1]


	#compute affinity matrix over landmark distances
	n = len(dist[0]) #num points
	sigma = .1*np.std(dist, ddof=1)
	Aff = np.exp(np.multiply(-0.5*(1/np.power(sigma, 2)), np.power(dist, 2)))

	#make aff matrix a stochastic operator
	Stoch = np.divide(Aff, np.sum(Aff, axis=0))
	Stoch[np.where(np.isnan(Stoch))[0]] = np.nanmin(Stoch)
	Y = np.multiply(np.dot(Stoch.T, evec2), np.power(t.T, 0.7))

	#for each datapoint, find closest landmark
	I = np.argmin(np.absolute(dist[0:landmarks.size, :]), axis=0)
	RNK = c_new[I]
	
	return RNK, pb, diffdists, Y 


def _muteCrossBranchVoting(W, RNK, trunk_id, landmarks, Y):
	#range between -1 and 1
	Y_scale = np.subtract(Y, np.median(Y[landmarks]))
	indices = np.where(Y_scale < 0)[0]
	if(len(indices) > 0):
		Y_scale[indices] = np.divide(Y_scale[indices], np.absolute(Y_scale[indices]).max())
	indices = np.where(Y_scale > 0)[0]
	if(len(indices) > 0):
		Y_scale[indices] = np.divide(Y_scale[indices], Y_scale[indices].max())
	Y_pos = np.absolute(Y_scale).T

	W_test = W
	b = np.std(Y_scale, ddof=1)

	for i in range(Y.size):
		crossb = np.where((np.sign(Y_scale[i]) != np.sign(Y_scale[landmarks])) == True)[0]
		temp = np.exp(np.divide(-0.5*np.power(Y_pos[landmarks[crossb]], 2), b))
		temp = [max(np.exp(np.divide(-0.5*np.power(Y_pos[int(i)], 2), b)), i) for i in temp]
		
		for j in range(len(crossb)):
		 	W_test[crossb[j]][i] = np.multiply(W[crossb[j]][i], temp[j])
	
	W = np.divide(W_test, np.sum(W_test, axis=0))
	return W
