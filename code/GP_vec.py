## inf_utility.py ##
'''
Code to go from training data to utility function using two independent GPs with Squared exponential kernel

'''
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import GP_deriv
import gpflow
import time
from utils import KL_div

def SE_covariance_mult_dim(x1,x2,hypInput,hypOutput):
	"""
	Calculates the covariance matrix between two vectors x1 and x2 according to the sqared exponential for multiple dimensions
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,d) and (M,d)
		hypInput: Input scale length (length scale) (d,)
		hypOutput: Output scale length
	---
	Output: 
		K = covariance shape (N,M)
	"""
	dim = np.shape(x1)[1]
	N = np.shape(x1)[0] 
	M = np.shape(x2)[0]
	X = np.zeros((N,M))
	
	for d in range(dim):
		X += (hypInput[d] ** -2) * np.power(np.subtract.outer(x1[:,d],x2[:,d]),2).reshape((N,M))

	K = (hypOutput**2)*np.exp(-0.5*X)
	return K

def SE_covariance_mult_dim_derivative_x1(x1,x2,hypInput,hypOutput):
	"""
	Calculates the first derivative the sqared exponential kernel in multiple dimensions w.r.t. x1
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,d) and (M,d)
		hypInput: Input scale length (length scale) (d,)
		hypOutput: Output scale length
	---
	Output: 
		K_d: covariance shape (d,N,M), where d is the derivative in that direction
	"""
	dim = np.shape(x1)[1]
	N = np.shape(x1)[0] 
	M = np.shape(x2)[0]
	X = np.zeros((N,M))
	K_d = np.zeros((dim,N,M))

	for d in range(dim):
		X += (hypInput[d] ** -2) * np.power(np.subtract.outer(x1[:,d],x2[:,d]),2).reshape((N,M))
	
	K = (hypOutput**2)*np.exp(-0.5*X)

	for d in range(dim):
		K_d[d,:,:] = - (hypInput[d] ** -2) * np.subtract.outer(x1[:,d],x2[:,d]).reshape((N,M)) * K

	return K_d

class Pred_vector_field:
	""" Go from training data to scalar field """
	
	def __init__(self, x_train, y_train, t_train, x_test, y_test, t_test, time_hyp=20,xy_hyp=1,var_hyp=1,noise_hyp=0.01):
		"""
		Need to define the constructor
		---
		Inputs:
			t_train: np.array of shape (N,1)
			t_test: np.array of shape (M,1)
			y_train: np.array of shape (N,1)
			x_train: np.array of shape (N,1)
			x_test: np.array of shape (Ntx,1)
			y_test: np.array of shape (Nty,1)
			time_hyp: Set the time hyperparameter 
		""" 
		t_tot = time.time()

		######################
		self.N_train = np.shape(t_train)[0]
		self.N_test = np.shape(t_test)[0]
		self.N_test_x = np.shape(x_test)[0]
		self.N_test_y = np.shape(y_test)[0]
		self.N_test_t = np.shape(t_test)[0]
		self.t_train = t_train
		self.t_test = t_test 
		self.x_train = x_train
		self.y_train = y_train
		self.x_test = x_test
		self.y_test = y_test
		self.time_hyp = time_hyp

		## Optimise the GPs
		t0 = time.time()
		######################
		K_y = GPy.kern.RBF(1)
		self.m_y = GPy.models.GPRegression(t_train,y_train,K_y)
		self.m_y.optimize()
		K_x = GPy.kern.RBF(1)
		self.m_x = GPy.models.GPRegression(t_train,x_train,K_x)
		self.m_x.optimize()
		######################
		print('Wall time for first layer GPs: %0.2fs' % (time.time() - t0))

		## Derivative GPs

		# Using hyperparameters from GPy independent
		hyp_x = [self.m_x.rbf.lengthscale,np.sqrt(self.m_x.rbf.variance),np.sqrt(self.m_x.Gaussian_noise.variance)]
		hyp_y = [self.m_y.rbf.lengthscale,np.sqrt(self.m_y.rbf.variance),np.sqrt(self.m_y.Gaussian_noise.variance)]
		
		## Initialise GPs
		self.GP_x = GP_deriv.DerivGP_2nd(t_train=self.t_train.reshape((self.N_train,)),t_test=self.t_train.reshape((self.N_train,)),
			y_train=x_train.reshape((self.N_train,)),hyp=hyp_x)
		self.GP_y = GP_deriv.DerivGP_2nd(t_train=self.t_train.reshape((self.N_train,)),t_test=self.t_train.reshape((self.N_train,)),
			y_train=y_train.reshape((self.N_train,)),hyp=hyp_y)

		### Second stage of predicting the vector field
		self.X,self.Y,self.T = np.meshgrid(self.x_test,self.y_test,self.t_test)
		y_t = self.Y.reshape((self.N_test_x*self.N_test_y*self.N_test_t))
		x_t = self.X.reshape((self.N_test_x*self.N_test_y*self.N_test_t))
		t_t = self.T.reshape((self.N_test_x*self.N_test_y*self.N_test_t))
		# Test points stacked 
		self.x_3d_test = np.vstack([x_t,y_t,t_t]).T

		#### WEIGHTS AND REMOVAL OF POINTS ####
		self.pred_dif_x = self.GP_x.pred_deriv_mean_xx.reshape((self.N_train,1))    # Second derivatives
		self.pred_dif_y = self.GP_y.pred_deriv_mean_xx.reshape((self.N_train,1))    # Second derivatives
		if norm:
			self.pred_dif_x = self.pred_dif_x * x_std
			self.pred_dif_y = self.pred_dif_y * y_std           

		# Training points
		self.N_train    = x_der.shape[0]
		self.pred_dif_x = x_der.reshape((self.N_train,1))
		self.pred_dif_x = np.array(self.pred_dif_x,dtype=np.float64)
		self.pred_dif_y = y_der.reshape((self.N_train,1))
		self.pred_dif_y = np.array(self.pred_dif_y,dtype=np.float64)
		self.x_3d_train = np.vstack([x,y,t]).T
		self.x_3d_train = np.array(self.x_3d_train,dtype=np.float64)

		# 3D GPs
		t0 = time.time()
		####################################
		K_3d_x = gpflow.kernels.RBF(input_dim=3, ARD=True)
		self.m_3d_x = gpflow.models.GPR(self.x_3d_train,self.pred_dif_x,K_3d_x)
		self.m_3d_x.compile()
		opt_3dx = gpflow.train.ScipyOptimizer()
		opt_3dx.minimize(self.m_3d_x)           
		
		K_3d_y = gpflow.kernels.RBF(input_dim=3, ARD=True)
		self.m_3d_y = gpflow.models.GPR(self.x_3d_train,self.pred_dif_y,K_3d_y)
		self.m_3d_y.compile()
		opt_3dy = gpflow.train.ScipyOptimizer()
		opt_3dy.minimize(self.m_3d_y)           
		###################################
		print('Wall time for second layer GPs: %0.2fs' % (time.time() - t0))

		# Test derivative field:
		t0 = time.time()
		####################################
		self.V,self.S = self.vector_field()
		####################################
		print('Wall time for predictive vector field: %0.2fs' % (time.time() - t0))

		## Probabilistic derivatives of vector field
		t0 = time.time()
		####################################
		(self.V_xx,self.V_xy,self.V_yx,self.V_yy), (self.S_xx,self.S_xy,self.S_yx,self.S_yy) = self.first_derivative()
		####################################
		print('Wall time for second order derivatives: %0.2fs' % (time.time() - t0))

		# Time order test points
		t0 = time.time()
		####################################
		self.loc_time_indexed, self.V_time_indexed, self.V_2d, self.S_time_indexed, self.S_2d = self.time_order()
		####################################
		print('Wall time for rearranging the arrays to be time ordered: %0.2fs' % (time.time() - t0))

		self.V_xx_time_indexed = self.V_2d[:,:,0]
		self.V_xy_time_indexed = self.V_2d[:,:,1]
		self.V_yx_time_indexed = self.V_2d[:,:,2]
		self.V_yy_time_indexed = self.V_2d[:,:,3]

		self.S_xx_time_indexed = self.S_2d[:,:,0]
		self.S_xy_time_indexed = self.S_2d[:,:,1]
		self.S_yx_time_indexed = self.S_2d[:,:,2]
		self.S_yy_time_indexed = self.S_2d[:,:,3]

		# Calculate div and curl
		####################################
		self.div = self.V_xx_time_indexed + self.V_yy_time_indexed
		self.div_var = self.S_xx_time_indexed + self.S_yy_time_indexed
		self.curl = self.V_yx_time_indexed - self.V_xy_time_indexed
		self.curl_var = self.S_yx_time_indexed + self.S_xy_time_indexed
		# Reorder for plotting
		self.loc_time_ind_mesh_x = []   
		self.loc_time_ind_mesh_y = []   
		self.div_mesh = []          
		self.div_var_mesh = []          
		self.curl_mesh = []         
		self.curl_var_mesh = []         
		self.S_mesh_x = []              
		self.S_mesh_y = []              
		
		t0 = time.time()
		####################################
		self.calc_div_curl()
		####################################
		print('Wall time for rearranging the arrays to be a mesh: %0.2fs' % (time.time() - t0))
		
		# KL-divergence at each point in the space:
		t0 = time.time()
		####################################

		# Prior: GP(0,K_te_te_der)
		hyp_x = [self.m_3d_x.kern.lengthscales.value,np.sqrt(self.m_3d_x.kern.variance.value),np.sqrt(self.m_3d_x.likelihood.variance.value)]
		hyp_y = [self.m_3d_y.kern.lengthscales.value,np.sqrt(self.m_3d_y.kern.variance.value),np.sqrt(self.m_3d_y.likelihood.variance.value)]
		prior = (hyp_x[1]/hyp_x[0][0])**2 + (hyp_y[1]/hyp_y[0][1])**2
		prior_var = prior * np.ones_like(self.div_mesh[0])

		N_0 = [np.zeros_like(self.div_mesh[0]),prior_var]

		self.KL_div_mesh = []
		self.KL_div_mesh_signed = []
		for t in range(self.N_test):
			# Posterior:
			N_1 = [self.div_mesh[t],self.div_var_mesh[t]]
			self.KL_div_mesh.append(KL_div(N_0,N_1))
			self.KL_div_mesh_signed.append(KL_div(N_0,N_1)*np.sign(N_1[0]))

		###################################
		print('Wall time for KL-divergence: %0.2fs' % (time.time() - t0))
		
		###################################
		print('Total wall time: %0.2fs' % (time.time() - t_tot))


	def vector_field(self):
		"""
		Calculates the predictive vector field from tthe optimised hyperparameters
		---
		Output: 
			V: vector field (2,N_x*N_y*N_t)
			S: standard deviation (2,N_x*N_y*N_t)
		"""

		prior_mu = 0
		data_mu = 0
		V = np.zeros((2,self.N_test_x*self.N_test_y*self.N_test_t))
		S = np.zeros((2,self.N_test_x*self.N_test_y*self.N_test_t))
		L = self.N_test_x*self.N_test_y*self.N_test_t

		# X direction:
		hyp_x = [self.m_3d_x.kern.lengthscales.value,np.sqrt(self.m_3d_x.kern.variance.value),np.sqrt(self.m_3d_x.likelihood.variance.value)]
		K_tr_tr_x = SE_covariance_mult_dim(self.x_3d_train,self.x_3d_train,hyp_x[0],hyp_x[1]) + np.identity(self.N_train)*hyp_x[2]**2
		K_te_tr_x = SE_covariance_mult_dim(self.x_3d_test,self.x_3d_train,hyp_x[0],hyp_x[1])
		K_te_te_x = SE_covariance_mult_dim(self.x_3d_test,self.x_3d_test,hyp_x[0],hyp_x[1])

		V[0,:] = GP_deriv.Pred_mean(prior_mu,K_tr_tr_x,K_te_tr_x,self.pred_dif_x,data_mu).reshape((L,))
		S[0,:] = np.diagonal(GP_deriv.Pred_var(K_te_te_x, K_tr_tr_x, K_te_tr_x))

		# Y direction:
		hyp_y = [self.m_3d_y.kern.lengthscales.value,np.sqrt(self.m_3d_y.kern.variance.value),np.sqrt(self.m_3d_y.likelihood.variance.value)]
		K_tr_tr_y = SE_covariance_mult_dim(self.x_3d_train,self.x_3d_train,hyp_y[0],hyp_y[1]) + np.identity(self.N_train)*hyp_y[2]**2
		K_te_tr_y = SE_covariance_mult_dim(self.x_3d_test,self.x_3d_train,hyp_y[0],hyp_y[1])
		K_te_te_y = SE_covariance_mult_dim(self.x_3d_test,self.x_3d_test,hyp_y[0],hyp_y[1])

		V[1,:] = GP_deriv.Pred_mean(prior_mu,K_tr_tr_x,K_te_tr_x,self.pred_dif_y,data_mu).reshape((L,))
		S[1,:] = np.diagonal(GP_deriv.Pred_var(K_te_te_y, K_tr_tr_y, K_te_tr_y))

		return V,S
	
	def first_derivative(self):
		"""
			Calculates the predictive mean derivative in a given direction
		---
		Output: 
			V_xx:  predictive mean in x direction derivative wrt x
			V_xy:  predictive mean in x direction derivative wrt y
			V_yx:  predictive mean in y direction derivative wrt x
			V_yy:  predictive mean in y direction derivative wrt y
			S_xx:  predictive variance in x direction derivative wrt x
			S_xy:  predictive variance in x direction derivative wrt y
			S_yx:  predictive variance in y direction derivative wrt x
			S_yy:  predictive variance in y direction derivative wrt y
		"""

		prior_mu = 0
		data_mu = 0

		# X direction:
		hyp_x = [self.m_3d_x.kern.lengthscales.value,np.sqrt(self.m_3d_x.kern.variance.value),np.sqrt(self.m_3d_x.likelihood.variance.value)]
		K_tr_tr_x = SE_covariance_mult_dim(self.x_3d_train,self.x_3d_train,hyp_x[0],hyp_x[1]) + np.identity(self.N_train)*hyp_x[2]**2
		K_te_tr_der_x = SE_covariance_mult_dim_derivative_x1(self.x_3d_test,self.x_3d_train,hyp_x[0],hyp_x[1])
		
		# Select appropriate lengthscales for the variance
		h_xx = [self.m_3d_x.kern.lengthscales.value[0],np.sqrt(self.m_3d_x.kern.variance.value),np.sqrt(self.m_3d_x.likelihood.variance.value)]
		h_xy = [self.m_3d_x.kern.lengthscales.value[1],np.sqrt(self.m_3d_x.kern.variance.value),np.sqrt(self.m_3d_x.likelihood.variance.value)]
		
		V_xx = GP_deriv.Pred_mean(prior_mu,K_tr_tr_x,K_te_tr_der_x[0],self.pred_dif_x,data_mu)
		S_xx = np.diag(GP_deriv.Pred_var_deriv(K_tr_tr_x,K_te_tr_der_x[0],K_te_tr_der_x[0].T,h_xx))
		V_xy = GP_deriv.Pred_mean(prior_mu,K_tr_tr_x,K_te_tr_der_x[1],self.pred_dif_x,data_mu)
		S_xy = np.diag(GP_deriv.Pred_var_deriv(K_tr_tr_x,K_te_tr_der_x[1],K_te_tr_der_x[1].T,h_xx))

		# Y direction:
		hyp_y = [self.m_3d_y.kern.lengthscales.value,np.sqrt(self.m_3d_y.kern.variance.value),np.sqrt(self.m_3d_y.likelihood.variance.value)]
		K_tr_tr_y = SE_covariance_mult_dim(self.x_3d_train,self.x_3d_train,hyp_y[0],hyp_y[1]) + np.identity(self.N_train)*hyp_y[2]**2
		K_te_tr_der_y = SE_covariance_mult_dim_derivative_x1(self.x_3d_test,self.x_3d_train,hyp_y[0],hyp_y[1])

		# Select appropriate lengthscales for the variance
		h_yx = [self.m_3d_y.kern.lengthscales.value[0],np.sqrt(self.m_3d_y.kern.variance.value),np.sqrt(self.m_3d_y.likelihood.variance.value)]
		h_yy = [self.m_3d_y.kern.lengthscales.value[1],np.sqrt(self.m_3d_y.kern.variance.value),np.sqrt(self.m_3d_y.likelihood.variance.value)]

		V_yx = GP_deriv.Pred_mean(prior_mu,K_tr_tr_y,K_te_tr_der_y[0],self.pred_dif_y,data_mu)
		S_yx = np.diag(GP_deriv.Pred_var_deriv(K_tr_tr_y,K_te_tr_der_y[0],K_te_tr_der_y[0].T,h_yx))
		V_yy = GP_deriv.Pred_mean(prior_mu,K_tr_tr_y,K_te_tr_der_y[1],self.pred_dif_y,data_mu)
		S_yy = np.diag(GP_deriv.Pred_var_deriv(K_tr_tr_y,K_te_tr_der_y[1],K_te_tr_der_y[1].T,h_yy))

		return (V_xx,V_xy,V_yx,V_yy),(S_xx,S_xy,S_yx,S_yy)

	def time_order(self):
		"""
		Time order vector field inputs and outputs
		---
		Output: 
			loc_time_ind: reordered locations in the x y space according to time np.array(self.N_test_t,self.N_test_x*self.N_test_y,2)
			V_time_ind: reordered derivatives according to time np.array(self.N_test_t,self.N_test_x*self.N_test_y,2)
			S_time_ind: reordered variance of derivatives according to time np.array(self.N_test_t,self.N_test_x*self.N_test_y,2)
			V_2d: partial derivatives reordered according to time. np.array(self.N_test_t,self.N_test_x*self.N_test_y,4) 
				  order V_xx, V_xy, V_yx, V_yy
			S_2d: partial derivatives variance reordered according to time. np.array(self.N_test_t,self.N_test_x*self.N_test_y,4) 
				  order V_xx, V_xy, V_yx, V_yy
		"""

		loc_time_indexed=np.zeros((self.N_test_t,self.N_test_x*self.N_test_y,2))
		V_time_indexed=np.zeros((self.N_test_t,self.N_test_x*self.N_test_y,2))
		S_time_indexed=np.zeros((self.N_test_t,self.N_test_x*self.N_test_y,2))
		V_2d=np.zeros((self.N_test_t,self.N_test_x*self.N_test_y,4))
		S_2d=np.zeros((self.N_test_t,self.N_test_x*self.N_test_y,4))
		for n in range(self.N_test_t):
			index = np.where(self.x_3d_test[:,2] == self.t_test[n])[0]
			loc_time_indexed[n] = self.x_3d_test[index,0:2]
			V_time_indexed[n,:,0] = self.V[0,index] # x value derivative
			V_time_indexed[n,:,1] = self.V[1,index] # y value derivative
			S_time_indexed[n,:,0] = self.S[0,index] # x value derivative
			S_time_indexed[n,:,1] = self.S[1,index] # y value derivative
			
			V_2d[n,:,0] = self.V_xx[index,0]
			V_2d[n,:,1] = self.V_xy[index,0]
			V_2d[n,:,2] = self.V_yx[index,0]
			V_2d[n,:,3] = self.V_yy[index,0]
			S_2d[n,:,0] = self.S_xx[index]
			S_2d[n,:,1] = self.S_xy[index]
			S_2d[n,:,2] = self.S_yx[index]
			S_2d[n,:,3] = self.S_yy[index]
		
		return loc_time_indexed, V_time_indexed, V_2d, S_time_indexed, S_2d

	def calc_div_curl(self):
		"""
		Calculate values for div and curl
		---
		Output: 
			meshed div and curl and locations ready for contour plots in lists ordered by test times
		""" 

		for n in range(self.N_test_t):
			x_mesh = self.loc_time_indexed[n,:,0].reshape((self.N_test_x,self.N_test_y))
			y_mesh = self.loc_time_indexed[n,:,1].reshape((self.N_test_x,self.N_test_y))
			S_x_mesh = self.S_time_indexed[n,:,0].reshape((self.N_test_x,self.N_test_y))
			S_y_mesh = self.S_time_indexed[n,:,1].reshape((self.N_test_x,self.N_test_y))

			self.loc_time_ind_mesh_x.append(x_mesh)
			self.loc_time_ind_mesh_y.append(y_mesh)
			self.S_mesh_x.append(S_x_mesh)
			self.S_mesh_y.append(S_y_mesh)
			d = self.div[n,:].reshape((self.N_test_x,self.N_test_y))
			c = self.curl[n,:].reshape((self.N_test_x,self.N_test_y))
			d_var = self.div_var[n,:].reshape((self.N_test_x,self.N_test_y))
			c_var = self.curl_var[n,:].reshape((self.N_test_x,self.N_test_y))
			self.div_mesh.append(d)
			self.curl_mesh.append(c)
			self.div_var_mesh.append(d_var)
			self.curl_var_mesh.append(c_var)


class Multiple_traj(Pred_vector_field):
	"""
	Class takes in multiple trajectories in the same area on the map and combines them to allow the parent class to be used.
	"""

	def __init__(self, x_train, y_train, t_train, x_test, y_test, t_test, time_hyp=20,xy_hyp=1,var_hyp=1,noise_hyp=0.01):
		"""
		Need to define the constructor
		---
		Inputs:
			t_train: np.array of shape, tuple of N_traj arrays
			t_test: np.array of shape (M,1)
			y_train: np.array of shape, tuple of N_traj arrays
			x_train: np.array of shape, tuple of N_traj arrays
			x_test: np.array of shape (Ntx,1)
			y_test: np.array of shape (Nty,1)
		
		"""
		t_tot = time.time()

		######################
		self.N_traj = np.shape(t_train)[0]
		self.N_train = np.shape(t_train)[0]
		self.N_test = np.shape(t_test)[0]
		self.N_test_x = np.shape(x_test)[0]
		self.N_test_y = np.shape(y_test)[0]
		self.N_test_t = np.shape(t_test)[0]
		self.t_train = t_train
		self.t_test = t_test 
		self.x_train = x_train
		self.y_train = y_train
		self.x_test = x_test
		self.y_test = y_test
		self.time_hyp = time_hyp

		### Compute derivative observations for each trajectory
		pred_der_x_list = []
		pred_der_y_list = []
		mx_list = []
		my_list = []
		self.x_norm_list = []
		self.y_norm_list = []

		## Optimise the GPs
		t0 = time.time()
		######################

		for tr in range(self.N_traj):
			prior_mu = 0
			data_mu = np.zeros_like(self.x_train[tr])
			n_traj = self.x_train[tr].shape[0]

			x_train = self.x_train[tr]
			y_train = self.y_train[tr]
			self.x_norm_list.append(x_train)
			self.y_norm_list.append(y_train)
			## Reshaping ##
			t_tr = np.array(self.t_train[tr],dtype=np.float64).reshape(n_traj,1)
			x_tr = np.array(x_train,dtype=np.float64).reshape(n_traj,1)
			y_tr = np.array(y_train,dtype=np.float64).reshape(n_traj,1)

			# First Layer: Creating the model (fitting to trajectories for each agent)
			K_y = gpflow.kernels.RBF(1)
			m_y = gpflow.models.GPR(t_tr,y_tr,K_y)
			m_y.compile()
			opt_y = gpflow.train.ScipyOptimizer()
			opt_y.minimize(m_y)         

			K_x = gpflow.kernels.RBF(1)
			m_x = gpflow.models.GPR(t_tr,x_tr,K_x)
			m_x.compile()
			opt_x = gpflow.train.ScipyOptimizer()
			opt_x.minimize(m_x)             

			mx_list.append(m_x)
			my_list.append(m_y)

			## Derivative GPs
			# Using hyperparameters from gpflow independent
			# X
			hyp_x = [m_x.kern.lengthscales.value,np.sqrt(m_x.kern.variance.value),np.sqrt(m_x.likelihood.variance.value)]
			K_tr_tr_x = GP_deriv.SE_covariance(self.t_train[tr],self.t_train[tr],hyp_x[0],hyp_x[1]) + np.identity(n_traj)*hyp_x[2]**2
			# Second derivative
			K_te_tr_der_x = GP_deriv.SE_covariance_derivative_xx(self.t_train[tr],self.t_train[tr],hyp_x[0],hyp_x[1])
			pred_x_der = GP_deriv.Pred_mean(prior_mu,K_tr_tr_x,K_te_tr_der_x,x_train,data_mu)
			pred_der_x_list.append(pred_x_der)

			# Y
			hyp_y = [m_y.kern.lengthscales.value,np.sqrt(m_y.kern.variance.value),np.sqrt(m_y.likelihood.variance.value)]
			K_tr_tr_y = GP_deriv.SE_covariance(self.t_train[tr],self.t_train[tr],hyp_y[0],hyp_y[1]) + np.identity(n_traj)*hyp_y[2]**2
			# Second derivative
			K_te_tr_der_y = GP_deriv.SE_covariance_derivative_xx(self.t_train[tr],self.t_train[tr],hyp_y[0],hyp_y[1])
			pred_y_der = GP_deriv.Pred_mean(prior_mu,K_tr_tr_y,K_te_tr_der_y,y_train,data_mu)
			pred_der_y_list.append(pred_y_der)

		######################
		print('Wall time for first layer GPs: %0.2fs' % (time.time() - t0))

		# assign to plot derivatives of trajectories
		self.pred_der_x_list = pred_der_x_list
		self.pred_der_y_list = pred_der_y_list
		self.deriv_mx_list = mx_list
		self.deriv_my_list = my_list

		x = self.x_train[0]
		y = self.y_train[0]
		t = self.t_train[0]
		x_der = pred_der_x_list[0]
		y_der = pred_der_y_list[0]
		for tr in range(1,self.N_traj):
			x = np.hstack([x,self.x_train[tr]])
			y = np.hstack([y,self.y_train[tr]])
			t = np.hstack([t,self.t_train[tr]])
			x_der = np.hstack([x_der,pred_der_x_list[tr]])
			y_der = np.hstack([y_der,pred_der_y_list[tr]])

		### Second stage of predicting the vector field
		self.X,self.Y,self.T = np.meshgrid(self.x_test,self.y_test,self.t_test)
		y_t = self.Y.reshape((self.N_test_x*self.N_test_y*self.N_test_t))
		x_t = self.X.reshape((self.N_test_x*self.N_test_y*self.N_test_t))
		t_t = self.T.reshape((self.N_test_x*self.N_test_y*self.N_test_t))
		# Test points stacked 
		self.x_3d_test = np.vstack([x_t,y_t,t_t]).T

		# Training points
		self.N_train = x_der.shape[0]
		self.pred_dif_x = x_der.reshape((self.N_train,1))
		self.pred_dif_x = np.array(self.pred_dif_x,dtype=np.float64)
		self.pred_dif_y = y_der.reshape((self.N_train,1))
		self.pred_dif_y = np.array(self.pred_dif_y,dtype=np.float64)
		self.x_3d_train = np.vstack([x,y,t]).T
		self.x_3d_train = np.array(self.x_3d_train,dtype=np.float64)

		# 3D GPs
		t0 = time.time()
		####################################
		K_3d_x = gpflow.kernels.RBF(input_dim=3, ARD=True)
		self.m_3d_x = gpflow.models.GPR(self.x_3d_train,self.pred_dif_x,K_3d_x)
		self.m_3d_x.compile()
		opt_3dx = gpflow.train.ScipyOptimizer()
		opt_3dx.minimize(self.m_3d_x)           
		
		K_3d_y = gpflow.kernels.RBF(input_dim=3, ARD=True)
		self.m_3d_y = gpflow.models.GPR(self.x_3d_train,self.pred_dif_y,K_3d_y)
		self.m_3d_y.compile()
		opt_3dy = gpflow.train.ScipyOptimizer()
		opt_3dy.minimize(self.m_3d_y)           
		###################################
		print('Wall time for second layer GPs: %0.2fs' % (time.time() - t0))

		# Test derivative field:
		t0 = time.time()
		####################################
		self.V,self.S = self.vector_field()
		####################################
		print('Wall time for predictive vector field: %0.2fs' % (time.time() - t0))

		## Probabilistic derivatives of vector field
		t0 = time.time()
		####################################
		(self.V_xx,self.V_xy,self.V_yx,self.V_yy), (self.S_xx,self.S_xy,self.S_yx,self.S_yy) = self.first_derivative()
		####################################
		print('Wall time for second order derivatives: %0.2fs' % (time.time() - t0))

		# Time order test points
		t0 = time.time()
		####################################
		self.loc_time_indexed, self.V_time_indexed, self.V_2d, self.S_time_indexed, self.S_2d = self.time_order()
		####################################
		print('Wall time for rearranging the arrays to be time ordered: %0.2fs' % (time.time() - t0))

		self.V_xx_time_indexed = self.V_2d[:,:,0]
		self.V_xy_time_indexed = self.V_2d[:,:,1]
		self.V_yx_time_indexed = self.V_2d[:,:,2]
		self.V_yy_time_indexed = self.V_2d[:,:,3]

		self.S_xx_time_indexed = self.S_2d[:,:,0]
		self.S_xy_time_indexed = self.S_2d[:,:,1]
		self.S_yx_time_indexed = self.S_2d[:,:,2]
		self.S_yy_time_indexed = self.S_2d[:,:,3]

		# Calculate div and curl
		self.div = self.V_xx_time_indexed + self.V_yy_time_indexed
		self.div_var = self.S_xx_time_indexed + self.S_yy_time_indexed
		self.curl = self.V_yx_time_indexed - self.V_xy_time_indexed
		self.curl_var = self.S_yx_time_indexed + self.S_xy_time_indexed
		# Reorder for plotting
		self.loc_time_ind_mesh_x = []
		self.loc_time_ind_mesh_y = []
		self.div_mesh = []      
		self.div_var_mesh = []      
		self.curl_mesh = []     
		self.curl_var_mesh = []     
		self.S_mesh_x = []          
		self.S_mesh_y = []          


		t0 = time.time()
		####################################
		self.calc_div_curl()
		####################################
		print('Wall time for rearranging the arrays to be a mesh: %0.2fs' % (time.time() - t0))
		
		# KL-divergence at each point in the space:
		t0 = time.time()
		####################################

		# Prior: GP(0,K_te_te_der)
		hyp_x = [self.m_3d_x.kern.lengthscales.value,np.sqrt(self.m_3d_x.kern.variance.value),np.sqrt(self.m_3d_x.likelihood.variance.value)]
		hyp_y = [self.m_3d_y.kern.lengthscales.value,np.sqrt(self.m_3d_y.kern.variance.value),np.sqrt(self.m_3d_y.likelihood.variance.value)]
		prior = (hyp_x[1]/hyp_x[0][0])**2 + (hyp_y[1]/hyp_y[0][1])**2
		prior_var = prior * np.ones_like(self.div_mesh[0])

		N_0 = [np.zeros_like(self.div_mesh[0]),prior_var]

		self.KL_div_mesh = []
		self.KL_div_mesh_signed = []
		for t in range(self.N_test):
			# Posterior:
			N_1 = [self.div_mesh[t],self.div_var_mesh[t]]
			self.KL_div_mesh.append(KL_div(N_0,N_1))
			self.KL_div_mesh_signed.append(KL_div(N_0,N_1)*np.sign(N_1[0]))

		###################################
		print('Wall time for KL-divergence: %0.2fs' % (time.time() - t0))
		
		###################################
		print('Total wall time: %0.2fs' % (time.time() - t_tot))