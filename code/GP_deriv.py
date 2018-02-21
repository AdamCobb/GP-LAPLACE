import numpy as np
import scipy.optimize

def SE_covariance(x1, x2, hypInput, hypOutput):
	"""
	Calculates the covariance matrix between two vectors x1 and x2 according to the sqared exponential
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,1) and (M,1)
		hypInput: Input scale length (length scale)
		hypOutput: Output scale length
	---
	Output: 
		K = covariance shape (N,M)
	"""
	K = (hypOutput**2)*np.exp(np.power(np.subtract.outer(x1,x2),2)*-0.5/(hypInput**2))
	return K

def SE_covariance_derivative_x1(x1, x2, hypInput, hypOutput):
	"""
	Calculates the derivative of the squared exponential covariance matrix with respect to x1
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,1) and (M,1)
		hypInput: Input scale length (length scale)
		hypOutput: Output scale length
	---
	Output: 
		K = covariance shape (N,M)
	"""
	X = -hypInput**(-2) * np.subtract.outer(x1,x2)
	K = X * (hypOutput**2)*np.exp(np.power(np.subtract.outer(x1,x2),2)*-0.5/(hypInput**2))
	return K

def SE_covariance_derivative_x2(x1, x2, hypInput, shypOutput):
	"""
	Calculates the derivative of the squared exponential covariance matrix with respect to x1
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,1) and (M,1)
		hypInput: Input scale length (length scale)
		hypOutput: Output scale length
	---
	Output: 
		K = covariance shape (N,M)
	"""
	X = hypInput**(-2) * np.subtract.outer(x1,x2)
	K = X * (hypOutput**2)*np.exp(np.power(np.subtract.outer(x1,x2),2)*-0.5/(hypInput**2))
	return K

def SE_covariance_derivative_xx(x1, x2, hypInput, hypOutput):
	"""
	Calculates the second derivative of the squared exponential covariance matrix with respect to x1
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,1) and (M,1)
		hypInput: Input scale length (length scale)
		hypOutput: Output scale length
	---
	Output: 
		K = covariance shape (N,M)
	"""
	d = np.power(np.subtract.outer(x1,x2),2)
	X = hypInput**(-2) * (d*hypInput**(-2) - 1)
	K = X * (hypOutput**2)*np.exp(d*-0.5/(hypInput**2))
	return K

def SE_covariance_derivative_x1x2(x1, x2, hypInput, hypOutput):
	"""
	Calculates the second derivative of the squared exponential covariance matrix with respect to x1
	---
	Inputs:
		x1,x2: vector location in the plane np.array shape (N,1) and (M,1)
		hypInput: Input scale length (length scale)
		hypOutput: Output scale length
	---
	Output: 
		K = covariance shape (N,M)
	"""
	N = np.shape(x1)[0]; M = np.shape(x2)[0]
	X = np.zeros((N,M))
	for i in range(N):
		for j in range(M):
			X[i,j] = x1[i]*x2[j]
	d = np.power(np.subtract.outer(x1,x2),2)
	K = (hypInput**-4) * X * (hypOutput**2)*np.exp(d*-0.5/(hypInput**2))
	return K

def Pred_mean(priorMu, K_tr_tr, K_te_tr, y_train, dataMu):
	"""
	Calculates the predictive mean of the GP
	---
	Inputs:
		priorMu: mean function as a function of the test inputs shape(x_te)
		K_tr_tr: covariance matrix (N,N) for N training points
		K_te_tr: covariance matrix (M,N) for N training points and M test points 
		y_train: training points corresponding to N inputs
		dataMu: mean function as a function of training inputs
	---
	Output: 
		y_test = predictive mean shape (M,)
	"""
	#Inverse using Cholesky Decomp 
	# LL' = K, so K = inv(L')*inv(L)
	sigma = np.eye(K_tr_tr.shape[0])*0.
	flag = False
	while flag == False:
		try: 
			cholX = scipy.linalg.cho_solve(scipy.linalg.cho_factor(K_tr_tr+sigma),(y_train-dataMu))
			flag = True
		except np.linalg.LinAlgError:
			print('LinAlgError so add to diagonal')
			# Add to diagonal to try and make positive semi definite
			sigma = np.eye(K_tr_tr.shape[0])*0.1

	y_test = priorMu + np.dot(K_te_tr,cholX)
	return y_test

def Pred_var(K_te_te, K_tr_tr, K_te_tr):
	"""
	Calculates the predictive covariance of the GP
	---
	Inputs:
		K_te_te: covariance matrix (M,M) for M testing points 
		K_tr_tr: covariance matrix (N,N) for N training points
		K_te_tr: covariance matrix (M,N) for N training points and M test points 
	---
	Output: 
		K_pred = predictive cov shape (M,M)
	"""
	K_tr_te = K_te_tr.T
	#Inverse using Cholesky Decomp 
	# LL' = K, so K = inv(L')*inv(L) 
	cholX = scipy.linalg.cho_solve(scipy.linalg.cho_factor(K_tr_tr),(K_tr_te))
	
	K_pred = K_te_te - np.dot(K_te_tr,cholX)
	
	return K_pred

def Pred_var_deriv(K_tr_tr, K_te_tr, K_tr_te, hyp):
	"""
	Calculates the predictive covariance of the derivative of the GP
	---
	Inputs:
		K_tr_tr: covariance matrix (N,N) for N training points
		K_te_tr: covariance matrix (M,N) for N training points and M test points
		K_tr_te: covariance matrix (N,M) for N training points and M test points
		hyp: np.array of hyperparameters 
			hyp[0] = Input scale length
			hyp[1] = Output scale length
			hyp[2] = Noise hyperparameter
	---
	Output: 
		K_pred_deriv = predictive derivative cov shape (M,M)
	"""
	M = np.shape(K_te_tr)[0]
	A = np.eye(M) * (hyp[1]/hyp[0])**2
	#Inverse using Cholesky Decomp 
	# LL' = K, so K = inv(L')*inv(L) 
	cholX = scipy.linalg.cho_solve(scipy.linalg.cho_factor(K_tr_tr),(K_tr_te))
	
	K_pred_deriv = A - np.dot(K_te_tr,cholX)
	
	return K_pred_deriv

def Pred_var_deriv_xx(K_tr_tr, K_te_tr, K_tr_te, hyp):
	"""
	Calculates the predictive covariance of the second derivative of the GP
	---
	Inputs:
		K_tr_tr: covariance matrix (N,N) for N training points
		K_te_tr: covariance matrix (M,N) for N training points and M test points
		K_tr_te: covariance matrix (N,M) for N training points and M test points
		hyp: np.array of hyperparameters 
			hyp[0] = Input scale length
			hyp[1] = Output scale length
			hyp[2] = Noise hyperparameter
	---
	Output: 
		K_pred_deriv_xx = predictive derivative cov shape (M,M)
	"""
	M = np.shape(K_te_tr)[0]
	A = np.eye(M) * (3*hyp[1]**2)/(hyp[0]**4) 
	#Inverse using Cholesky Decomp 
	# LL' = K, so K = inv(L')*inv(L) 
	cholX = scipy.linalg.cho_solve(scipy.linalg.cho_factor(K_tr_tr),(K_tr_te))
	
	K_pred_deriv_xx = A - np.dot(K_te_tr,cholX)
	
	return K_pred_deriv_xx

def Loglik(K_tr_tr, y_train, dataMu):
	"""
	Calculate log likelihood for optimisation of hyperparameters
	---
	Inputs:
		K_tr_tr: covariance matrix (N,N) for N training points
		y_train: training points corresponding to N inputs
		dataMu: mean function as a function of training inputs
	---
	Output: 
		ll = log likelihood 
	"""
	cholX = scipy.linalg.cho_solve(scipy.linalg.cho_factor(K_tr_tr),(y_train-dataMu))
	(sign, logdet) = np.linalg.slogdet(K_tr_tr)
	ll =( -0.5 * np.dot(np.transpose(y_train-dataMu),cholX)) - 0.5*logdet - len(y_train) * 0.5 * np.log(2*np.pi)
	
	return ll

def optNLLfun(hyp, x_train, y_train, dataMu, verbose):
	"""
	Calculate negative log likelihood for optimisation of hyperparameters (This is the function to be optimised)
	---
	Inputs:
		x_train: input training points corresponding to N inputs
		y_train: training points corresponding to N inputs
		dataMu: mean function as a function of training inputs
		hyp: np.array of hyperparameters - hyp[0] = Input scale length, hyp[1] = Output scale length
										   hyp[2] = Noise hyperparameter
		verbose: Boolean - True if print NLL, False otherwise
	---
	Output: 
		ll = log likelihood 
	"""
	K_tr_tr = SE_covariance(x_train,x_train,hyp[0],hyp[1]) + np.identity(len(x_train))*hyp[2]**2
	NLL = - Loglik(K_tr_tr,y_train,dataMu)
	if verbose:    
		print ('Negative LL = ', NLL)    
	return NLL

def var_bounds(var, mean, std_multiple=1.96):
	"""
	Calculate the upper and lower bounds of the std 
	---
	Inputs:
		var: input covariance (M,M) for M training points
		std_multiple: float multiple of the standard deviations e.g. 1.96 == 95% confidence interval
		mean: mean function (M,)
	---
	Output: 
		conf_interval: upper row is upper bound and lower row is lower bound. np.array of size (2,M) 
	"""
	Var = np.diagonal(var)
	Sigma = np.sqrt(Var)
	ConfInterv = std_multiple*Sigma
	
	UpperConf = mean+ConfInterv
	LowerConf = mean-ConfInterv
	
	conf_interval = np.vstack([UpperConf,LowerConf])
	
	return conf_interval


class DerivGP:
	"""GP class for calculating derivatives """
	def __init__(self, t_train, t_test, y_train, hyp=[1,1,0.1]):
		"""
		Deriv GP.
		---
		Inputs:
			t_train: np.array of shape (N,)
			t_test:  np.array of shape (M,)
			y_train: np.array of shape (N,)
		"""

		self.t_train = t_train
		self.t_test = t_test
		self.y_train = y_train
		self.hyp = hyp
		# Prior means: (Zero)
		self.dataMu = np.zeros_like(t_train)
		self.priorMu = np.zeros_like(t_test)

		## Covariances:
		self.K_tr_tr = SE_covariance(self.t_train,self.t_train,self.hyp[0],self.hyp[1]) + np.identity(len(self.t_train))*self.hyp[2]**2
		self.K_te_tr = SE_covariance(self.t_test,self.t_train,self.hyp[0],self.hyp[1])
		self.K_te_te = SE_covariance(self.t_test,self.t_test,self.hyp[0],self.hyp[1])
		# Derivative of GP
		self.K_te_tr_der = SE_covariance_derivative_x1(self.t_test,self.t_train,self.hyp[0],self.hyp[1])
		self.K_tr_te_der = SE_covariance_derivative_x2(self.t_train,self.t_test,self.hyp[0],self.hyp[1])

		# Predictions
		self.pred_mean 		   = Pred_mean(self.priorMu,self.K_tr_tr,self.K_te_tr,self.y_train,self.dataMu)
		self.pred_deriv_mean   = Pred_mean(self.priorMu,self.K_tr_tr,self.K_te_tr_der,self.y_train,self.dataMu)
		self.var 			   = Pred_var(self.K_te_te, self.K_tr_tr, self.K_te_tr)
		self.pred_bounds       = var_bounds(var=self.var,mean=self.pred_mean,std_multiple=1.96)
		self.var_deriv 		   = Pred_var_deriv(self.K_tr_tr,self.K_te_tr_der,self.K_tr_te_der,self.hyp)
		self.pred_deriv_bounds = var_bounds(var=self.var_deriv,mean=self.pred_deriv_mean,std_multiple=1.96)
		self.num_deriv         = np.gradient(self.y_train) # Numerical derivative (2nd order accurate differences and 1st order at edges)

	def apply_deriv_GP(self):
		"""
		Calculates the derivaive GP given data  
		---
		Inputs:
			self: from the class
		---
		Output: 
			Updates the GP by optimising the hyperparameters and calculates the pred. mean and pred. variance for both the data
			and the derivative process.	    
		"""

		## Optimisation of GP-LL to get hyperparameters
		hyp_x = np.zeros((3,1))
		hyp_x[0] = self.hyp[0]
		hyp_x[1] = self.hyp[1]
		hyp_x[2] = self.hyp[2]

		#hyperparameters must be positive
		#Bound noise to prevent covariance from becoming singular
		#Limit length scale to avoid dividing by zero
		bnds =((0.,None),(0.,None),(0.,None))

		#Optimise but limit on 100
		hyp_opt_x = scipy.optimize.minimize(optNLLfun, hyp_x, method='L-BFGS-B',
						args = (self.t_train,self.y_train,self.dataMu,False), bounds=bnds, options={'maxiter':50})

		## Make covariance with new hyperparameters
		self.hyp[0] = hyp_opt_x['x'][0]
		self.hyp[1] = hyp_opt_x['x'][1]
		self.hyp[2] = hyp_opt_x['x'][2]
		# Add noise to covariance of training
		self.K_tr_tr = SE_covariance(self.t_train,self.t_train,self.hyp[0],self.hyp[1]) + np.identity(len(self.t_train))*self.hyp[2]**2
		self.K_te_tr = SE_covariance(self.t_test,self.t_train,self.hyp[0],self.hyp[1])
		self.K_te_te = SE_covariance(self.t_test,self.t_test,self.hyp[0],self.hyp[1])

		# Final predictions
		self.pred_mean = Pred_mean(self.priorMu,self.K_tr_tr,self.K_te_tr,self.y_train,self.dataMu)
		self.var = Pred_var(self.K_te_te, self.K_tr_tr, self.K_te_tr)
		self.pred_bounds = var_bounds(var=self.var,mean=self.pred_mean,std_multiple=1.96)

		# Derivative of GP
		self.K_te_tr_der = SE_covariance_derivative_x1(self.t_test,self.t_train,self.hyp[0],self.hyp[1])
		self.K_tr_te_der = SE_covariance_derivative_x2(self.t_train,self.t_test,self.hyp[0],self.hyp[1])
		# Final predictions of deriv of GP
		self.pred_deriv_mean = Pred_mean(self.priorMu,self.K_tr_tr,self.K_te_tr_der,self.y_train,self.dataMu)
		self.var_deriv = Pred_var_deriv(self.K_tr_tr,self.K_te_tr_der,self.K_tr_te_der,self.hyp)
		self.pred_deriv_bounds = var_bounds(var=self.var_deriv,mean=self.pred_deriv_mean,std_multiple=1.96)

	def predict(self,test):
		"""
		Calculates the derivative GP given data  
		---
		Inputs:
			self: from the class
			test: test points np.array(T,)
		---
		Output: 
			pred: predictive mean (T,)
		"""
		priorMu = np.zeros_like(test)
		K_te_tr = SE_covariance(test,self.t_train,self.hyp[0],self.hyp[1])
		return Pred_mean(priorMu,self.K_tr_tr,K_te_tr,self.y_train,self.dataMu)

# Inherit from the DerivGP
class DerivGP_2nd(DerivGP):
	def __init__(self,t_train,t_test,y_train,hyp=[1,1,0.1]):
		super().__init__(t_train, t_test, y_train, hyp)

		# second derivative of GP
		self.K_te_tr_der_xx = SE_covariance_derivative_xx(self.t_test,self.t_train,self.hyp[0],self.hyp[1])
		self.K_tr_te_der_xx = SE_covariance_derivative_xx(self.t_train,self.t_test,self.hyp[0],self.hyp[1])
		# Final predictions of deriv of GP
		self.pred_deriv_mean_xx = Pred_mean(self.priorMu,self.K_tr_tr,self.K_te_tr_der_xx,self.y_train,self.dataMu)
		self.var_deriv_xx = Pred_var_deriv_xx(self.K_tr_tr,self.K_te_tr_der_xx,self.K_tr_te_der_xx,self.hyp)
		self.pred_deriv_bounds_xx = var_bounds(var=self.var_deriv_xx,mean=self.pred_deriv_mean_xx,std_multiple=1.96)

	def update_deriv_GP_xx(self):
		"""
		Calculates the 2nd derivaive GP given data  (No optimisation)
		---
		Inputs:
			self: from the class
		---
		Output: 
			Given updated hyperparameters it calculates the 2nd derivative process
			Updates the GP by optimising the hyperparameters and calculates the pred. mean and pred. variance for both the data
			and the 2nd derivative process.	    
		"""

		# second derivative of GP
		self.K_te_tr_der_xx = SE_covariance_derivative_xx(self.t_test,self.t_train,self.hyp[0],self.hyp[1])
		self.K_tr_te_der_xx = SE_covariance_derivative_xx(self.t_train,self.t_test,self.hyp[0],self.hyp[1])
		# Final predictions of deriv of GP
		self.pred_deriv_mean_xx = Pred_mean(self.priorMu,self.K_tr_tr,self.K_te_tr_der_xx,self.y_train,self.dataMu)
		self.var_deriv_xx = Pred_var_deriv_xx(self.K_tr_tr,self.K_te_tr_der_xx,self.K_tr_te_der_xx,self.hyp)
		self.pred_deriv_bounds_xx = var_bounds(var=self.var_deriv_xx,mean=self.pred_deriv_mean_xx,std_multiple=1.96)