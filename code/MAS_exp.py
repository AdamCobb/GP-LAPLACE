import numpy as np
from scipy.stats import multivariate_normal

class Agent:
	""" Agent class that moves according to a vector field """
	def __init__(self, args):
		self.position        = args.get('position', np.array([.0,.0]))     # starting position of agent
		self.velocity        = args.get('velocity', np.array([.0,.0]))     # starting velocity of agent
		self.acceleration    = args.get('acceleration', np.array([.0,.0])) # starting acceleration of agent
		self.derivative_func = args.get('der_func', lambda x: x)           # vector field function
		self.t_step          = args.get('time_step', 0.1)    			     # time step (increment)
		self.noise_var       = args.get('noise_var', 0.1)    				 # noise model for agent
	    
	def dyn_sys_1st(self):
		""" Step through to the next time step in the first order dynamical system dynamical system """
		self.velocity = self.derivative_func(self.position).reshape((2,)) + np.random.normal(0,self.noise_var,2)
		self.position = self.position + self.velocity * self.t_step
	    
	def dyn_sys_2nd(self):
		""" Step through to the next time step in the second order dynamical system dynamical system """
		self.acceleration = self.derivative_func(self.position).reshape((2,))
		self.velocity = self.velocity + self.acceleration * self.t_step
		self.position = self.position + self.velocity * self.t_step
	    
	def get_position(self):
		""" Returns current agent position """
		return self.position

	def get_velocity(self):
		""" Returns current agent velocity """
		return self.velocity

	def get_acceleration(self):
		""" Returns current agent acceleration """ 
		return self.acceleration

	def update_derivative_func(self,new_func):
		""" Updates the derivative function """
		self.derivative_func = new_func

class ContinuousEnvironment:
	""" Potential field in which agents operate """

	def __init__(self, N=200, M=10, D=2, size=3):
		self.D = D # Dimension of world
		self.N = N # Time resolution
		self.M = M # Space resolution

		# Create potential field
		self.t_train, self.t_step = np.linspace(0,2*np.pi,self.N,retstep=True)
		self.x_test			      = np.linspace(-size,size,self.M) # Test point in x direction
		self.y_test			  	  = np.linspace(-size,size,self.M) # Test point in x direction
		self.X_test, self.Y_test  = np.meshgrid(self.x_test,self.y_test) # Grid of testpoints
		self.scale                = 100 # Potential scale. Strength 

		# Frames		
		self.utility_frames       = []
		self.vector_field_frames  = []
		self.div_frames           = []
		self.curl_frames          = []

		# Agents
		self.agents = []
	   
	def initialise_agents(self, agent_positions):
		""" Creates an agent at each position given """
		self.agent_starting_positions = agent_positions # For resetting the environment

		# Create agents at given positions
		for pos in agent_positions:
		    agent_args = {'position':np.array(pos), 'velocity':np.array([0.,0.]) ,
		        'acceleration':np.array([0.,0.]) , 'der_func': lambda x: x , 'time_step': self.t_step, 'noise_var': 0.01}
		    agent = Agent(agent_args)
		    self.agents.append(agent)

		# Store list of each agent movements, velocities, and accelerations
		self.n_agents         = len(self.agents)
		self.agent_pos 		  = np.zeros((self.N,2,self.n_agents)) 
		self.agent_vel        = np.zeros((self.N,2,self.n_agents)) 
		self.agent_acc        = np.zeros((self.N,2,self.n_agents))

	def reset(self):
		""" Resets all agents to their starting positions and clears all frames """
		self.utility_frames        = []
		self.vector_field_frames   = []
		self.div_frames            = []
		self.curl_frames           = []			

		self.agent_pos 		  = np.zeros((self.N,2,self.n_agents))
		self.agent_vel        = np.zeros((self.N,2,self.n_agents))
		self.agent_acc        = np.zeros((self.N,2,self.n_agents))	

		for i, agent in enumerate(self.agents):
			agent.position = np.array(self.agent_starting_positions[i])
			agent.velocity = np.array([0., 0,])
			agent.acceleration = np.array([0., 0.])

	def calculate_trajectories(self, means, covs):
		""" Steps through the environment, creating trajectories of length N """
		for t in range(self.N):
			self.step(t, means, covs)

	def Gaussian_deriv(self, mean, covariance, t):
		"""
	    Produce a time dependent derivative function of a Gaussian
	    	Derivatives -p(x)∑^-1 (x-mu)
		---
		Inputs:
		    t: time
		    mean: function of t 
		    covariance: funcion of t
		---
		Output: 
		    der: function that returns derivatives (returns (D,))
		"""
		cov = covariance(t)
		mu  = mean(t)

		rv   = multivariate_normal(mu.reshape((self.D,)), cov)
		prec = np.linalg.inv(cov)

		# Derivative function
		der = lambda x: (- rv.pdf(x) * np.dot(prec,x.reshape((self.D,1))-mu.reshape((self.D,1)))) * self.scale
		return der

	def Gaussian_2_deriv(self, mean, covariance, t):
		"""
		Produce a time dependent 2nd derivative function of a Gaussian
			Hessian p(x)(∑^-1 (x-mu)(x-mu)^T ∑^-1 - ∑-1)
		---
		Inputs:
		    t: time
		    mean: function of t 
		    covariance: funcion of t
		---
		Output: 
		    der_2: function that returns derivatives (returns (D,))
		"""
		cov = covariance(t)
		mu  = mean(t)

		rv   = multivariate_normal(mu.reshape((self.D,)), cov)
		prec = np.linalg.inv(cov)

		# Derivative function
		der2 = lambda x: (rv.pdf(x) * (np.dot(prec,np.dot(x.reshape((self.D,1))-mu,np.dot((x.reshape((self.D,1))-mu).T,prec)))-prec)) * self.scale   
		return der2

	def step(self, t, mean_list, cov_list):
		"""
		Step through the environment for all agents
		---
		Inputs:
		    t: index of t_train
		    mean_list: list of mean functions for mulivariate gaussian
		    cov_list: list of covariance functions for mulivariate gaussian  
		"""
		time = self.t_train[t]
		
		c = lambda x: np.sum([self.Gaussian_deriv(mu,cov,time)(x) for mu, cov in zip(mean_list, cov_list)],axis=0) 

		for n,agent in enumerate(self.agents):
			agent.update_derivative_func(c)
			agent.dyn_sys_2nd()
			self.agent_pos[t,:,n] = agent.get_position()
			self.agent_vel[t,:,n] = agent.get_velocity()
			self.agent_acc[t,:,n] = agent.get_acceleration()

	def potential(self, mean_list, cov_list):
		"""
		Calculate true div grad and curl.
		---
		Inputs:
		    mean_list: list of mean functions for mulivariate gaussian
		    cov_list: list of covariance functions for mulivariate gaussian  
		"""
		for t in self.t_train:
			util = np.empty(self.X_test.shape)
			grad_util = np.empty(self.X_test.shape + (2,))
			div_grad_util = np.empty(self.X_test.shape)
			curl_grad_util = np.empty(self.X_test.shape)
			for i in range(self.M):
				for j in range(self.M):
					x_ij = np.array([self.X_test[i,j],self.Y_test[i,j]])
					
					potential   = lambda x: np.sum([multivariate_normal.pdf(x, mu(t).reshape((self.D,)), cov(t)) for mu, cov in zip(mean_list, cov_list)],axis=0) * self.scale
					potent_der  = lambda x: np.sum([self.Gaussian_deriv(mu,cov,t)(x) for mu, cov in zip(mean_list, cov_list)],axis=0)
					potent_2der = lambda x: np.sum([self.Gaussian_2_deriv(mu,cov,t)(x) for mu, cov in zip(mean_list, cov_list)],axis=0)

					util[i,j] = potential(x_ij)
					grad_util[i,j,:] = potent_der(x_ij).reshape(2,)#-(d1 + d2).reshape((2,))
					div_grad_util[i,j] = potent_2der(x_ij)[0,0] + potent_2der(x_ij)[1,1]
					curl_grad_util[i,j] = potent_2der(x_ij)[0,1] - potent_2der(x_ij)[1,0]

			self.utility_frames.append(util)
			self.vector_field_frames.append(grad_util)
			self.div_frames.append(div_grad_util)
			self.curl_frames.append(curl_grad_util)    
        

	def get_trajectories(self):
		""" Creates and returns a copy of each agent's trajectory """
		t_test = np.copy(self.t_train[0::1])
		x_train = []
		y_train = []
		t_t = []
		for a in range(self.n_agents):
		    x_train.append(self.agent_pos[:,0,a])
		    y_train.append(self.agent_pos[:,1,a])
		    t_t.append(self.t_train)
		    
		t_train = t_t
		x_test = np.copy(self.x_test)
		y_test = np.copy(self.y_test)

		return x_train, y_train, t_train, x_test, y_test, t_test