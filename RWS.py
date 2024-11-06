import numpy as np
import pywt
from scipy import stats
# Class that defines the Extended Besov prior with all its functionalities
class Besov_Bernoulli_prior:
    def __init__(self,alpha,beta,theta,gamma,dist,wavelet,J,nu=0,mu=0):
        self.alpha = alpha
        self.beta = beta
        self.theta = theta
        self.gamma = gamma
        self.dist = dist
        self.wavelet = wavelet
        self.J = J
        self.nu = nu
        self.mu = mu
        self.L1 = 0
        self.L2 = 0
        self.phi = 0
        self.psi = 0
        self.initialize_wavelet_basis()

    
    def initialize_wavelet_basis(self):
        #Initializing the wavelet basis using wavefun function and making/storing a linear interpolator
        [phi, psi, x] = pywt.Wavelet(self.wavelet).wavefun(20)
        self.L1 = x[0]
        self.L2 = x[-1]
        self.phi = lambda y:  np.interp(y, x , phi, left=0, right=0)
        self.psi = lambda y:  np.interp(y, x , psi, left=0, right=0)

    def scale_weight(self,k):
        #Scaling function weight depending of the sumation index k
        return (1+np.abs(k))**self.gamma
    
    def wavelet_weight(self,j,k):
        # wavelet function weight given summation indicies j and k
        return 2**(self.alpha*j)*(j+1)**self.theta*(1+np.abs(k)/2**j)**self.beta
    
    def k_range(self,a,b,j):
        # Range of translations leading to an overlap between the support of the wavelet and the selected domain
        return np.arange(np.ceil(2**j*a-self.L2),np.floor(2**j*b-self.L1))
    

    def sample(self,x):
        #Computing a sample of the truncated random wavelet expansion
        k_range = self.k_range(x[0],x[-1],0)
        sample = np.zeros(len(x))
        for i in range(len(k_range)):
            Prob_Bernoulli = 2**(0*self.mu)*(1+np.abs(k_range[i])/2**0)**self.nu
            sample+= self.scale_weight(k_range[i])*stats.bernoulli.rvs(Prob_Bernoulli)*self.dist.rvs()*self.phi(x-k_range[i])

        for j in  range(self.J):
            k_range = self.k_range(x[0],x[-1],j)
            for i in range(len(k_range)):
                Prob_Bernoulli = 2**(j*self.mu)*(1+np.abs(k_range[i])/2**j)**self.nu
                sample+= self.wavelet_weight(j,k_range[i])*stats.bernoulli.rvs(Prob_Bernoulli)*self.dist.rvs()*2**(j/2)*self.psi(2**j*x-k_range[i])

        return sample  