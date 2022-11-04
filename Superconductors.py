import numpy as np
import matplotlib.pyplot as plt

def modulus(x):
	return np.sqrt(x.real**(2) + x.imag**(2))

x = np.arange(-1,1,0.01)

y = []
for i in x:
	y.append(complex(np.cos(i),np.sin(i)))

y = np.array(y)

def y_plot(x,t1,t2):
	return t1*x+t2

# plt.plot(x,modulus(y_plot(y,1,1)),label='t1=1')
# plt.plot(x,modulus(y_plot(y,2,1)),label='t1=2')
# plt.plot(x,modulus(y_plot(y,3,1)),label='t1=3')
# plt.plot(x,modulus(y_plot(y,4,1)),label='t1=4')
# plt.title('t2=1')
# plt.ylabel('\u03B5$_{1}$',fontsize=20)
# plt.xlabel('k/\u03C0',fontsize=20)
# plt.legend()
# plt.show()

# plt.plot(x,modulus(y_plot(y,1,2)),label='t1=1')
# plt.plot(x,modulus(y_plot(y,2,2)),label='t1=2')
# plt.plot(x,modulus(y_plot(y,3,2)),label='t1=3')
# plt.plot(x,modulus(y_plot(y,4,2)),label='t1=4')
# plt.title('t2=2')
# plt.ylabel('\u03B5$_{1}$',fontsize=20)
# plt.xlabel('k/\u03C0',fontsize=20)
# plt.legend()
# plt.show()

# t1 = 1
# t2 = 0.5

# print(modulus(y_plot(complex(np.cos(2*np.pi*0/4),np.sin(2*np.pi*0/4)),t1,t2)))
# print(modulus(y_plot(complex(np.cos(2*np.pi*1/4),np.sin(2*np.pi*1/4)),t1,t2)))
# print(modulus(y_plot(complex(np.cos(2*np.pi*2/4),np.sin(2*np.pi*2/4)),t1,t2)))
# print(modulus(y_plot(complex(np.cos(2*np.pi*3/4),np.sin(2*np.pi*3/4)),t1,t2)))

# print(-modulus(y_plot(complex(np.cos(2*np.pi*0/4),np.sin(2*np.pi*0/4)),t1,t2)))
# print(-modulus(y_plot(complex(np.cos(2*np.pi*1/4),np.sin(2*np.pi*1/4)),t1,t2)))
# print(-modulus(y_plot(complex(np.cos(2*np.pi*2/4),np.sin(2*np.pi*2/4)),t1,t2)))
# print(-modulus(y_plot(complex(np.cos(2*np.pi*3/4),np.sin(2*np.pi*3/4)),t1,t2)))


# def matrix_H(t1,t2,N,bc):
# 	H = np.zeros((N,N),dtype='float64')
# 	H.ravel()[1::(2*N+2)]=-t1
# 	H.ravel()[2+N::(2*N+2)]=-t2
# 	H.ravel()[N::2*N+2]=-t1
# 	H.ravel()[2*N+1::2*N+2]=-t2
# 	H[N-1,0]=-t2
# 	if bc == 'periodic':
# 		H[0,N-1]=-t2
# 	return H

# H_mat=matrix_H(t1,t2,8,'periodic')
# H_mat_1=matrix_H(t1,t2,8,'not periodic')

# print(H_mat)

# print(np.linalg.eigvals(H_mat))
# print(np.linalg.eigvals(H_mat_1))

# L = 65

# N = np.arange(2,L,1)

# y_plot_per = []
# for n in N:
# 	y_plot_per.append(np.linalg.eigvals(matrix_H(t1,t2,n,'periodic')))

# for i in range(L-2):
# 	for j in range(0,i+2,1):
# 		plt.plot(i+2,y_plot_per[i][j],'o',color='blue')

# y_plot_ope = []
# for n in N:
# 	y_plot_ope.append(np.linalg.eigvals(matrix_H(t1,t2,n,'not periodic')))

# for i in range(L-2):
# 	for j in range(0,i+2,1):
# 		plt.plot(i+2,y_plot_ope[i][j],'o',color='black')

# plt.title('Blue = Periodic, Black = Open')


def matrix_Kitaev(t,del_0,mu,N):
	H_t = np.zeros((N,N))
	H_del = np.zeros((N,N))
	H_mu = np.zeros((N,N))
	H_t.ravel()[1::(2*N+2)]=-t
	H_t.ravel()[N::2*N+2]=-t
	H_t.ravel()[2+N::(2*N+2)]=-t
	H_t.ravel()[2*N+1::2*N+2]=-t
	H_del.ravel()[1::(2*N+2)]=del_0
	H_del.ravel()[N::2*N+2]=-del_0
	H_del.ravel()[2+N::(2*N+2)]=del_0
	H_del.ravel()[2*N+1::2*N+2]=-del_0
	for i in range(N):
		H_mu[i,i] = -mu
	C = H_t+H_mu
	S = H_del
	S_dag = S.conj().T
	H = np.block([[C,S],
		          [S_dag,-C]])
	return H

t = 1
del_0 = 2
mu = 1

# print(matrix_Kitaev(t,del_0,mu,8))

# print(np.linalg.eigvals(matrix_Kitaev(t,del_0,mu,8)))

L = 24

N = np.arange(2,L,1)

y_plot = []
for n in N:
	y_plot.append(np.linalg.eigvals(matrix_Kitaev(t,del_0,mu,n)))

y_plot = np.array(y_plot)

print(y_plot[0])

for i in range(L-2):
	for j in range(0,2*(i+2),1):
		plt.plot(i+2,y_plot[i][j],'o',color='blue')

plt.ylabel('E')
plt.xlabel('N')
plt.title('t='+str(t)+', \u0394$_{0}$='+str(del_0)+', \u03BC='+str(mu))
plt.show()