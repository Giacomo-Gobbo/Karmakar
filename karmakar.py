import numpy as np
import numpy.linalg as LA

def add_costraint(self):
    pass

def is_feasible(A, b, c):
    
  
    return True




def karmarkar(A, b, c, x0, gamma):
    iter = 24
    v = np.zeros((iter+1, len(b)))
    x = np.zeros((iter+1, len(x0)))

    x[0] = x0+1    
    x[1] = x0    

    k = 1

    while k < iter and abs(np.sum(x[k-1] - x[k])) > 1e-16:
        print(x[k])
        v[k] = b-A@x[k]
        Dv = np.diag(v[k]) + np.diag(np.array([1e-16] * len(b)))
        
        try:
            Dvinvert = LA.inv(Dv.T@Dv)@Dv.T
        except:
            try:
                Dvinvert = LA.inv(Dv)
            except:           
                return "unbounded1", None

        Dvinvert = np.power(Dvinvert, 2)
        tmp = A.T@Dvinvert@A
        
        maximum = np.max(tmp)
        
        tmp /= maximum
        tmp += np.diag(np.array([1e-16] * len(tmp)))

        try:
            hx = LA.inv(tmp.T@tmp)@tmp.T@c
        except:
            try:
                hx = LA.inv(tmp)@c
            except:
                return "unbounded2", None
        
        hv = -A@hx

        if (hv >= 0).all():
            return "unbounded3", None
        
        alpha = gamma*np.min([-v[k, i] / hv[i] for i in range(len(hv)) if hv[i] < 0])
        x[k+1] = x[k] + alpha*hx

        k += 1
    
    return "bounded", x


if __name__ == "__main__":
    gamma = 0.5
    tipo = type(np.array([]))

# -----------------------------------------------------------------------
    c = -np.array([-20, -24])
    A = np.array([[3, 6], [4, 2]])
    b = np.array([60, 32])
    x0 = np.array([1, 1])

    print('----------------------------------------')
    bound, x = karmarkar(A, b, c, x0, gamma)
    print(bound)
    if type(x) == tipo:
        print(c@x[-1])

# -----------------------------------------------------------------------
    c = -np.array([-1, -1, 0])
    A = np.array([[1, -1, 0], [1, 1, 1]])
    b = np.array([0,2])
    x0 = np.array([2.1, 2.1, 1.8])

    print('----------------------------------------')
    bound, x = karmarkar(A, b, c, x0, gamma)
    print(bound)
    if type(x) == tipo:
        print(c@x[-1])

# -----------------------------------------------------------------------
    A = np.eye(2)
    b = np.array([2, 1])
    c = np.array([2, 2])
    x0 = np.array([1, 0.1])

    print('----------------------------------------')
    bound, x = karmarkar(A, b, c, x0, gamma)
    print(bound)
    if type(x) == tipo:
        print(c@x[-1])

# -----------------------------------------------------------------------
    A = np.array([[2, -3], [1, 5]])
    b = np.array([2, 1])
    c = np.array([2, 2])
    x0 = np.array([0.5, 0])

    print('----------------------------------------')
    bound, x = karmarkar(A, b, c, x0, gamma)
    print(bound)
    if type(x) == tipo:
        print(c@x[-1])

# -----------------------------------------------------------------------
    A = np.array([[1, -1], [1, 1]])
    b = np.array([2, 6])
    c = np.array([2, 2])
    x0 = np.array([0.5, 0])

    print('----------------------------------------')
    bound, x = karmarkar(A, b, c, x0, gamma)
    print(bound)
    if type(x) == tipo:
        print(c@x[-1])

# -----------------------------------------------------------------------
    A = np.array([[1, -1], [1, 1]])
    b = np.array([2, 6])
    c = np.array([1, 0.5])
    x0 = np.array([0.5, 0])

    print('----------------------------------------')
    bound, x = karmarkar(A, b, c, x0, gamma)
    print(bound)
    if type(x) == tipo:
        print(c@x[-1])

    # c = [-20, -24]
    # A = [[3, 6], [4, 2]]
    # b = [60, 32]
    # print(is_feasible(c, A, b, [], []))
    
    # c = [-1, -1, 0]
    # A = [[1, -1, 0], [1, 1, 1]]
    # b = [0,2]
    # print(is_feasible(c, A, b, [], []))
    
    # c = [2, 2]
    # A = [[1, -1], [1, 1]]
    # b = [2, 1]
    # print(is_feasible(c, A, b, [], []))
    
    # c = [2, 2]
    # A = [[1, -1], [1, 1]]
    # b = [2, 6]
    # print(is_feasible(c, A, b, [], []))
    
    # c = [1, 0.5]
    # A = [[1, -1], [1, 1]]
    # b = [2, 6]
    # print(is_feasible(c, A, b, [], []))

    
    # c = [1, -2]
    # A = [[1, 1], [1, 1]]
    # b = [-1, 6]
    # print(is_feasible(c, A, b, [], []))
    