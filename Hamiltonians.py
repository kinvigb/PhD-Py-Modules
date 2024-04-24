import numpy as np



def oned_h(structure ,N , inpt=None, oupt=None,  A=(1/50), B=(1/93), Zin=(1/50), Zout=(1/50)):
    
    M = np.zeros((N+1,N+1), dtype = complex)
    for i, val in enumerate(structure[-1]):
        M[i,i+1] = A if val == 'A' else B
        M[i+1,i] = A if val == 'A' else B
        
    H = np.zeros((N+1,N+1), dtype = complex)
    for i in range(N):
        i=int(i)
        if i>0 and i<N-1: 
            H[i,i+1] = M[i,i+1]/(np.sqrt(M[i+1,i+2]+M[i+1,i])*np.sqrt(M[i,i+1]+M[i,i-1]))
            H[i+1,i] = H[i,i+1]
        elif i==0:
            H[i,i+1] = M[i,i+1]/(np.sqrt(M[i+1,i+2]+M[i+1,i])*np.sqrt(M[i,i+1]))
            H[i+1,i] = H[i,i+1]
        elif i==N-1:
            H[i,i+1] = M[i,i+1]/(np.sqrt(M[i,i+1])*np.sqrt(M[i+1,i]+M[i,i-1]))
            H[i+1,i] = H[i,i+1]
    
            
    if inpt and oupt != None:
        if inpt == 0:     
            H[inpt,inpt] = 1j * Zin / (np.sqrt(M[inpt,inpt+1]))**2
        elif inpt == N:
            H[inpt,inpt] = 1j * Zin / (np.sqrt(M[inpt,inpt-1]))**2
        elif inpt == N-1:
            H[inpt,inpt] = 1j * Zin / (np.sqrt(M[inpt,inpt-1]+M[inpt,inpt+1])*np.sqrt(M[inpt+1,inpt]))
        else:
            H[inpt,inpt] = 1j * Zin / (np.sqrt(M[inpt,inpt-1]+M[inpt,inpt+1])*np.sqrt(M[inpt+1,inpt+2]+M[inpt+1,inpt]))
            
            
        if oupt == 0:     
            H[oupt,oupt] = 1j * Zin / (np.sqrt(M[oupt,oupt+1]))**2
        elif oupt == N:
            H[oupt,oupt] = 1j * Zin / (np.sqrt(M[oupt,oupt-1]))**2
        elif oupt == N-1:
            H[oupt,oupt] = 1j * Zin / (np.sqrt(M[oupt,oupt-1]+M[oupt,oupt+1])*np.sqrt(M[oupt+1,oupt]))
        else:
            H[oupt,oupt] = 1j * Zin / (np.sqrt(M[oupt,oupt-1]+M[oupt,oupt+1])*np.sqrt(M[oupt+1,oupt+2]+M[oupt+1,oupt]))
    
        
    return M, H

def bdg_h(i, N, u, t, d):
    Hbdg = np.zeros((N, N))
    ipos = int((i/2))   #nth term
    ineg = int((i/2)-1) #n-1 term
    print(N)
    print(i)
    if i==0:
        sig_n = (abs(u) + abs(t[ipos]) + abs(d[ipos]))**(-1/2)
        sig_n1 = (abs(u) + abs(t[ipos]) + abs(t[ipos+1]) + abs(d[ipos])+ abs(d[ipos+1]))**(-1/2)
        sig_neg = 0
    elif i==N-4:
        sig_n = (abs(u) + abs(t[ineg]) + abs(t[ipos]) + abs(d[ineg]) + abs(d[ipos]))**(-1/2)
        sig_n1 = (abs(t[ipos]) + abs(d[ipos]))**(-1/2)
        sig_neg = (abs(u) + abs(t[ineg-1]) + abs(t[ineg]) + abs(d[ineg-1]) + abs(d[ineg]))**(-1/2)
        
    elif i==N-2:
        sig_n = (abs(u) + abs(t[ineg]) + abs(d[ineg]))**(-1/2)
        sig_n1 = 0
        sig_neg = (abs(u) + abs(t[ineg-1]) + abs(t[ineg]) + abs(d[ineg-1]) + abs(d[ineg]))**(-1/2)
    else:
        sig_n = (abs(u) + abs(t[ineg]) + abs(t[ipos]) + abs(d[ineg]) + abs(d[ipos]))**(-1/2)
        sig_n1 = (abs(u) + abs(t[ipos]) + abs(t[ipos+1]) + abs(d[ipos])+ abs(d[ipos+1]))**(-1/2)
        sig_neg = (abs(u) + abs(t[ineg-1]) + abs(t[ineg]) + abs(d[ineg-1]) + abs(d[ineg]))**(-1/2)
        
        
        
        
    Hbdg[i, i] = sig_n*u*sig_n
    Hbdg[i+1, i+1] = -1*sig_n*u*sig_n

    if i+3 <= N:
        Hbdg[i, i+2] = sig_n*t[ipos]*sig_n1          #first row 3rd column ...
        Hbdg[i+1, i+3] = -1*sig_n*t[ipos]*sig_n1

        Hbdg[i, i+3] = sig_n*d[ipos]*sig_n1
        Hbdg[i+1, i+2] = -1*sig_n*d[ipos]*sig_n1
    if i-2 >= 0:
        Hbdg[i, i-2] = sig_n*t[ineg]*sig_neg
        Hbdg[i+1, i-1] = -1*sig_n*t[ineg]*sig_neg

        Hbdg[i, i-1] = -1*sig_n*d[ineg]*sig_neg
        Hbdg[i+1, i-2] = sig_n*d[ineg]*sig_neg

    i += 2
    
    return Hbdg
