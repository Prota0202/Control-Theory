import package_DBR
import os

from importlib import reload
package_DBR = reload(package_DBR)

from package_DBR import *

#----------------------------------------------------------------------------------------------

def LeadLag_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method="EBD"):
    
    """
    The function "LEADLAG_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "FO_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]  slide 110 111
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+(Tlead/Ts))*MV[-1])-((Tlead/Ts)*MV[-2])))
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + K*Kp*(((Tlead/Ts)*MV[-1])+(((1-(Tlead/Ts))*MV[-2]))))           
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+(Tlead/Ts))*MV[-1])-((Tlead/Ts)*MV[-2])))
    else:
        PV.append(Kp*MV[-1])
        

    
#----------------------------------------------------------------------------------------------   

# PID functions

def Proportional_action(MVP, Kc, E, method):
    """
    Action proportionnel du PID
    MVP
    """
    if len(MVP)==0:
        MVP.append((Kc*E[-1]))
    else:
        if method == 'EBD-EBD':
            MVP.append(Kc*E[-1])
        elif method == 'TRAP-TRAP':  
            MVP.append(Kc*E[-1])        
        else: # meme méthode
            MVP.append(Kc*E[-1])
    return MVP

def Intergral_action(MVI, Kc, Ts, Ti, E, method):
    """
    Action intégrale du PID
    MVI
    """
    # MV[k] is MV[-1] and MV[k-1] is MV[-2]
    if len(MVI)==0:
        MVI.append(((Kc*Ts)/Ti)*(E[-1]))
    else:
        if method == 'EBD-EBD':
            MVI.append(MVI[-1] + ((Kc*Ts)/Ti)*(E[-1]))
        elif method == 'TRAP-TRAP':
            MVI.append(MVI[-1] + ((Kc*Ts)/(2*Ti))*(E[-1] + E[-2]))        
        else: # EBD
            MVI.append(MVI[-1] + ((Kc*Ts)/Ti)*(E[-1]))
    return MVI

def Derivative_action(MVD,Tfd, Ts, Kc, Td, E, method):
    """
    Action derivative du PID
    MVD
    """
    if len(MVD)==0:
        MVD.append(((Kc*Td)/(Tfd+Ts))*(E[-1]))
    else:
        if method == 'EBD-EBD':
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1] + ((Kc*Td)/(Tfd+Ts))*(E[-1] - E[-2])) 
        elif method =='TRAP-TRAP':
            MVD.append((((Tfd-(Ts/2))/(Tfd+(Ts/2)))*MVD[-1] + ((Kc*Td)/(Tfd+(Ts/2)))*(E[-1] - E[-2])))     
        else:  # EBD
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1] + ((Kc*Td)/(Tfd+Ts))*(E[-1] - E[-2]))
    return MVD

def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVmin, MVmax, MV, MVP, MVI, MVD, E, ManFF=False, PVinit=0, method="EBD-EBD"):
    """
    SP = setpoint vector
    PV = process value vector
    Man = manual controller mode vector : bool
    MVMan = Manual value for MV vector
    MVFF = feedforward vector
    Kc = controller gain
    Ti= integral time constant
    Td = derivative time constant
    alpha = Tfd = alpha*Td = where Tfd is derivative filter time constant

    MVMin = min value for MV
    MVMax = max value for MV

    MV = Maniplated value vector
    MVP = proportional part of MV vector
    MVI =  integral part of MV vector
    MVD = derivative part of MV vector
    E = control error vector

    ManFF = activated FF in manuel mode
    PVInit = initial value of PV
    Method : discreditisation value for PV
        EBD-EBD: Euler Backward difference
        TRAP-TRAP: Trapezoïdal method

    The function PID_RT appends new values to the vector MV, MVP, MVI, MVD based on PID algorithm, controller mode and FF
    """

    # Initialisation 
    Tfd = alpha*Td
    if len(PV) == 0:
        E.append(SP[-1] - PVinit)
    else: 
        E.append(SP[-1] - PV[-1])
    
    # Calcul des 3 parties
    MVP = Proportional_action(MVP, Kc, E, method="EBD-EBD")
    MVI = Intergral_action(MVI, Kc, Ts, Ti, E, method="EBD-EBD")
    MVD = Derivative_action(MVD,Tfd, Ts, Kc, Td, E, method="EBD-EBD")

    # Mode manuel + anti-wind up (integrator help)
    if Man[-1] == True:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] 
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]
            
    # Saturation
    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) > MVmax:
        MVI[-1] = MVmax - MVP[-1] - MVD[-1] - MVFF[-1]
    elif (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) < MVmin:
        MVI[-1] = MVmin - MVP[-1] - MVD[-1] - MVFF[-1]

    # Ajout sur MV
    MV.append(MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1])
    
    
    
#----------------------------------------------------------------------------------------------
    
def IMCTuning(K, Tlag1, Tlag2=0, theta=0, gamma = 0.5, process="FOPDT-PI"):
    """
    IMCTuning computes the IMC PID tuning parameters for FOPDT and SOPDT processes.
    K: process gain (Kp)
    Tlag1: first (main) lag time constant [s]
    Tlag2: second lag time constant [s]
    theta: delay [s]
    gamma: used to computed the desired closed loop time constant Tclp [s] (range [0.2 -> 0.9])
    process:
        FOPDT-PI: First Order Plus Dead Time for P-I control (IMC tuning case G)
        FOPDT-PID: First Order Plus Dead Time for P-I-D control (IMC tuning case H)
        SOPDT :Second Order Plus Dead Time for P-I-D control (IMC tuning case I)
        
    return : PID controller parameters Kc, Ti and Td
    """
    Tclp = gamma*Tlag1 
    if process=="FOPDT-PI":
        Kc = (Tlag1/(Tclp+theta))/K
        Ti = Tlag1
        Td = 0
    elif process=="FOPDT-PID":
        Kc= ((Tlag1 + theta/2)/(Tclp + theta/2))/K
        Ti = Tlag1 + theta/2
        Td = (Tlag1*theta)/(2*Tlag1+theta)
    elif process=="SOPDT": 
        Kc = ((Tlag1 + Tlag2)/(Tclp + theta))/K
        Ti = (Tlag1 +Tlag2)
        Td = ((Tlag1*Tlag2))/(Tlag1+Tlag2)
    else:
        Kc = (Tlag1/(Tclp+theta))/K
        Ti = Tlag1
        Td = 0
    return (Kc, Ti, Td)


#----------------------------------------------------------------------------------------------

def Margin(Ps,C,omega,Show=True):
    """
    Calcule la marge de gain et la marge de phase. Elles nous permettent d'analyser la robustesse du PID.
    :Ps : Process
    :C: Fonction de transfert du Contrôleur 
    :omega : vecteur de la fréquence 
    :show : autorise l'affichage graphique

    """
    # Initialisation des paramètres
    s = 1j*omega
    Kc = C.parameters['Kc']
    Ti = C.parameters['Ti']
    Td = C.parameters['Td']
    Tfd = C.parameters['Tfd']
    
    # Calcul du Controller 
    Cs = Kc*(1 + 1/(Ti*s)+ (Td*s)/(Tfd*s +1))
    
    # Loop gain L(s) = P(s)C(s)
    Ls =Cs*Ps 

    # Plot de L(s)
    if Show ==True:
        fig, (ax_freq, ax_time) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Amplitude
        ax_freq.semilogx(omega,20*np.log10(np.abs(Ls)),label='L(s)')
        gain_min = np.min(20*np.log10(np.abs(Ls)/5))
        gain_max = np.max(20*np.log10(np.abs(Ls)*5))
        ax_freq.set_xlim([np.min(omega), np.max(omega)])
        ax_freq.set_ylim([gain_min, gain_max])
        ax_freq.set_ylabel('Amplitude |P| [db]')
        ax_freq.set_title('Bode plot of P')
        ax_freq.legend(loc='best')
            
        # Phase
        ax_time.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)),label='L(s)')
        ax_time.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_time.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_time.set_ylabel(r'Phase $\angle P$ [°]')
        ax_time.legend(loc='best')
        ax_freq.axhline(y=0,color='black')
        ax_time.axhline(y=-180,color='black')

    # Crossover frequency
    i = 0
    for value in Ls:   # slide 69
        i+=1
        dB = 20*np.log10(np.abs(value))
        if dB < 0.05 and dB > -0.05:
            OmegaC =  omega[i-1]
            PhaseC = np.angle(value,deg=True)
            break        

    # Ultimate Frequency
    n = 0
    for value in Ls:
        n+=1
        deg = np.angle(value,deg=True)
        if deg < -179.5 and deg > -180.5:
            OmegaU = omega[n-1]
            u_freq = 20*np.log10(np.abs(value))
            break
    
    # Affichage graphique
    if Show ==True:
        ax_freq.plot([OmegaU,OmegaU],[0,u_freq]) 
        ax_time.plot([OmegaC,OmegaC],[PhaseC,-180])
    print('Gain margin :',-u_freq,'dB at the ultimate frequency :',OmegaU,'rad/s')
    print('Phase margin : ',PhaseC +180,'° at the crossover frequency :',OmegaC,'rad/s')
    
    # Save picture
    nameFile = 'Plots/Gain_phase_margin'

    if not os.path.exists('Plots'):
        os.makedirs('Plots')

    #plt.savefig(nameFile + '.png',transparent=True)
    #plt.savefig(nameFile + '.pdf',transparent=True)  


