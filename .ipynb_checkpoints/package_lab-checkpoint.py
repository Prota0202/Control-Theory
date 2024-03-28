import math
import numpy as np
from package_DBR import Process, Bode
import matplotlib.pyplot as plt


def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD-EBD'):
    '''
The function "PID_RT" needs to be included in a "for or while loop". 

:SP: SP (or SetPoint) vector 
:PV: PV (or Process Value) vector 
:Man: Man (or Manual controller mode) vector (True or False) 
:MVMan: MVMan (or Manual value for MV) vector 
:MVFF: MVFF (or Feedforward) vector 

:Kc: controller gain 
:Ti: integral time constant [s] 
:Td: derivative time constant [s] 
:alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s] 
:Ts: sampling period [s] 

:MVMin: minimum value for MV (used for saturation and anti wind-up) 
:MVMax: maximum value for MV (used for saturation and anti wind-up) 

:MV: MV (or Manipulated Value) vector 
:MVP: MVP (or Propotional part of MV) vector 
:MVI: MVI (or Integral part of MV) vector 
:MVD: MVD (or Derivative part of MV) vector 
:E: E (or control Error) vector 

:ManFF: Activated FF in manual mode (optional: default boolean value is False) 
:PVInit: Initial value for PV (optional: default value is 0): used if PID_RT is ran first in the squence and no value of PV is available yet. 

:method: discretisation method (optional: default value is 'EBD') 
    EBD-EBD: EBD for integral action and EBD for derivative action 
    EBD-TRAP: EBD for integral action and TRAP for derivative action 
    TRAP-EBD: TRAP for integral action and EBD for derivative action 
    TRAP-TRAP: TRAP for integral action and TRAP for derivative action 

The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD". The appended values are based on the PID algorithm, the controller mode, and feedforward. Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up. 
    '''
    if len(PV) == 0:
        E.append(SP[-1] - PVInit)
    else:
        E.append(SP[-1] - PV[-1])

    methodI, methodD = method.split('-')

    ''' Initialisation de MVI'''
    if len(MVI) == 0:
        MVI.append((Kc*Ts/Ti)*E[-1])
    else:
        if methodI == 'TRAP':
            MVI.append(MVI[-1] + (0.5*Kc*Ts/Ti)*(E[-1]+E[-2]))
        else:
            MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])

    '''Initialisation de MVD : voir slide 193  '''
    Tfd = alpha * Td
    if Td > 0:
        if len(MVD) != 0:
            if len(E) == 1:
                MVD.append((Tfd / (Tfd + Ts)) *
                           MVD[-1] + ((Kc * Td) / (Tfd + Ts)) * (E[-1]))
            else:
                MVD.append((Tfd / (Tfd + Ts)) *
                           MVD[-1] + ((Kc * Td) / (Tfd + Ts)) * (E[-1] - E[-2]))
        else:
            if len(E) == 1:
                MVD.append((Kc * Td) / (Tfd + Ts) * (E[-1]))
            else:
                MVD.append((Kc * Td) / (Tfd + Ts) * (E[-1] - E[-2]))

    '''Actualisation de MVP'''
    MVP.append(E[-1] * Kc)

    '''Activation Feedforward'''
    if ManFF:
        MVFFI = MVFF[-1]
    else:
        MVFFI = 0
    '''Mode manuel et anti-wind-up'''

    if Man[-1]:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFFI

    '''Limitation de MV'''

    MV_TEMP = MVP[-1] + MVI[-1] + MVD[-1] + MVFFI

    if MV_TEMP >= MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFFI
        MV_TEMP = MVMax

    if MV_TEMP <= MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFFI
        MV_TEMP = MVMin

    MV.append(MV_TEMP)


def LeadLag_RT(MV, Kp, T_lead, T_lag, Ts, PV, PVInit=0, method='EBD'):
    """
    The function "LeadLag_RT" needs to be included in a "for or while loop".

    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method

    The function "LeadLag_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """
    if (T_lag != 0):
        K = Ts/T_lag
        if len(PV) == 0:
            PV.append(PVInit)
        else:  # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K)) * PV[-1] + ((K*Kp)/(1+K)) *
                          ((1+(T_lead/Ts)) * MV[-1] - (T_lead/Ts) * MV[-2]))
            elif method == 'EFD':
                PV.append((1-K) * PV[-1] + K * Kp * ((T_lead/Ts)
                          * MV[-1] + (1 - (T_lead/Ts)) * MV[-2]))
            elif method == 'TRAP':
                PV.append((1/(2 * T_lag+Ts)) * ((2 * T_lag-Ts) *
                          PV[-1]+(2 * T_lead+Ts) * Kp * MV[-1] + (Ts-2*T_lead) * Kp * MV[-2]))
            else:
                PV.append((1/(1+K)) * PV[-1] + ((K*Kp)/(1+K)) *
                          ((1+(T_lead/Ts)) * MV[-1] - (T_lead/Ts) * MV[-2]))
    else:
        PV.append(Kp*MV[-1])


def IMCTuning(Kp, Tlag1, Tlag2=0, theta=0, gamma=0, process="FOPDT", model="classic", Tg=0, Tu=0, a=0, t1=0, t2=0):
    """
    The function "imc_tuning" is only for first and second order systems.
    :Kp: process gain
    :Tlag1: first (or main) lag time constant [s] used in your process
    :Tlag2: second lag time constant [s] used in your process
    :theta: delay [s] used in your process
    :gamma : constant used to get the closed loop time constant
    :process: process order (ex : FOPDT first order system wuth delay)
    :model: broida_simple or broida_complex for FOPDT
    :Tg:
    :Tu:
    :a:
    :t1: time for 28% of PV (100% being the steady state)
    :t2: time for 44% of PV
    :return: imc tuning parameters respectively: 
        - Kc: controller gain
        - Ti: reset time 
        - Td: derivative time       
    The function "imc_tuning" returns the parameteres that you will use in your PID depending on your process parameters
    """

    if (process == "FOPDT"):
        if (model == "broida_simple"):
            Tlag1 = Tg
            theta = Tu
        elif (model == "broida_complex"):
            Tlag1 = 5.5*(t2 - t1)
            theta = (2.8*t1) - (1.8*t2)

        Tc = gamma * Tlag1
        Kc = ((Tlag1 + theta/2) / (Tc + theta/2)) / Kp
        Ti = Tlag1 + theta/2
        Td = (Tlag1*theta) / (2*Tlag1 + theta)

    elif (process == "SOPDT"):
        if (model == "vdG"):
            Tlag1 = Tg * ((3*a*math.exp(1) - 1) / (1 + a*math.exp(1)))
            Tlag2 = Tg * ((1 - a*math.exp(1)) / (1 + a*math.exp(1)))
            theta = Tu - ((Tlag1*Tlag2) / (Tlag1 + 3*Tlag2))

        Tc = gamma * Tlag1
        Kc = ((Tlag1 + Tlag2) / (Tc + theta)) / Kp
        Ti = Tlag1 + Tlag2
        Td = (Tlag1*Tlag2) / (Tlag1 + Tlag2)

    print(f"Kp : {Kp}, Tlag1 : {Tlag1}, Tlag2 : {Tlag2}, thetha : {theta}, Kc : {Kc}, Ti : {Ti}, Td : {Td}")

#     else :
#         an = Tu / Tg
#         table = [ [0.0, 1.0],  [0.10, 2.72], [0.22, 3.69],[0.32, 4.46],[0.41, 5.12],[0.49, 5.70],[0.57, 6.23] ]

#         for i in range(len(table)):
#             if ( table[i][0] <= an < table[i+1][0]):
#                 n = i + 1
#                 bn = table[i][1]

#         Tlag1 = Tg / bn
#         Tuth = an * Tg
#         theta = Tu - Tuth

    return Kc, Ti, Td


class PID:
    def __init__(self, parameters):
        self.parameters = parameters
        self.parameters['Kc'] = float(parameters['Kc']) if 'Kc' in parameters else 1.0
        self.parameters['alpha'] = float(parameters['alpha']) if 'alpha' in parameters else 0.0
        self.parameters['Ti'] = float(parameters['Ti']) if 'Ti' in parameters else 0.0
        self.parameters['Td'] = float(parameters['Td']) if 'Td' in parameters else 0.0



def StabilityMargins(P: Process, C: PID, omega):
    """
    The function "stability_margins" needs to have 2 processes object in paramaters.

    :P: the system process
    :C: the controller 
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".

    The function "stability_margins" generates the bodes plots of the Loop gain and gives the phase and gain margins.
    """
    s = 1j*omega

    Ptheta = np.exp(-P.parameters['theta']*s)
    PGain = P.parameters['Kp']*np.ones_like(Ptheta)
    PLag1 = 1/(P.parameters['Tlag1']*s + 1)
    PLag2 = 1/(P.parameters['Tlag2']*s + 1)
    PLead1 = P.parameters['Tlead1']*s + 1
    PLead2 = P.parameters['Tlead2']*s + 1

    Ps = np.multiply(Ptheta, PGain)
    Ps = np.multiply(Ps, PLag1)
    Ps = np.multiply(Ps, PLag2)
    Ps = np.multiply(Ps, PLead1)
    Ps = np.multiply(Ps, PLead2)

    integration_action = 1 / (C.parameters['Ti']*s)
    Tfd = C.parameters['Td']*C.parameters['alpha']
    derivative_action = C.parameters['Td']*s / (1 + Tfd * s)
    Cs = np.multiply(C.parameters['Kc'],
                     (1 + integration_action + derivative_action))

    Ls = np.multiply(Ps, Cs)

    gain_values = 20*np.log10(np.abs(Ls))
    phase_values = (180/np.pi)*np.unwrap(np.angle(Ls))

    for i in range(len(gain_values)):
        if gain_values[i] > 0 and gain_values[i + 1] < 0:
            x_gain = i
            y_gain = gain_values[i]
            Wc = round(omega[i], 2)
            phase_margin = round(abs(phase_values[i] + 180), 2)
            break
    
    for i in range(len(phase_values)):
        if phase_values[i] > -180 and phase_values[i + 1] < -180:
            x_phase = i
            y_phase = phase_values[i]
            Wu = round(omega[i], 2)
            gain_margin = round(abs(gain_values[i]), 2)
            break


    fig, (ax_gain, ax_phase) = plt.subplots(2, 1)
    fig.set_figheight(12)
    fig.set_figwidth(22)

    # Gain part
    ax_gain.axhline(y=0, color='b', linestyle='-')
    ax_gain.semilogx(omega, 20*np.log10(np.abs(Ls)), label='L(s)')
    gain_min = np.min(20*np.log10(np.abs(Ls)/5))
    gain_max = np.max(20*np.log10(np.abs(Ls)*5))
    ax_gain.vlines(omega[x_phase], gain_min,
                   gain_values[x_phase], color='r', linestyle='--')
    ax_gain.vlines(omega[x_gain], gain_min, 0, color='blue', linestyle='--')

    #ax_gain.vlines(omega[x_phase],gain_values[x_phase],0, color = 'g')

    ax_gain.annotate(
        '', xy=(omega[x_phase], gain_values[x_phase]), xycoords='data',
        xytext=(omega[x_phase], 0), textcoords='data',
        arrowprops={'arrowstyle': '<->'})

    ax_gain.annotate(
        f"Gain margin = {gain_margin}dB at {Wu}rad/s ", xy=(omega[x_phase], gain_values[x_phase] / 2), xycoords='data',
        xytext=(5, 0), textcoords='offset points')

    ax_gain.set_xlim([np.min(omega), np.max(omega)])
    ax_gain.set_ylim([gain_min, gain_max])
    ax_gain.set_ylabel('Amplitude |L| [db]')
    ax_gain.set_title('Bode plot of L : ' +
                      f"Gain margin = {gain_margin}dB at {Wu}rad/s " + f"Phase margin = {phase_margin}° at {Wc}rad/s ")
    ax_gain.legend(loc='best')

    # Phase part
    ax_phase.axhline(y=-180, color='b', linestyle='-')
    ax_gain.vlines(y_gain, -180, y_phase)
    ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)), label='L(s)')
    ax_phase.set_xlim([np.min(omega), np.max(omega)])
    ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ls))) - 10
    ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ls))) + 10
    ax_phase.vlines(omega[x_gain], phase_values[x_gain],
                    ph_max, color='blue', linestyle='--')
    ax_phase.vlines(omega[x_phase], -180, ph_max, color='r', linestyle='--')
    #ax_phase.vlines(omega[x_gain],-180 ,phase_values[x_gain], color = 'g')

    ax_phase.annotate(
        '', xy=(omega[x_gain], -180), xycoords='data',
        xytext=(omega[x_gain], phase_values[x_gain]), textcoords='data',
        arrowprops={'arrowstyle': '<->'})

    ax_phase.annotate(
        f"Phase margin = {phase_margin}° at {Wc}rad/s ", xy=(omega[x_gain], (-180 + phase_values[x_gain]) / 2), xycoords='data',
        xytext=(5, 0), textcoords='offset points')

    ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
    ax_phase.set_ylabel(r'Phase $\angle L$ [°]')
    ax_phase.legend(loc='best')