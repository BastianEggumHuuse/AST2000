import numpy as np
import matplotlib.pyplot as plt

def LorentzFactor(v):
    return 1/(1-v**2)**0.5

def ToMeters(seconds):
    return seconds * 3 * 10**8

def ToSeconds(meters):
    return meters / (3 * 10**8)

l_y = 9.46 * 10**15

def Part2():

    v = 0.99
    l = LorentzFactor(v)
    x = 6.3*10**9
    t = 6.36*10**9

    print(f"t = {t/(60*60*24*365.25)}")

    t_m = t/l
    print(f"{t_m}")
    print(f"{t_m/(60*60*24*365.25):.3f} Y\n") # Works!!

    # Switching system:
    x,x_m = 0,x
    t,t_m = t_m,t
    v = -v

    t_m_m = t/l

    print(f"{t_m_m}")
    print(f"{t_m_m/(60*60*24*365.25):.3f} Y\n") # Works!!

def Part345():

    c = 3*10**8
    v = 0.99
    l = LorentzFactor(v)
    x_b = 6.3*10**9
    t_b = 6.36*10**9

    t_m = -v*l*x_b + (x_b*l)/v
    print(f"tm_b : {t_m:.3e} Seconds")
    print(f"tm_b : {t_m/(60*60*24*365.25):.3f} Years\n") # Works!!

    t_bm = x_b/v - v*x_b
    print(f"t_bm : {t_bm:.3e} Seconds")
    print(f"t_bm : {t_bm/(60*60*24*365.25):.3f} Years\n") # Works!!

    print(f"tmm_bmm : {(x_b/(v*l))/(60*60*24*365.25):.3f} Years")
    print(f"t_bmm : {((x_b/v + x_b*v))/(60*60*24*365.25):.3f} Years\n")

    g = -0.1/c
    print(f"g : {g:.3e} 1/s")
    print(f"t_b : {t_b:.3e} s")
    t_turn = t_b - v/g
    print(f"t_turn : {t_turn:.3e} s")
    print(f"t_turn : {t_turn/(60*60*24*365.25):.3f} Y\n")

    def t_y_m_acc(t_y):

        v_y = v + g * (t_y - t_b)
        x_y = x_b + v*(t_y - t_b) + 0.5*g*(t_y-t_b)**2

        return t_y - x_y * v_y
    
    def t_y_m_const(t_y):

        v_y = v
        x_y = x_b + v*(t_y - t_b)

        return t_y - x_y * v_y
    
    def t_y_m(t_y):
        
        if t_y > t_b:
            return t_y_m_acc(t_y)
        return t_y_m_const(t_y)

    t_turn_m = t_y_m_acc(t_turn)
    print(f"t_turn_m : {t_turn_m:.3e} s")
    print(f"t_turn_m : {t_turn_m/(60*60*24*365.25):.3f} Y\n")

    dt_turn = t_turn - t_b
    print(f"dt_turn: {(dt_turn)/(60*60*24*365.25):.3f} Y\n")
    t_destiny_return = t_turn + dt_turn
    t_destiny_return_m = t_y_m_acc(t_destiny_return)
    print(f"t_destiny: {(t_destiny_return)/(60*60*24*365.25):.3f} Y")
    print(f"t_destiny: {(t_destiny_return_m)/(60*60*24*365.25):.3f} Y")
    t_final = 202 * 2 + (dt_turn/(60*60*24*365.25))*2
    print(f"t_final: {(t_final):.3f} Y\n")

    # Plotting:

    T_values = np.linspace(0,t_turn,10000)
    T_m_Values = np.zeros_like(T_values)
    for i in range(len(T_values)):
        T_m_Values[i] = t_y_m(T_values[i])

    # plt.plot(T_values/(60*60*24*365.25),T_m_Values/(60*60*24*365.25))
    # plt.grid()
    # plt.xlabel(r"$t_Y$ [år]")
    # plt.ylabel(r"$t_{Y'}$ [år]")
    # plt.show()

    # Measuring Lisa's total time
    g = abs(g)
    print(v)
    print(g*c)
    print()

    T_m = ((v*np.sqrt(1-v**2) + np.arcsin(v))/(2*g)) / (60*60*24*365.25)
    T_total = (2 * T_m) + 57
    print(f"T_m: {(T_m):.3f} Y")
    print(f"T_total: {(T_total):.3f} Y\n")

#Part2()
#Part345()