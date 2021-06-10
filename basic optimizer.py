import numpy as np


def isbetter(new_eval, current_eval):
    if new_eval > current_eval:
        return 1
    else:
        return 0


def Bew(params, component, which_one):
    if which_one == -1:
        return f(params)
    else:
        z=params
        z[which_one] = component
        return f(z)


def I_neu(d, params):
    state=summon_results(params)
    S = state[0:n - 1]
    I = state[n:2 * n - 1]
    R = state[2 * n:3 * n - 1]
    D = R * np.array(p)
    sum=0
    for i in range(d):
        sum=sum+I[i]+R[i]+D[i]-I[i-1]


def f(params):
    number_of_days = 100
    sum=0
    for d in range(number_of_days):
        sum = sum + I_neu(d, params)**2
    return sum


def imporve(params, component,which_one, sign, eval):

    i = 0
    d = 0
    last_sucess = 0
    new_eval = eval
    value = component
    o_err=-1
    while i - last_sucess < 15 and o_err < 14:
        if i - last_sucess >8:
            o_err = o_err + 1
            last_sucess = 0
            i = 0
            component = value
        d = i * 10**(-o_err) * sign
        component = value + d
        eval=new_eval
        new_eval = Bew(params, component, which_one)
        i = i + 1
        if new_eval < eval:
            last_sucess=last_sucess + 1
            value=component
        print(component, new_eval,eval)
    return value, new_eval


def improve_component(params, which_one):
    if which_one ==-1:
        component = params
    else:
        component = params[which_one]
    component_pos, eval1 = imporve(params, component,which_one, 1, Bew(params, component, which_one))
    component_neg, eval2 = imporve(params, component_pos,which_one, -1, eval1)
    if eval1<eval2:
        return component_pos
    else:
        return component_neg







def optimize(parameters):
    improved_paramters=np.array(parameters)
    for i in range(len(parameters)):
        improved_paramters[i] = improve_component(parameters, i)
    return improve_component(parameters, -1)


