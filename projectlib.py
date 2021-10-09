#Project Library

#recommended modules (auto-import not possible?)
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt, matplotlib.cm as cm, matplotlib.animation as animation
from IPython.display import HTML
import math
from datetime import datetime as datt, date as date

#Section 1: Networks

def undirect_adjacency(matrix):
    """Turns adjacency matrix of a directed network into adjacency matrix of corresponding undirected network (add ij and ji)

    Args:
        matrix (array): Adjacency matrix of a directed network
    """    
    #Puts sum of ij and ji cell in both ij and ji cell.
    matrix += matrix.transpose()

def unweight_adjacency(matrix):
    """Turns adjacency matrix of a weighted network into simple adjacency matrix of corresponding unweighted network

    Args:
        matrix (array): Adjacency matrix of a weighted network
    """    
    #Replaces non-zero cell entries with ones.
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][j] != 0:
                matrix[i][j] = 1

#Section 2.1: RKI Data Import
#For definition of intern region numbers -> see ..\Media\region 38 with border (cut).png

def region_setup(mode):
    """Assembles necessary AdmUnitIDs for future computations and labels for the chosen Region (12 or 38)
    Args:
        mode (integer): clarifies chosen project Region: 1 for Göttingen only, 12 for Region 12 or 38 for Region 38
    Returns:
        tuple:
            array: contains RKI LK IDs of the respective intern regions (therefore of length mode), index is intern region number
            dictionary: relates intern region number with name of LK
            dictionary: relates intern region number with abbreviation (KFZ-style, no alternative in Germany) of LK
    """    
    lk_hildesheim_id = 3254
    lk_holzminden_id = 3255
    lk_goslar_id = 3153
    lk_höxter_id = 5762
    lk_northeim_id = 3155
    lk_göttingen_id = 3159
    lk_harz_id = 15085
    lk_kassel_id = 6633
    sk_kassel_id = 6611
    lk_werrameißnerkreis_id = 6636
    lk_eichsfeld_id = 16061
    lk_nordhausen_id = 16062
    
    #from here on: Region38
    region_hannover_id = 3241
    lk_schaumburg_id = 3257
    lk_hamelnpyrmont_id = 3252
    lk_lippe_id = 5766
    sk_bielefeld_id = 5711
    lk_gütersloh_id = 5754
    lk_paderborn_id = 5774
    lk_soest_id = 5974
    lk_hochsauerlandkreis_id = 5958
    lk_waldeckfrankenberg_id = 6635
    lk_schwalmederkreis_id = 6634
    lk_hersfeldrotenburg_id = 6632
    lk_wartburgkreis_id = 16063
    sk_eisenach_id = 16056
    lk_unstruthainichkreis_id = 16064
    lk_kyffhäuserkreis_id = 16065
    lk_mansfeldsüdharz_id = 15087
    lk_salzlandkreis_id = 15089
    sk_magdeburg_id = 15003
    lk_börde_id = 15083
    lk_helmstedt_id = 3154
    sk_wolfsburg_id = 3103
    lk_wolfenbüttel_id = 3158
    sk_braunschweig_id = 3103
    sk_salzgitter_id = 3102
    lk_peine_id = 3157

    if(mode == 1):
        region_ids = [lk_göttingen_id]
        region_names = ['Göttingen']
        labels = {}
        short_labels = {}
    elif(mode == 12):
        region_ids = [lk_hildesheim_id,lk_holzminden_id,lk_goslar_id,lk_höxter_id,lk_northeim_id,lk_göttingen_id,lk_harz_id,lk_kassel_id,sk_kassel_id,lk_werrameißnerkreis_id,lk_eichsfeld_id,lk_nordhausen_id]
        region_names = ['Hildesheim','Holzminden','Goslar', 'Höxter','Northeim','Göttingen','Harz','Kassel (LK)','Kassel (SK)', 'Werra-Meißner-Kreis','Eichsfeld','Nordhausen']
        labels = {}
        short_labels = {}
        for intern_region_number in range(mode):
            labels[intern_region_number] = region_names[intern_region_number]
    elif(mode == 38):
        region_ids = [lk_hildesheim_id,lk_holzminden_id,lk_goslar_id,lk_höxter_id,lk_northeim_id,lk_göttingen_id,lk_harz_id,lk_kassel_id,sk_kassel_id,lk_werrameißnerkreis_id,lk_eichsfeld_id,lk_nordhausen_id,region_hannover_id,lk_schaumburg_id,lk_hamelnpyrmont_id,lk_lippe_id,sk_bielefeld_id,lk_gütersloh_id,lk_paderborn_id,lk_soest_id,lk_hochsauerlandkreis_id,lk_waldeckfrankenberg_id,lk_schwalmederkreis_id,lk_hersfeldrotenburg_id,lk_wartburgkreis_id,sk_eisenach_id,lk_unstruthainichkreis_id,lk_kyffhäuserkreis_id,lk_mansfeldsüdharz_id,lk_salzlandkreis_id,sk_magdeburg_id,lk_börde_id,lk_helmstedt_id,sk_wolfsburg_id,lk_wolfenbüttel_id,sk_braunschweig_id,sk_salzgitter_id,lk_peine_id]
        region_names = ['Hildesheim','Holzminden','Goslar', 'Höxter','Northeim','Göttingen','Harz','Kassel (LK)','Kassel (SK)', 'Werra-Meißner-Kreis','Eichsfeld','Nordhausen','Hannover','Schaumburg','Hameln-Pyrmont','Lippe','Bielefeld','Gütersloh','Paderborn','Soest','Hochsauerlandkreis','Waldeck-Frankenberg','Schwalm-Eder-Kreis','Hersfeld-Rotenburg','Wartburgkreis','Eisenach','Unstrut-Hainich-Kreis','Kyffhäuserkreis','Mansfeld-Südharz','Salzlandkreis','Magdeburg','Börde','Helmstedt','Wolfsburg','Wolfenbüttel','Braunschweig','Salzgitter','Peine']
        region_short_names =['HI','HOL', 'GS', 'HX', 'NOM', 'GÖ', 'HZ', 'KS-S', 'KS-L', 'ESW', 'EIC', 'NDH', 'H', 'SHG', 'HM', 'LIP', 'BI', 'GT', 'PB', 'SO', 'HSK', 'KB', 'HR', 'HEF', 'WAK', 'EA', 'UH', 'KYF', 'MSH', 'SLK', 'MD', 'BK', 'HE', 'WOB', 'WF', 'BS', 'SZ', 'PE']
        labels = {}
        short_labels = {}
        for intern_region_number in range(mode):
            labels[intern_region_number] = region_names[intern_region_number]
            short_labels[intern_region_number] = region_short_names[intern_region_number]
    else:
        print("mode err: invalid")
    return region_ids,labels,short_labels

def import_rki_data(region_ids, n):
    """Imports RKI History Data from ../External Data, manipulates it to be (locally) handled in context of the regarded region
    Args:
        region_ids (array): contains RKI AdmUnitIDs of respective LKs (this array may be output of region_setup)
        n (integer): setting for desired n-day-incidence (usually 7)
    Returns:
        tuple:
            array: of the format [intern_region_number][setting][day(since beginning, day zero=2020/03/01)]
                where setting can be:
                    0: new cases
                    1: cumulative no. of cases
                    2: n-day-incidence per 100,000 inhabitants
                    3: new deaths
                    4: cumulative no. of deaths (=current)
                    5: new recoverd
                    6: cumulative no. of recovered (=current)
                    7: current no. of susceptible
                    8: current no. of infected
            array: SIRD compartment distribution of the format [intern_region_number][compartment][day(since beginning, day zero=2020/03/01)]
                warning: invalid/useless for last 14 days of input data!
                where compartment can be:
                    0: S
                    1: I
                    2: R
                    3: D
            array: contains population sizes, index is intern region number
    """    
    #import case data etc. from RKI_History
    rki = pd.read_csv('External Data/RKI_History.csv', sep = ',', header = 'infer')
    rki = rki.sort_values(by = 'Datum')

    #write data in arrays
    lk = np.array(rki['AdmUnitId'])
    lk_comp_num = len(lk)
    day = np.array(rki['Datum'])
    case = np.array(rki['AnzFallErkrankung'])
    case_add = np.array(rki['AnzFallMeldung']) #for additions, may be unused
    case_cum = np.array(rki['KumFall'])

    #import population data and write in array
    pop = pd.read_csv('External Data/RKI_Corona_Landkreise.csv', sep = ',', header = 'infer')
    lk_popcalc = np.array(pop['AdmUnitId'])
    lk_popsize = np.array(pop['EWZ'])
    lk_num = len(lk_popsize)

    #find out number of time steps
    time = len([a for a in lk if a == 0])

    #similar process for RKI_COVID19 (contains information on new deaths and recovered)
    #all variables that were initially dependent on this database carry number 2 in the name

    rki2 = pd.read_csv('External Data/RKI_COVID19_cut_Region38.csv', sep = ',', header = 'infer')
    rki2 = rki2.sort_values(by = 'Meldedatum')
    rki2 = rki2[(rki2['Meldedatum'] >= '2020/03/01')]
    lk2 = np.array(rki2['IdLandkreis'])
    newdead2 = np.array(rki2['AnzahlTodesfall'])
    newrec2 = np.array(rki2['AnzahlGenesen'])
    newcase2 = np.array(rki2['AnzahlFall'])

    #indicator of case type (see documentation on https://www.arcgis.com/home/item.html?id=f10774f1c63e40168479a1feb6c7ca74)
    newdead2check = np.array(rki2['NeuerTodesfall'])
    newrec2check = np.array(rki2['NeuGenesen'])
    newcase2check = np.array(rki2['NeuerFall'])

    #convert dates to days from beginning
    date2 = np.array(rki2['Meldedatum'])
    rki2len = len(lk2)
    date2_day = np.zeros((rki2len))
    for i in range(rki2len):
        if(date2[i] == "2020/03/01"):
            date2_day[i][:10] = 0
        else:
            date2_day[i] = (datt.strptime(date2[i][:10], '%Y/%m/%d') - datt(2020,3,1)).days
        date2_day[i]=int(date2_day[i])

    dimensions = 9 #as shown in description

    #create and fill regional array

    #set up regional case array and regional compartment dirtribution array
    region_num = len(region_ids)
    required_duration = int(max(time,date2_day[-1] +1)) #may be necessary if files are not updated simultaneously
    region_cases = np.zeros((region_num,dimensions, required_duration))
    region_popsize = np.zeros((region_num))
    region_compartment_distribution = np.zeros((region_num,4, required_duration))

    testcount=0
    for intern_region_number in range(region_num): #for rki and rki2, two different processes
        current_time = 0
        
        for i in range(lk_comp_num):
            if lk[i] == region_ids[intern_region_number] and current_time < time:
                region_cases[intern_region_number][0][current_time] = case[i] + case_add[i]
                #region_cases[intern_region_number][1][current_time] = case_cum[i]
                current_time += 1  
        for k in range(rki2len):
            if lk2[k] == region_ids[intern_region_number]: #follow documentation on https://www.arcgis.com/home/item.html?id=f10774f1c63e40168479a1feb6c7ca74
                testcount+=1
                #if(newcase2check[k]==1 or newcase2check[k]==0):
                    #region_cases[intern_region_number][0][date2[k]]+=newcase2[k]
                if(newdead2check[k] == 1 or newdead2check[k] == 0):
                    region_cases[intern_region_number][3][int(date2_day[k])] += newdead2[k]
                if((newrec2check[k] == 1 or newrec2check[k] == 0) and (date2_day[k] + 14 < date2_day[-1] +1 )):
                    region_cases[intern_region_number][5][int(date2_day[k])+14] += newrec2[k]
        for j in range(lk_num):
            if lk_popcalc[j] == region_ids[intern_region_number]:
                region_popsize[intern_region_number] = lk_popsize[j]
                region_cases[intern_region_number][2] = n_day_incidence(region_cases[intern_region_number][0],region_popsize[intern_region_number],n)
        region_cases[intern_region_number][1] = cumulate_data(region_cases[intern_region_number][0])
        region_cases[intern_region_number][4] = cumulate_data(region_cases[intern_region_number][3])
        region_cases[intern_region_number][6] = cumulate_data(region_cases[intern_region_number][5])
        region_cases[intern_region_number][7] = region_popsize[intern_region_number] * np.ones((required_duration)) - region_cases[intern_region_number][1]
        region_cases[intern_region_number][8] = region_cases[intern_region_number][1] - region_cases[intern_region_number][4] - region_cases[intern_region_number][6]
        region_compartment_distribution[intern_region_number][0] = region_cases[intern_region_number][7] / region_popsize[intern_region_number]
        region_compartment_distribution[intern_region_number][1] = region_cases[intern_region_number][8] / region_popsize[intern_region_number]
        region_compartment_distribution[intern_region_number][2] = region_cases[intern_region_number][6] / region_popsize[intern_region_number]
        region_compartment_distribution[intern_region_number][3] = region_cases[intern_region_number][4] / region_popsize[intern_region_number]
    return region_cases, region_compartment_distribution, region_popsize


def update_rki_data_arrays(mode,n):
    """Saves results from import_rki_data in arrays directly in the Internal Data directory, warning: invalid/useless for last 14 days of input data"""

    rki_region_cases, rki_region_compartment_distribution,pop = import_rki_data(region_setup(mode)[0],n)
    if(mode == 12):
        np.save('Internal Data/rki_region_cases12.npy',rki_region_cases)
        np.save('Internal Data/rki_region_compartment_distribution12.npy',rki_region_compartment_distribution)
    elif(mode == 38):
        np.save('Internal Data/rki_region_cases38.npy',rki_region_cases)
        np.save('Internal Data/rki_region_compartment_distribution38.npy',rki_region_compartment_distribution)
    else:
        print('mode error: invalid argument')

def initial_compartment_distribution(mode, date):
    """provides the requested initial compartment distribution, warning: invalid/useless for last 14 days of input data

    Args:
        mode (integer): regarded region ()
        date (string): date that is to be considered day zero in the format: "YYYY/MM/DD"

    Returns:
        array: initial compartment distribution:
            indices: 0:S, 1:I, 2:R, 3:D (as usual)
    """    
    day=int((datt.strptime(date[:10], '%Y/%m/%d') - datt(2020,3,1)).days)
    if(mode == 12):
        return np.load('Internal Data/rki_region_compartment_distribution12.npy')[:,:,day]
    elif(mode == 38):
        return np.load('Internal Data/rki_region_compartment_distribution38.npy')[:,:,day]
    else:
        print("mode error: invalid argument")

def save_relevant_timeline():
    """returns SIRD distribution of Region 38 of each day from Jul 24, 2020 to Nov 1, 2020

    Returns:
        see above
    """ 

    day_zero=int((datt(2020,7,24) - datt(2020,3,1)).days)
    output=np.zeros((38,4,100))
    for d in range(100):
        for n in range(38):
            for m in range(4):
                output[n][m][d] = np.load('Internal Data/rki_region_compartment_distribution38.npy')[n,m,day_zero+d]
    np.save('Internal Data/timeline.npy',output)

def save_initial_compartment_distribution(mode, date):
    #name = date+".txt"
    file = open("initial_data.txt","w")
    np.savetxt(file, initial_compartment_distribution(mode, date))
    file.close()

def save_popdata():
    file = open("popdata38.txt", "w")
    np.savetxt(file, import_rki_data(region_setup(38)[0], 7)[2])
    file.close()

#Section 2.2: Data Manipulation

def cumulate_data(case_array):
    """ "integrates" new cases etc. over time and replaces new cases with accumulated cases

    Args:
        case_array (array): development of new cases over time

    Returns:
        array: development of cumulated cases over time
    """    
    cumulated = case_array.copy()
    sum = 0
    for i in range(len(cumulated)):
        sum += cumulated[i]
        cumulated[i] = sum
    return cumulated

def n_day_moving_average(case_array,n):
    """Computes moving average of any time-dependent quantity (mostly cases) over n days in one cell
            Note (1): undefined for first (n-1) days, set to zero
            Note (2): this function will most likely only be used by n_day_incidence

    Args:
        case_array (array): temporal development of a quantity
        n (integer): days over which the average is computed (usually 7)

    Returns:
        array: moving average of the quantity in case_array over n days in one cell, computed for each day (=each entry)
    """    
    output=np.zeros((len(case_array)))
    for i in range(len(case_array)):
        if(i>=n-1):
            isum=0
            for j in range(n):
                isum+=case_array[i-j]
            output[i]=isum
        else:
            output[i]=0 #or undefined, not sure tbh
    return output/n

def n_day_incidence(case_array,pop,n):
    """calculates n-day-incidence per 100,000 inhabitants in one cell

    Args:
        case_array (array): temporal development of a quantity
        pop (integer): population size of the regarded cell
        n (integer): setting for desired n-day-incidence (usually 7)

    Returns:
        array: contains n-day-incidence per 100,000 inhabitants, computed for each day (=each entry)
    """    
    #calculates n-day-incidence per 100,000 inhabitants
    return n*n_day_moving_average(case_array,n)/pop*100000


#Section 2.3: Commuter related functions

def remaining_population(comFrom, N, i):
    """Function to calculate the remaining population in a cell during commuting stage = N_i^{rest} in document

    Args:
        comFrom (callable(j)): function that returns commuters from cell.
        N (array): array with population of every cell
        i (integer): the cell of which the remaining population is wanted.

    Returns:
        float: the remaining population in cell i.
    """
    return N[i] - np.sum(comFrom(i))



def contacts_infected(comTo, comFrom, N, i, infected, dimension):
    """Function to calculate the relative number of infected of the population from a cell during commuting stage = I_i^{pen} in document

    Args:
        comTo (callable(j)): function that returns array of commuters to cell.
        comFrom (callable(j)): function that returns array of commuters from a cell.
        N (array): array with population of every cell
        i (integer): the cell of which the number of infectious contacts is wanted.
        infected (array): array with number of infected from every cell
        dimension (integer): number of cells

    Returns:
        float: infectious contacts of cell i inhabitants
    """

    # remaining population in cell i
    Nrest = remaining_population(comFrom, N, i)

    # array for commuters to cell i
    cto = comTo(i)
    
    #variables for sum of commuters & sum of infectious commuters
    sto = 0
    sto_inf = 0

    for j in range(dimension):
        # adding commuters from all cells to i
        sto += cto[j]

        # adding infectious commuters from all cells to i
        sto_inf += cto[j] * infected[j]
    
    # absolute number of infectious contacts
    contacts = Nrest * infected[i] + sto_inf
    
    # normalizing
    Ipen = contacts / (Nrest + sto)

    return Ipen


def effective_infected(comTo, comFrom, N, i, infected, dimension):
    """Function to calculate the effective number of infected in a cell = I_i^{eff} in document

    Args:
        comTo (callable(j)): function that returns array of commuters to cell.
        comFrom (callable(j)): function that returns array of commuters from a cell.
        N (array): array with population of every cell
        i (integer): the cell of which the effective population is wanted.
        infected (array): array with number of infected from every cell
        dimension (integer): number of cells

    Returns:
        float: effective infected in cell i
    """

     # remaining population in cell i
    Nrest = remaining_population(comFrom, N, i)

    # array for commuters from cell i
    cfrom = comFrom(i)

    # array to be filled with the infectious contacts
    Ipen = np.zeros(dimension)

    # variable for sum of infected commuters
    spen_inf = 0

    for j in range(dimension):
        # infectious contacts of cell j inhabitants
        Ipen[j] = contacts_infected(comTo, comFrom, N, j, infected, dimension)

        # adding commuters infected in other cells
        spen_inf += cfrom[j] * Ipen[j]

    #  effective infected in cell i, normalized
    Ieff = (Nrest * Ipen[i] + spen_inf) / N[i]

    return Ieff


def periodic_heaviside(t, t0):
    """Function to make a periodic heaviside. The period is 1.

    Args:
        t (float): current time
        t0 (float): time that divides the heaviside, i.e. t >= t0 -> 1 and t < t0 -> 0 

    Raises:
        ValueError: To ensure that the give time is positive
        ValueError: To ensure that t0 is between 0 and 1

    Returns:
        float: either 1 or 0 depending on t relativ to t0
    """
    # returning error in case t isn't positive
    if t < 0 or t0 < 0:
        raise ValueError("Time must be positive")
    
    # the t0 can only be between 0 and 1. Thus, an error will appear if this is not fullfilled
    elif t0 > 1:
        raise ValueError("t0 must be between 0 and 1")
    
    # in case of t being between 0 and 1 there is no need of further manipulation
    elif t < 1 and t >= 0:
        return np.heaviside(t - t0, 1)
    
    # for t > 1 it needs to be modified to fit the periodic character of the function
    else:
        # ignoring everything that isn't a decimal in t, as the periode is 1
        string = str(t).split(".")
        # the new time is thusly 
        t = float("0." + string[1])
        return np.heaviside(t - t0, 1)


#Section 2.4: Class

class Data:
    def __init__(self, alpha, beta, dimension, name, commuters_day, N):
        "Dimension gives the dimension of the system, i.e. number of cells"
        # ------- TO DO ---------
        # Depending on the system this has to be adjusted
        self.dimension = dimension
        
        "Alpha describes the rate of infection"
        self.alpha = alpha # Dimensions day^(-1)

        "Beta describes the recovery rate"
        self.beta = beta # Dimensions day^(-1)

        "p is the probability to die from the disease"
        self.p = 2.64e-2
        
        """
        Commuters is an array describing the commuters from and to the different cells. 
        The structure is: The entries in column i are how many commute from i to the cell of the row.
                          Thus, the entries of the row i describe how many commute from different cells to cell i.
        """

        self.commuters = np.array([np.ones(self.dimension) for i in range(self.dimension)])
        "load the txt with given name"
        file = np.loadtxt(name, delimiter="\t")
        "fill commuters with the entries from the file. Note the break statement so it doesn't go to far"
        i = 0
        for row in file:
            if i >= self.dimension:
                break
            self.commuters[i] = np.asarray(row[0:self.dimension])
            i += 1
        
        "Part of day that commuters are in other cells. This is the same value for all cells"
        self.commuters_day = commuters_day

        "Array of population in every cell"
        # ------- TO DO ---------
        # figure out how to best do this
        self.N = N

        
        
    def commutersFrom(self, cfrom):
        """Function to extract the commuters from a cell

        Args:
            cfrom (integer): the cell 

         Returns:
            array: array containing commuters from cfrom
        """
        # ------- TO DO ---------
        # optimize the algorithm
        "Sorting out the necessary entries into array a."
        a = np.zeros(self.dimension)
        for i in range(self.dimension):
            a[i] = self.commuters[i][cfrom]
            
        return a 
    
    def commutersTo(self, cto):
        """Function to return an array with number of commuters coming to cto from other cells.

        Args:
            cto (integer): the cell to which the commuters travel

        Returns:
            array: array containing the commuters to cto.
        """
        return self.commuters[cto][:]



#Section 2.5: ODE functions

def simple_system(timestep, functions, data):
    """Function to calculate the time derivatives of the SIRD modell without commuters. The system may be one specific district (LK), several districts combined or the whole area, this is up to the user (input).

        Args:
            timestep (float): time of timestep. t=1 is 1 day
            functions (array or list): array with the scalar functions (compartments) of the system. Structure [S, I, R, D]
            data (class): the class with the relevant data.

        Returns:
            list: list containing the derivatives of functions at time t.
    """
    S,I,R,D = functions
    dfuncdt = [0,0,0,0]

    # dS/dt
    dfuncdt[0] =  - data.alpha * S * I

    # dI/dt \propto -dSdt
    dfuncdt[1] = - dfuncdt[0] - data.beta * I

    # dR/dt
    dfuncdt[2] = (1 - data.p) * data.beta * I

    # dD/dt
    dfuncdt[3] = data.p * data.beta * I

    return dfuncdt   


def function_of_system(timestep, functions,  data, method, t0):
    """Function to calculate the time derivatives of the dynamic SIRD modell with commuters. This function uses multiple approachs for different variant of the system based on prior assumptions. 

        Args:
            functions (array): array with the functions that define the system of equations. Structure [S_1, ... , S_n, I_1, .... , I_n, R_1, ... , R_n, D_1, ... , D_n]
            timestep (float): time of timestep. t=1 is 1 day
            data (class): the class with the relevant data
            method (string): dictates what method is to be used. The choices are "constant", "heaviside"
            t0 (array): array with len = 3 filled with times at which the heaviside functions in the "heaviside" method switch. 
                Note: all values of t0 have to be between 0 and 1.

        Returns:
            array: array containing the derivatives of functions at time t.
    """
    # array with infected
    Infected = functions[data.dimension:2*data.dimension]

    # initializing the return array with the time derivatives of functions
    dfuncdt = np.zeros(data.dimension*4)

    "~~~~~~~~~~~~~~~~~ See the PDF / LaTeX for further explanation of the equations ~~~~~~~~~~~~~~~~~~~~~~~~"

    # ------------------------------------
    # -     Without commuting effects    -
    # ------------------------------------
    
    if method == 'simple':

        "Please use the function simple_system for the method 'simple'"
    
        # dS/dt
        dfuncdt[0:data.dimension] =  - data.alpha * functions[:data.dimension] * Infected

        # dI/dt \propto -dSdt
        dfuncdt[data.dimension:2*data.dimension] = - dfuncdt[0:data.dimension] - data.beta * Infected

        # dR/dt
        dfuncdt[2*data.dimension:3*data.dimension] = (1 - data.p) * data.beta * Infected

        # dD/dt
        dfuncdt[3*data.dimension:] = data.p * data.beta * Infected
            
        
    
    
    # making array with the effective infected Ieff
    Ieff = np.array([effective_infected(data.commutersTo, data.commutersFrom, data.N, i, Infected, data.dimension) for i in range(data.dimension)])
        
    # ------------------------------------
    # -      Constant Coefficients       -
    # ------------------------------------

    if method == "constant":

        t_out = data.commuters_day #time in another vertex
            
        # the for loop to fill dfuncdt

        for i in range(data.dimension):
            # see LaTeX for equations
            # dS/dt
            dfuncdt[i] =  - (1 - t_out) * data.alpha *  functions[i] * functions[i + data.dimension] - t_out * data.alpha * functions[i] * Ieff[i]

            # dI/dt \propto -dSdt
            dfuncdt[i + data.dimension] = - dfuncdt[i] - data.beta * functions[i + data.dimension]

            # dR/dt
            dfuncdt[i + 2 * data.dimension] = (1 - data.p) * data.beta * functions[i + data.dimension]

            # dD/dt
            dfuncdt[i + 3 * data.dimension] = data.p * data.beta * functions[i + data.dimension]


    # ------------------------------------
    # -      Heaviside Coefficients      -
    # ------------------------------------

    elif method == "heaviside":


        # the for loop to fill dfuncdt
        for i in range(data.dimension):
            # see LaTeX for equations
            # dS/dt
            dfuncdt[i] = - data.alpha * (periodic_heaviside(timestep, t0[0]) + periodic_heaviside(timestep, t0[2]) - 1) * functions[i] * functions[i + data.dimension] - (periodic_heaviside(timestep, t0[1]) - periodic_heaviside(timestep, t0[2])) * data.alpha * functions[i] * Ieff[i]

            # dI/dt \propto -dSdt
            dfuncdt[i + data.dimension] = - dfuncdt[i] - data.beta * functions[i + data.dimension]

            # dR/dt
            dfuncdt[i + 2 * data.dimension] = (1 - data.p) * data.beta * functions[i + data.dimension]

            # dD/dt
            dfuncdt[i + 3 * data.dimension] = data.p * data.beta * functions[i + data.dimension]

    
    # returning the derivatives at time t
    return dfuncdt


def system_function(functions, timestep, data, method, t0):
    # Function for using scipy.integrate.odeint, as the arguments are switched

    return function_of_system(timestep, functions, data, method, t0)


def rk4(func,t,y,dt):
    # Runge-Kutta Method Single Time Step
    """
    Peforms a single time step of the Runge-Kutta 4th Order Method.
    The below function finds the ki value for func and return the value to move Yn+1
    func is the rhs of Y'
    
    Recall Rk4 Equations :
    k1 = h*func(tn,yn)
    k2 = h*func(tn+h/2,yn+k1/2)
    k3 = h*func(tn+h/2,yn+k2/2)
    k4 = h*func(tn,yn+k3)
    Where func is the function of the derivative(s) dy/dt
    Yn+1 = Yn + 1/6*(k1+k2+k3+k4)
    """

    k1 = dt*func(t,y)
    k2 = dt*func(t+dt/2,y+k1/2)
    k3 = dt*func(t+dt/2,y+k2/2)
    k4 = dt*func(t+dt,y+k3)

    return (1/6)*(k1+k2+k3+k4)


#Section 2.6: Jacobian

def jacobi(data, comp, model, t0):


    # TODO need to finish this function. Remaining are almost all non-trivial blocks

    # initializing the matrix that will be returned

    #jac = np.array([np.zeros(data.dimension * 4) for i in range(data.dimension * 4)])

    # initializing zero blocks of jacobian

    """
    zeroBlock = np.array([np.zeros(data.dimension) for i in range(data.dimension)])
    SS = np.array([np.zeros(data.dimension) for i in range(data.dimension)])
    SI = SS
    IS = SS
    II = SS
    RI = SS
    DI = SS
    SR = SS
    SD = SS
    IR = SS
    ID = SS
    RS = SS
    RR = SS
    RD = SS
    DS = SS
    DR = SS
    DD = SS




    if model == "constant":
        for i in range(data.dimension):

            #variable for sum in SS
            m = 0
            for k in range(data.dimension):
                m += data.commutersFrom(i)[k] * effective_infected(data.commutersTo, data.commutersFrom, data.N, k, comp[data.dimension:2 * data.dimension], data.dimension)
            
            SS[i, i] = - ((1 - data.commuters_day)*data.alpha*comp[data.dimension + i] + \
                    data.commuters_day*data.alpha  / data.N[i] * (effective_population(data.commutersFrom, data.N, i) * \
                    effective_infected(data.commutersTo, data.commutersFrom, data.N, i, comp[data.dimension:2 * data.dimension], data.dimension) + m))
            
            IS[i][i] = - SS[i][i]
            for k in range(data.dimension):
                
        
        
        
        
        IS = SS

    elif model == "heaviside":

    
    
    jac = np.bmat([[SS, SI, SR, SD], [IS, II, IR, ID], [RS, RI, RR, RD], [DS, DI, DR, DD]])

    """
    return jac


#Section 2.7: Data Visualization








#Section 2.8: Comparing Data



def squared_variance(numerical_data, external_data):
    """Simple function to calculate the variance, or deviation, of two datasets

    Args:
        numerical_data (array): numerical dataset
        external_data (array): dataset to compare with the numerical dataset 

    Raises:
        ValueError: arrays have to be the same length

    Returns:
        float: the deviation of the datasets
    """
    
    sum = 0

    if len(numerical_data) != len(external_data):
        raise ValueError("Arrays have to be of same length")

    for i in range(len(numerical_data)):
        sum += (numerical_data[i] - external_data[i])**2
    
    return sum


