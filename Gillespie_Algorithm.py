import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#A Gillespie Approach to Acid-Base Bromophenol Blue Indicator Equilibrium
#The optimal concentration of bromophenol blue is 2.4*10**-5 M
#The pH of the solutions are 2.0, 4.0 and 13.0 for acid, alkali, and buffer respectively
#The Acid-Base bromophenol blue indicator equilibrium is given by the equation HIn- ⇌ H+ + In2-
#BPB indicator is anionic by nature


Acid_spectra = np.loadtxt("acidspec1.txt", delimiter=",", skiprows=6)
Alkali_spectra = np.loadtxt("alkalispec1.txt", delimiter=",", skiprows=6)
Buffer_spectra = np.loadtxt("bufferspec1.txt", delimiter=",", skiprows=6)

Wavelength_acid = Acid_spectra[:,0]
Absorbance_acid = Acid_spectra[:,1]

Wavelength_alkali = Alkali_spectra[:,0]
Absorbance_alkali = Alkali_spectra[:,1]

Wavelength_buffer = Buffer_spectra[:,0]
Absorbance_buffer = Buffer_spectra[:,1]

plt.plot(Wavelength_acid, Absorbance_acid, linestyle = '-', color = 'Red')
plt.xlim(425,750)
plt.ylim(-0.5,1.0)
plt.xlabel("Wavelength / nm")
plt.ylabel("Absorbance")
plt.title("A graph of absorbance vs wavelength for an acidic solution")
plt.show()

plt.plot(Wavelength_alkali, Absorbance_alkali, linestyle = '-', color = 'Blue')
plt.xlim(425,750)
plt.ylim(-0.5,2.0)
plt.xlabel("Wavelength / nm")
plt.ylabel("Absorbance")
plt.title("A graph of absorbance vs wavelength for a basic solution")
plt.show()

plt.plot(Wavelength_buffer, Absorbance_buffer, linestyle = '-', color = 'Green')
plt.xlim(425,750)
plt.ylim(-0.5,1.0)
plt.xlabel("Wavelength / nm")
plt.ylabel("Absorbance")
plt.title("A graph of absorbance vs wavelength for a buffer solution")
plt.show()

plt.plot(Wavelength_acid, Absorbance_acid, linestyle = '-', color = 'Red', label = "Acid")
plt.plot(Wavelength_alkali, Absorbance_alkali, linestyle = '-', color = 'Blue', label = "Alkali")
plt.plot(Wavelength_buffer, Absorbance_buffer, linestyle = '-', color = 'Green', label = "Buffer")
plt.xlim(425,750)
plt.ylim(-0.5,2.0)
plt.xlabel("Wavelength / nm")
plt.ylabel("Absorbance")
plt.title("A graph of absorbance vs wavelength for all three solutions")
plt.legend()
plt.show()

def find_absorbance(data):
    absorbance = data[:, 1]
    peak_abs = float("-inf")
    for i in range(len(absorbance)):
        if absorbance[i] != absorbance[i]:
            continue
        if data[i, 0] < 525 or data[i, 0] > 625:
            continue
        if absorbance[i] > peak_abs:
            peak_abs = absorbance[i]
    return peak_abs

print("The value of peak absorbance for each environment:")
print("Acid - "+str(find_absorbance(Acid_spectra)))
print("Alkali - "+str(find_absorbance(Alkali_spectra)))
print("Buffer - "+str(find_absorbance(Buffer_spectra)))

def find_wavelength(data):
    absorbance = data[:, 1]
    wavelength = data[:, 0]
    peak_abs = float("-inf")
    lambda_max = 0
    for i in range(len(absorbance)):
        if absorbance[i] != absorbance[i]:
            continue
        if data[i, 0] < 525 or data[i, 0] > 625:
            continue
        if absorbance[i] > peak_abs:
            peak_abs = absorbance[i]
            lambda_max = wavelength[i]
    return lambda_max

print("The value of λmax for each environment:")
print("Acid - "+str(find_wavelength(Acid_spectra)))
print("Alkali - "+str(find_wavelength(Alkali_spectra)))
print("Buffer - "+str(find_wavelength(Buffer_spectra)))

def abs_at(data, target):
    wavelength = data[:, 0]
    absorbance = data[:, 1]
    for i in range(len(wavelength)):
        if abs(target - wavelength[i]) < 0.1:
            return absorbance[i]

iso_acid = abs_at(Acid_spectra, 589)
iso_alkali = abs_at(Alkali_spectra, 589)
iso_buffer = abs_at(Buffer_spectra, 589)
print("The absorbance value for each environment at target (isoelectric) point:")
print("Acid - "+str(iso_acid))
print("Alkali - "+str(iso_alkali))
print("Buffer - "+str(iso_buffer))

def fraction(abs1, abs2, abs3):
    f = (abs3 - abs1) / (abs2 - abs1)
    return f

bpb_f_dep = fraction(iso_acid, iso_alkali, iso_buffer)
print("Fracton of molecules in the deprotonated form to 3 d.p")
print(round(bpb_f_dep, 3))

def henderson_hasselbach(pH, f):
    pKa = pH - np.log10(f/(1-f))
    return pKa

bpb_pKa = henderson_hasselbach(4.0, bpb_f_dep)
print("pKa of bromophenol blue indicator at pH where HIn and In- simultaneously exist to 3 d.p")
print(round(bpb_pKa, 3))

def Ka_calc(pKa):
    Ka = 10**-pKa
    return Ka

Ka_of_eqm = Ka_calc(bpb_pKa)
print("Ka (dissociation constant) of acid-base bromophenol blue indicator equilibrium")
print(Ka_of_eqm)

bpb_f_p = 1 / (1 + (Ka_of_eqm/10**-4))
print("Fraction of molecules in the protonated form to 3 d.p")
print(round(bpb_f_p, 3))

def Gillespie(N, t_max, forward_k, reverse_k):
    t = 0
    N_prot = N
    N_deprot = 0
    time = []
    deprot_mols = []
    while t < t_max:
        a1 = forward_k * N_prot
        a2 = reverse_k * N_deprot
        a0 = a1 + a2
        if a0 == 0:
            break
        τ = np.random.exponential(1/a0)
        choose = np.random.uniform(0, 1)
        if choose < a1/a0:
            N_prot -= 1
            N_deprot += 1
        else:
            N_prot += 1
            N_deprot -= 1
        t += τ
        time.append(t)
        deprot_mols.append(N_deprot)
    return time, deprot_mols

ten_molecules = Gillespie(10, 1000, 1, 1.486)
hundred_molecules = Gillespie(100, 1000, 1, 1.486)
thousand_molecules = Gillespie(1000, 1000, 1, 1.486)
ten_thousand_molecules = Gillespie(10000, 1000, 1, 1.486)

plt.scatter(ten_molecules[0], ten_molecules[1], s = 2, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Deprotonated molecule count with time where N=10")
plt.show()

plt.scatter(hundred_molecules[0], hundred_molecules[1], s = 2, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Deprotonated molecule count with time where N=100")
plt.show()

plt.scatter(thousand_molecules[0], thousand_molecules[1], s = 1, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Deprotonated molecule count with time where N=1000")
plt.show()

plt.scatter(ten_thousand_molecules[0], ten_thousand_molecules[1], s = 1, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Deprotonated molecule count with time where N=10000")
plt.show()

def stationary_count_split(data, N, window, flat, repeat):
    previous_mean = 0
    current_mean = 0
    persist_count = 0
    for i in range(2*window, len(data)):
        current_mean = sum(data[i-window:i]) / window
        previous_mean = sum(data[i-2*window:i-window]) / window
        derivative = current_mean - previous_mean
        if abs(derivative/N) < flat:
            persist_count += 1
        else:
            persist_count = 0
        if persist_count >= repeat:
            return i
    return len(data) // 2

#I want to make sure this function works with a small test

test = [0,1,1,2,2,3,3,4,4,4,3,4,5,4,4,3,4,5,4,4,4,3,4,5,4,3,4,4,5,4]

print("Test run of my stationary_count_split function:")
print(stationary_count_split(test, 10, 2, 0.1, 3))

#This should return a value of 12 (it did for me)

def stationary_count(N, t_max, forward_k, reverse_k, window, flat, repeat):
    time, deprot_mols = Gillespie(N, t_max, forward_k, reverse_k)
    split = stationary_count_split(deprot_mols, N, window, flat, repeat)
    return time[split:], deprot_mols[split:]

still_ten_molecules = stationary_count(10, 1000, 1, 1.486, 2, 0.5, 5)
plt.scatter(still_ten_molecules[0], still_ten_molecules[1], s = 2, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Stationary deprotonated molecule count with time where N=10")
plt.show()

still_hundred_molecules = stationary_count(100, 1000, 1, 1.486, 10, 0.01, 10)
plt.scatter(still_hundred_molecules[0], still_hundred_molecules[1], s = 2, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Stationary deprotonated molecule count with time where N=100")
plt.show()

still_thousand_molecules = stationary_count(1000, 1000, 1, 1.486, 50, 0.0001, 25)
plt.scatter(still_thousand_molecules[0], still_thousand_molecules[1], s = 2, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Stationary deprotonated molecule count with time where N=1000")
plt.show()

still_ten_thousand_molecules = stationary_count(10000, 1000, 1, 1.486, 500, 0.000001, 50)
plt.scatter(still_ten_thousand_molecules[0], still_ten_thousand_molecules[1], s = 2, alpha = 0.5)
plt.xlabel("time")
plt.ylabel("Deprotonated molecule count")
plt.title("Stationary deprotonated molecule count with time where N=10000")
plt.show()

def binomial_plot(N):
    return np.arange(-0.5, N + 1.5, 1)

#Due to the distribution being binomial in nature i.e discrete, some values will not be integers

plt.hist(still_ten_molecules[1], bins = binomial_plot(10), density = True)
plt.xlabel("Number of deprotonated molecules")
plt.ylabel("Probability density")
plt.title("A distribution of deprotonated molecules when N=10")
plt.show()

plt.hist(still_hundred_molecules[1], bins = binomial_plot(100), density = True)
plt.xlabel("Number of deprotonated molecules")
plt.ylabel("Probability density")
plt.title("A distribution of deprotonated molecules when N=100")
plt.show()

plt.hist(still_thousand_molecules[1], bins = 40, density = True)
plt.xlabel("Number of deprotonated molecules")
plt.ylabel("Probability density")
plt.title("A distribution of deprotonated molecules when N=1000")
plt.show()

plt.hist(still_ten_thousand_molecules[1], bins = 40, density = True)
plt.xlabel("Number of deprotonated molecules")
plt.ylabel("Probability density")
plt.title("A distribution of deprotonated molecules when N=10000")
plt.show()

ten_mols_fractions = [mol/10 for mol in still_ten_molecules[1]]
hundred_mols_fractions = [mol/100 for mol in still_hundred_molecules[1]]
thousand_mols_fractions = [mol/1000 for mol in still_thousand_molecules[1]]
ten_thousand_mols_fractions = [mol/10000 for mol in still_ten_thousand_molecules[1]]

shared_edges = np.arange(0.2, 0.6, 0.005)

plt.hist(ten_mols_fractions, bins = shared_edges, alpha = 0.5, density = True, label = "N=10")
plt.hist(hundred_mols_fractions, bins = shared_edges, alpha = 0.5, density = True, label = "N=100")
plt.hist(thousand_mols_fractions, bins = shared_edges, alpha = 0.5,  density = True, label = "N=1000")
plt.hist(ten_thousand_mols_fractions, bins = shared_edges, alpha = 0.5, density = True, label = "N=10000")
plt.legend()
plt.xlabel("Deprotonated molecule in proportion to N")
plt.ylabel("Probability density")
plt.title("A distribution of deprotonated molecules in proportion to N at different N")
plt.show()

#I now want to plot a 3D model of these distributions for a clearer view of how they compare

height_ten, edges_ten = np.histogram(ten_mols_fractions, bins = shared_edges, density = True)
height_hundred, edges_hundred = np.histogram(hundred_mols_fractions, bins = shared_edges, density = True)
height_thousand, edges_thousand = np.histogram(thousand_mols_fractions, bins = shared_edges, density = True)
height_ten_thousand, edges_ten_thousand = np.histogram(ten_thousand_mols_fractions, bins = shared_edges, density = True)

centres_of_bars = (shared_edges[:-1] + shared_edges[1:]) / 2
heights_of_bars = [height_ten, height_hundred, height_thousand, height_ten_thousand]
depth_of_bar = [np.log10(10), np.log10(100), np.log10(1000), np.log10(10000)]

ax = plt.axes(projection='3d')
mpl.use('TkAgg')
for y_values, z_values in zip(depth_of_bar, heights_of_bars):
    ax.plot(centres_of_bars, [y_values]*len(centres_of_bars), z_values)
ax.set_xlabel("Deprotonated molecule in proportion to N")
ax.set_ylabel("Log10(N)")
ax.set_zlabel("Probability Density")
plt.show()

#This 3D plot CAN be interactive in PyCharm, just remember to turn off 'Show in window' in python plots


