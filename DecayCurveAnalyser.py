####################################################
#               Decay Curve analyser
# Python scripts used for standard analysis in Luminescence Experiments

# Features:
#   - Decay Curves from a Time-Correlated Single Photon Counter
#       - Plot and Analyser
#       - Calculations of the Lifetime with some methods *Under Construction*
#
# Time-correlated photon counter home-assembled

# By Allison Pessoa, Nano-Optics Laboratory.
# UFPE, Brazil.
####################################################

import numpy as np

import scipy.integrate as integrate
import scipy.signal as signal

import uncertainties
from uncertainties.umath import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class DecayCurve():
# Contains all information about each individual decay curve
    def __init__(self, file_name):
        self.filename = file_name
        (self.time, self.intensity) = self._load_file()

    def _load_file(self):
        #Loads and processes the decay curves .txt files "time,intensity\n"
        decay_curve_file = open(self.filename)
        time = []
        intensity = []

        for line in decay_curve_file.readlines():
            time.append(float(line.split(',')[0]))
            intensity.append(float(line.split(',')[1]))

        decay_curve_file.close()
        return (time,intensity)

    def get_decay_curve(self, normalized=True):
        #Returns a tuple with the decay curve data: (time, intensity)
        #If normalized is True, intensity will be normalized to values between 0 and 1
        intensity = []

        if normalized == True:
            for j in range(len(self.intensity)):
                intensity.append((self.intensity[j]-min(self.intensity))/(max(self.intensity)-min(self.intensity)))
        else:
            intensity = self.intensity

        return (self.time, intensity)

    def plot_decay_curve(self, normalized=False, v_lines = None, h_lines = None):
        #Plots the decay curves. Returns the matplotlib object for any alterations, if necessary.
        #If normalized is True, intensity will be normalized to values between 0 and 1

        #Setups for decay curves plots
        fig, ax = plt.subplots()
        ax.set_xlabel('Time [$\mu$s]', fontsize='x-large')
        ax.set_ylabel('Counts [a.u]', fontsize='x-large')

        ax.grid(which='both')
        ax.tick_params(direction='in',which='both')

        time, intensity = self.get_decay_curve(normalized=normalized)
        ax.plot(time, intensity)

        if v_lines != None:
            for value in v_lines:
                ax.axvline(value, ymax = 1, linestyle='--', color='black')
        if h_lines != None:
            for value in h_lines:
                ax.axhline(value, xmax = 1, linestyle='--', color='black')

        return fig, ax

    # UNDER CONSTRUCTION
    def normalize_and_get_lifetime(self, time_interval, n=20, plot = False):
        #Calculates the lifetime of a decay curve, in a defined time interval
        #time_interval: two-dimensional list
        #n: length of the list created backward from time_interval[1]
        #   the mean value os this list will be used as the baseline
        #plot: baseline corrected decay curve plot

        (time, intensity) = self.get_decay_curve(normalized=False)

        index_of_separations = []
        for value in time_interval:
            nearest_time_value = min(time, key=lambda x: abs(x-value))
            index_of_separations.append(time.index(nearest_time_value))

        intensity_to_integrate = intensity[index_of_separations[0]:index_of_separations[1]]
        time_to_integrate = time[index_of_separations[0]:index_of_separations[1]]
        time_to_integrate = [t - min(time_to_integrate) for t in time_to_integrate]

        baseline_list = intensity_to_integrate[:-n:-1] #reverse. length n-1 backward
        self.baseline = sum(baseline_list)/len(baseline_list)
        intensity_bsl_corrected = [I-self.baseline for I in intensity_to_integrate]

        area = integrate.trapezoid(intensity_bsl_corrected, time_to_integrate)

        self.lifetime = area/intensity_bsl_corrected[0] #Method of area under curve

        if plot == True:
            self._plot_lftm_corrected(time_to_integrate, intensity_bsl_corrected)

        return(self.lifetime)

    def _plot_lftm_corrected(self, time, intensity, fitting_curve = None):
        #Plots the decay curves corrected (dark counts subtracted +
        #   initial decay time set to zero)

        #Setups for decay curves plots
        fig, ax = plt.subplots()
        ax.set_xlabel('Time [$\mu$s]', fontsize='x-large')
        ax.set_ylabel('Counts [a.u]', fontsize='x-large')

        ax.grid(which='both')
        ax.tick_params(direction='in',which='both')

        ax.plot(time, intensity)
        if fitting_curve != None:
            ax.plot(time, fitting_curve)

        ax.axhline(max(intensity)/(np.e), xmax = 0.7, linestyle='--', color='black')


class LifeTime():
# Evaluates te LifeTime vs. Temperature dependence
    def __init__(self, decay_curve_set, temperature_set, particle_name):
        self.decay_curve_set = decay_curve_set
        self.temperatures = temperature_set
        self.particle_name = particle_name

    def calculate(self, time_interval, n=20, normalized=False, plot=True):
    # Calculates and plot the LT vs. T dependence
        self.lifetimes = []
        for curve in self.decay_curve_set:
            self.lifetimes.append(curve.get_lifetime(time_interval, n=n, plot=False))

        if plot == True:
            self._plot()

    def _plot(self):
        #Plots the dependence of the LIR with the temperature
        fig, ax = plt.subplots(constrained_layout=True)
        ax.scatter(self.temperatures, self.lifetimes, color='000000')

        ax.set_xlabel('Temperature (K)$', size='x-large')
        ax.set_ylabel('Lifetime ($\mu$s)', size='x-large')
        ax.tick_params(direction='in',which='both')
        ax.grid()
        ax.set_title('Lifetime dependence - Particle '+self.particle_name, size='x-large')
