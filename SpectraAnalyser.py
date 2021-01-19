####################################################
#               Spectrum analyser
# Python script used for analyzing spectral data
#
# Used and tested with data from iDus 401
# (obtained through Andor SOLIS for Spectroscopy)

# By Allison Pessoa, Nano-optics laboratory.
# UFPE, Brazil.
####################################################

import numpy as np

import scipy.integrate as integrate

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class Plots():
# Standard plot setup class
    @staticmethod
    def get_spectrum_baseplot():
        #Setups for spectrum plot
        fig, ax = plt.subplots()
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        ax.set_prop_cycle(color=CB_color_cycle)
        ax.set_ylabel('Intensity (a.u)', size='x-large')
        ax.set_xlabel('Wavelength (nm)', size='x-large')
        ax.grid(True)
        ax.tick_params(direction='in',which='both')

        return fig, ax

    @staticmethod
    def get_power_baseplot():
        #Setups for power-law dependency plots
        fig, ax = plt.subplots()
        ax.set_xlabel('Power Excitation [W]', fontsize='large')
        ax.set_ylabel('Signal Integral [a.u]', fontsize='large')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(which='both')
        ax.tick_params(direction='in',which='both')

        return fig, ax


class Spectrum():
# Contains all information about each individual spectrum
    def __init__(self, file_name, plot = True):
        #file_name: .asc file spectral data. File must contain two columns.
        #           first one is wavelength (decimal sep as ',')
        #           second one is intensity (counts/second, decimal sep as '.')
        #plot: Show (True) or not the spectrum after initializing the object

        self.filename = file_name
        (self.wavelength, self.intensity) = self._load_file()

    def _load_file(self):
        spectrum_file = open(self.filename)
        wavelength = []
        intensity = []

        for line in spectrum_file.readlines():
            wavelength.append(float(line.split('\t')[0].replace(',','.')))
            intensity.append(float(line.split('\t')[1][:-1]))

        spectrum_file.close()
        return (wavelength,intensity)

    def get_spectrum(self, normalized=False):
        #Returns a tuple with the spectral data: (wavelength, intensity)
        #If normalized is True, intensity will be normalized to values between 0 and 1
        intensity = []

        if normalized == True:
            for j in range(len(self.intensity)):
                intensity.append(self.intensity[j] + abs(min(self.intensity)))
                    #(self.intensity[j])/(max(self.intensity)-min(self.intensity)))
        else:
            intensity = self.intensity

        return (self.wavelength, intensity)

    def get_area_under_spectrum(self, integration_limits, normalized=False):
        #Returns a list with the areas under the spectrum plot between regions
        #                   defined by integration_limits
        #Integration_limits: list of wavelengths separating the data regions
        #                   to be integrated. Mininum two values
        #If normalized is True, intensity will be normalized to values between 0 and 1
        (wavelength, intensity) = self.get_spectrum(normalized=normalized)

        index_of_separations = []
        for value in integration_limits:
            nearest_wvlth_value = min(wavelength, key=lambda x: abs(x-value))
            index_of_separations.append(wavelength.index(nearest_wvlth_value))

        area = []
        index = index_of_separations
        for i in range(len(index_of_separations)-1):
            intensity_to_integrate = intensity[index[i]:index[i+1]]
            wavelength_to_integrate = wavelength[index[i]:index[i+1]]
            area.append(integrate.trapezoid(intensity_to_integrate,
                                            wavelength_to_integrate))

        return(area)

    def plot_spectrum(self, normalized = False):
        #Plots spectrum. Returns the matplotlib object for any alterations, if necessary.
        #If normalized is True, intensity will be normalized to values between 0 and 1
        fig, ax = Plots.get_spectrum_baseplot()
        wavelength, intensity = self.get_spectrum(normalized = normalized)
        ax.plot(wavelength, intensity)
        ax.set_title(self.filename)

        return(fig, ax)

class PowerDependence():
# Calculates the linear dependency between excitation power and intensity
# over a defined wavelength region
    def __init__(self, spectra_set, power_values):
        #spectra_set: list (length M) of Spectrum objects. Mininum two
        #power_values: list (length M) of sample excitation power in each spectrum
        #              (order is important).
        self.spectra_set = spectra_set
        self.power = power_values
        self.regions_label = None

    def add_spectrum(self, spectrum, power_value):
        #spectrum: Spectrum object
        #power_value: sample excitation power
        self.spectra_set.append(spectrum)
        self.power.append(power_value)

    def calculate_and_plot(self, integration_limits, regions_label=None, normalized=False):
        #Evaluates the linear fitting
        #Integration_limits: list of wavelengths separating the data regions
        #                   to be integrated. Mininum two values
        #regions_label: List of strings. Identifies each integration regions.
        #If normalized is True, intensity will be normalized to values between 0 and 1
        if regions_label == None:
            regions_label = range(len(integration_limits)-1)
        fig, ax = Plots.get_power_baseplot()

        areas = []#[[a,b,c], [c,d,e], ...] ; [a,b,c] = same power, different
                  #intergation limits
        for spectrum in self.spectra_set:
            areas.append(spectrum.get_area_under_spectrum(integration_limits,
                            normalized=normalized))

        log_areas = np.log10(areas)
        log_power = np.log10(self.power)

        for i in range(len(areas[0])): #for each region of interation
            linear_fit, cov_matrix = np.polyfit(log_power, log_areas[:,i], 1,
                            cov=True) #([a[0]*x + a[1]],residual error,...)
            fitting_line = [10**(x*linear_fit[0]+linear_fit[1]) for x in log_power]
                            #Int = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
            uncertainties = np.sqrt(np.diag(cov_matrix))

            ax.plot(self.power, np.array(areas)[:,i], color = '#4aff00',
                                marker = 'o', linestyle='None', label=regions_label[i])
            ax.plot(self.power, fitting_line, color = '#4aff00', linestyle='--',
                                label = "Slope: {:.1f} \u00B1 {:.1}"
                                .format(linear_fit[0], uncertainties[0])+" $[W^{-1}]$")
        ax.legend()

        self._plot_all(integration_limits, normalized=normalized)

    def _plot_all(self, integration_limits = None, normalized=False):
        #Plots all spectra in only one graph.
        fig, ax = Plots.get_spectrum_baseplot()
        for i in range(len(self.spectra_set)):
            wavelength, intensity = self.spectra_set[i].get_spectrum(normalized=normalized)
            ax.plot(wavelength,intensity, alpha=1-0.15*i,
                    label='Power = {:.1e} W'.format(self.power[i]))
        if integration_limits != None:
            for value in integration_limits:
                ax.axvline(value, ymax = 0.6, linestyle='--', color='black')
        ax.legend()

##### %%

particle01 = [Spectrum('Experiments/7_18-01-2021/0,11uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,18uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,27uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,34uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,40uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,50uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,65uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,87uW.asc'),
              Spectrum('Experiments/7_18-01-2021/1uW.asc'),
              Spectrum('Experiments/7_18-01-2021/2uW.asc'),
              Spectrum('Experiments/7_18-01-2021/3uW.asc'),
              Spectrum('Experiments/7_18-01-2021/4,3uW.asc'),
              Spectrum('Experiments/7_18-01-2021/5uW.asc')
             ]

power_particle01 = [x*10**(-6) for x in [0.11, 0.18, 0.27, 0.34, 0.4, 0.5, 0.65, 0.87, 1, 2, 3, 4.3, 5]]
power_dependence_particle01 = PowerDependence(particle01, power_particle01)
power_dependence_particle01.calculate_and_plot([535,555],regions_label=['544 nm'])

# %%

particle01 = [Spectrum('Experiments/7_18-01-2021/0,11uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,18uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,27uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,34uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,40uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,50uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,65uW.asc'),
              Spectrum('Experiments/7_18-01-2021/0,87uW.asc')
             ]

power_particle01 = [x*10**(-6) for x in [0.11, 0.18, 0.27, 0.34, 0.4, 0.5, 0.65, 0.87]]
power_dependence_particle01 = PowerDependence(particle01, power_particle01)
power_dependence_particle01.calculate_and_plot([535,555],regions_label=['544 nm'])


# %%
particle01 = [Spectrum('Experiments/6_16-01-2021/Particula1/0,3uW.asc'),
              Spectrum('Experiments/6_16-01-2021/Particula1/0,4uW.asc'),
              Spectrum('Experiments/6_16-01-2021/Particula1/0,5uW.asc'),
              Spectrum('Experiments/6_16-01-2021/Particula1/0,08uW.asc'),
              Spectrum('Experiments/6_16-01-2021/Particula1/0,13uW.asc'),
              Spectrum('Experiments/6_16-01-2021/Particula1/0,28uW.asc'),
              Spectrum('Experiments/6_16-01-2021/Particula1/0,48uW.asc')
             ]

power_particle01 = [x*10**(-6) for x in [0.3, 0.4, 0.5,0.08,0.13,0.28,0.48]]
power_dependence_particle01 = PowerDependence(particle01, power_particle01)
power_dependence_particle01.calculate_and_plot([535,555],regions_label=['544 nm'])

# %%
particle01 = [#Spectrum('Experiments/8_18-01-2021/0,08u.asc'),
              #Spectrum('Experiments/8_18-01-2021/0,20u.asc'),
              #Spectrum('Experiments/8_18-01-2021/0,30u.asc'),
              #Spectrum('Experiments/8_18-01-2021/0,40u.asc'),
              #Spectrum('Experiments/8_18-01-2021/0,51u.asc'),
              #Spectrum('Experiments/8_18-01-2021/0,60u.asc'),
              #Spectrum('Experiments/8_18-01-2021/0,70u.asc'),
              Spectrum('Experiments/8_18-01-2021/0,81u.asc'),
              Spectrum('Experiments/8_18-01-2021/1,0u.asc'),
              Spectrum('Experiments/8_18-01-2021/2,0u.asc'),
              Spectrum('Experiments/8_18-01-2021/2,9u.asc'),
              Spectrum('Experiments/8_18-01-2021/4,3u.asc')
             ]

power_particle01 = [x*10**(-6) for x in [0.81,1.0,2.0,2.9,4.3]]
power_dependence_particle01 = PowerDependence(particle01, power_particle01)
power_dependence_particle01.calculate_and_plot([517,535,554],regions_label=['525 nm','544 nm'])

# %%
particle01 = [
              Spectrum('Experiments/9_18-01-2021/0,1u.asc'),
              Spectrum('Experiments/9_18-01-2021/0,22u.asc'),
              Spectrum('Experiments/9_18-01-2021/0,31u.asc'),
              Spectrum('Experiments/9_18-01-2021/0,38u.asc'),
              Spectrum('Experiments/9_18-01-2021/0,48u.asc'),
              Spectrum('Experiments/9_18-01-2021/0,64u.asc')
             ]

power_particle01 = [x*10**(-6) for x in [0.1,0.22,0.31,0.38,0.48,0.64]]
power_dependence_particle01 = PowerDependence(particle01, power_particle01)
power_dependence_particle01.calculate_and_plot([516,535,570],regions_label=['525 nm','544 nm'])
