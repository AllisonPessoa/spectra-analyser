import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams.update({'errorbar.capsize': 2})

class Pot_Particula():
    def __init__(self,str_number,int_time,norm,plot):
        self.number = str_number
        self.int_time = int_time
        self.qtd_spectra = 0
        (self.verdes, self.vermelho) = self._loadFiles()
        self.pot = self._calculatePower()
        if norm == True:
            self._normalizeIntensity()
        else:
            pass
        self.printPowerDependence(535, plot)

    def _loadFiles(self):
        ##OPENS AND COMPILATE
        file_verdes = []
        file_vermelho = []

        for i in range(6):
            file_verdes.append(open('p'+self.number+'_verdes_OD_4_2+'+str(round(0.2*i,1)).replace('.','_')+'.asc'))
            file_vermelho.append(open('p'+self.number+'_vermelho_OD_4_2+'+str(round(0.2*i,1)).replace('.','_')+'.asc'))

        self.qtd_spectra = len(file_verdes)

        lambda_verdes = []
        int_verdes = []
        lambda_vermelho = []
        int_vermelho = []

        for i in range(self.qtd_spectra):
            lambda_verdes.append([]) #Creates a matrix [index temperature][wavelength]
            lambda_vermelho.append([])
            int_verdes.append([]) #Creates a matrix [index temperature][intensity]
            int_vermelho.append([])

            for line in file_verdes[i].readlines():
                lambda_verdes[i].append(float(line.split('\t')[0].replace(',','.')))
                int_verdes[i].append(float(line.split('\t')[1][:-1]))

            for line in file_vermelho[i].readlines():
                lambda_vermelho[i].append(float(line.split('\t')[0].replace(',','.')))
                int_vermelho[i].append(float(line.split('\t')[1][:-1]))

        for i in range(self.qtd_spectra):
            file_verdes[i].close()
            file_vermelho[i].close()

        return ([lambda_verdes,int_verdes],[lambda_vermelho,int_vermelho])

    def _normalizeIntensity(self):
        int_verdes = self.verdes[1]
        int_vermelho = self.vermelho [1]
        norm_verdes = []
        norm_vermelho = []
        for i in range(self.qtd_spectra):
            norm_verdes.append([]) #Creates a matrix [index power][wavelength]
            norm_vermelho.append([])
            for j in range(len(int_verdes[i])):
                norm_verdes[i].append((int_verdes[i][j])/(max(int_verdes[i])-min(int_verdes[i])))
                norm_vermelho[i].append((int_vermelho[i][j])/(max(int_verdes[i])-min(int_verdes[i])))

        (self.verdes, self.vermelho) = (([self.verdes[0],norm_verdes],[self.vermelho[0],norm_vermelho]))

    def _calculatePower(self):
        #CALCULATES THE POWER EXCITATION
        filtros = [4.2,4.4,4.6,4.8,5,5.2] #Ordering is important-> from largest to smallest power
        pot_0 = 2.9*10**(-6)/(10**(4.2-3)) #Power (W), because the power was measured with an OD3.0 filter and 4.2 is our reference
        pot = [(pot_0/(10**(i-filtros[0]))) for i in filtros]

        return pot

    def _calculate_integral(self,x,y): #x, y = [[arrays with x/y data]...[]]
        area = 0
        for i in range(len(y)-1):
            b1 = y[i+1]
            b2 = y[i]
            h = (x[i+1]-x[i])
            area += ((b1+b2)*h)/2
        return area

    def printSpectra(self):
        lambda_verdes = self.verdes[0]
        int_verdes = self.verdes[1]
        #PLOTS NORMALIZED SPECTRA
        fig, ax = plt.subplots()
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        ax.set_prop_cycle(color=CB_color_cycle)
        for i in range(self.qtd_spectra):
            ax.plot(lambda_verdes[i], int_verdes[i],  alpha=1-0.15*i, label='Power = {:.1e} W'.format(self.pot[i]))
        ax.set_xlim([510,535])
        ax.set_ylim([0,0.2])
        ax.set_title('p'+self.number)
        ax.set_ylabel('Intensity (a.u)', size='x-large')
        ax.set_xlabel('Wavelength (nm)', size='x-large')
        ax.grid(True)
        ax.legend()
        ax.tick_params(direction='in',which='both')
        #ax.set_yticklabels([])
        fig.set_size_inches(6,4)
        fig.savefig('p'+self.number+'_Max e Min Spectra.png', dpi = 600)

    def printPowerDependence(self,sep_verdes,plot):
        lambda_verdes = self.verdes[0]
        int_verdes = self.verdes[1]
        lambda_vermelho = self.vermelho[0]
        int_vermelho = self.vermelho[1]

        value_in_list_535 = min(lambda_verdes[0], key=lambda x:abs(x-sep_verdes))
        index_sep_535 = lambda_verdes[0].index(value_in_list_535)

        value_in_list_555 = min(lambda_verdes[0], key=lambda x:abs(x-555))
        index_sep_555 = lambda_verdes[0].index(value_in_list_555)

        self.integral_525=[]
        self.integral_547=[]
        self.integral_660=[]

        for i in range(self.qtd_spectra):
            self.integral_525.append(self._calculate_integral(lambda_verdes[i][:index_sep_535] , int_verdes[i][:index_sep_535]))
            self.integral_547.append(self._calculate_integral(lambda_verdes[i][index_sep_535:index_sep_555] , int_verdes[i][index_sep_535:index_sep_555]))
            self.integral_660.append(self._calculate_integral(lambda_vermelho[i] , int_vermelho[i]))

        #LINEAR FIT
        log_int_525 = [np.log10(a) for a in self.integral_525]
        log_int_547 = [np.log10(a) for a in self.integral_547]
        log_int_660 = [np.log10(a) for a in self.integral_660]
        log_pot = [np.log10(a) for a in self.pot]

        self.fit_525 = np.polyfit(log_pot, log_int_525, 1, full=True) #([a[0]*x + a[1]],residual error, rank, singular values, rcond)
        self.fit_547 = np.polyfit(log_pot, log_int_547, 1, full=True)
        self.fit_660 = np.polyfit(log_pot, log_int_660, 1, full=True)

        self.y_fit_525 = [10**(x*self.fit_525[0][0]+ self.fit_525[0][1]) for x in log_pot] #Int = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
        self.y_fit_547 = [10**(x*self.fit_547[0][0]+ self.fit_547[0][1]) for x in log_pot]
        self.y_fit_660 = [10**(x*self.fit_660[0][0]+ self.fit_660[0][1]) for x in log_pot]

        if plot==True:
            #PLOTS
            fig, ax = plt.subplots()
            ax.plot(self.pot, self.integral_525, color = '#4aff00', marker = 'o', linestyle='None', label="525 nm")
            ax.plot(self.pot, self.y_fit_525, color = '#4aff00', linestyle='--', label = "Slope: {:.1f} \u00B1 0.1".format(self.fit_525[0][0])+" $[W^{-1}]$")

            ax.plot(self.pot, self.integral_547, color = '#8fff00', marker = '^', linestyle='None', label="547 nm")
            ax.plot(self.pot, self.y_fit_547, color = '#8fff00', linestyle='--', label = "Slope: {:.1f} \u00B1 0.1".format(self.fit_547[0][0])+" $[W^{-1}]$")

            ax.plot(self.pot, self.integral_660, color = '#ff0000', marker = 'v', linestyle='None', label="660 nm")
            ax.plot(self.pot, self.y_fit_660, color = '#ff0000', linestyle='--', label = "Slope: {:.1f} \u00B1 0.1".format(self.fit_660[0][0])+" $[W^{-1}]$")

            plt.title('p'+self.number)
            plt.xlabel('Power Excitation [W]', fontsize='large')
            plt.ylabel('Signal Integral [a.u]', fontsize='large')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(labelspacing=0.1)
            plt.grid(which='both')

            ax.tick_params(direction='in',which='both')
            ax.set_xlim(9*10**(-9))

            plt.savefig('P'+self.number+' - Estudo de Potência.png', dpi=600)

# %%
class Pot_Estat_Espectra():
    def __init__(self,PTC):
        self.PTC = PTC
        self.calculate()

    def calculate(self):
        qtd_particles = len(self.PTC)
        #Media das integrais:
        todos_verdes = []
        todos_vermelho = []
        for p in self.PTC:
            todos_verdes.append(p.verdes[1])
            todos_vermelho.append(p.vermelho[1])

        qtd_part = len(todos_verdes)
        qtd_pot = len(todos_verdes[0])
        qtd_numbers = len(todos_verdes[0][0])

        media_verdes = []
        media_vermelho = []
        for k in range(qtd_pot):
            media_verdes.append([])
            media_vermelho.append([])
            for i in range(qtd_numbers):
                sum = 0
                sum2 = 0
                for j in range(qtd_particles):
                    sum += todos_verdes[j][k][i]
                    sum2 += todos_vermelho[j][k][i]
                media_verdes[k].append(sum/qtd_particles)
                media_vermelho[k].append(sum2/qtd_particles)

        self.plot_spectra_verdes(media_verdes, media_vermelho)
        #self.plot_spectra_vermelho(media_vermelho)

#  PLOTS
    def plot_spectra_verdes(self, media_verdes , media_vermelho):
        join_x = self.PTC[0].verdes[0][0]+self.PTC[0].vermelho[0][0]
        join_y_max = media_verdes[0]+media_vermelho[0]
        join_y_min = media_verdes[-1]+media_vermelho[-1]
        fig, ax = plt.subplots(constrained_layout=True)
        ax.plot(join_x, join_y_min, label='Power = {:.1e} W'.format(self.PTC[0].pot[-1]), color='b')
        ax.plot(join_x, join_y_max, label='Power = {:.1e} W'.format(self.PTC[0].pot[0]), color='r')

        ax.plot([535,535],[0,max(media_verdes[0])/2.1], '--', color='k')
        ax.plot([505,505],[0,max(media_verdes[0])/2.1], '--', color='k')
        ax.plot([555,555],[0,max(media_verdes[0])/2.1], '--', color='k')

        ax.plot([630,630],[0,max(media_verdes[0])/2.1], '--', color='k')
        ax.plot([680,680],[0,max(media_verdes[0])/2.1], '--', color='k')

        index = self.PTC[0].verdes[0][0].index(535.00842)
        #ax.fill_between(self.PTC[0].verdes[0][0][:index], media_verdes[0][:index], color = '#4aff00')
        #ax.fill_between(self.PTC[0].vermelho[0][0], media_vermelho[0], color = '#ff0000')
        #ax.fill_between(self.PTC[0].verdes[0][0][index:], media_verdes[0][index:], color = '#8fff00')

        ax.set_ylabel('Normalized Intensity', size=10)
        ax.set_ylim([0, max(media_verdes[0])*1.1])
        ax.set_xlabel('Wavelength [nm]',size=10)
        #ax.legend()

        ax.set_xticks(np.arange(500, 700, 20))
        ax.set_xticklabels(np.arange(500, 700, 20),rotation=45)
        ax.tick_params(direction='in',which='both', labelsize=7)

        ax.annotate(r'$^4H_{11/2}$',
                xy=(0.12, 0.23), xycoords='axes fraction',
                xytext=(0.03, 0.6),
                arrowprops={'arrowstyle': '->'}, va='center')
        ax.annotate(r'$^4S_{3/2}$',
                xy=(0.45, 0.65), xycoords='axes fraction',
                xytext=(0.22, 0.85))
        ax.annotate(r'$^4F_{9/2}$',
                xy=(0.75, 0.35), xycoords='axes fraction',
                xytext=(0.65, 0.6),
                arrowprops={'arrowstyle': '->'}, va='center')
        ax.annotate(r'$^2H_{9/2} \rightarrow ^4I_{13/2}$',
                xy=(0.30, 0.13), xycoords='axes fraction',
                xytext=(0.35, 0.35),
                arrowprops={'arrowstyle': '->'}, va='center')

        ax.annotate(r'$\lambda_0$', xy=(0.35, 0.46), xycoords='axes fraction',
                xytext=(0.02, 0.45))
        ax.annotate(r'$\lambda_1$', xy=(0.35, 0.46), xycoords='axes fraction',
                xytext=(0.15, 0.45))
        ax.annotate(r'$\lambda_2$', xy=(0.35, 0.46), xycoords='axes fraction',
                xytext=(0.27, 0.45))
        ax.annotate(r'$\lambda_3$', xy=(0.35, 0.46), xycoords='axes fraction',
                xytext=(0.60, 0.45))
        ax.annotate(r'$\lambda_4$', xy=(0.35, 0.46), xycoords='axes fraction',
                xytext=(0.90, 0.45))
        #ax.annotate('Max. Power', xy=(0.45, 0.2), xycoords='axes fraction')
        insec = False
        if insec == True:
            rect = [0.565,0.71,0.55,0.40]
            ax1 = self.add_subplot_axes(ax,rect)
            join_x0 = self.PTC[0].verdes[0][0] + self.PTC[0].vermelho[0][0]
            join_y0 = media_verdes[-1] + media_vermelho[-1]
            ax1.plot(join_x0, join_y0, 'k', linewidth=0.3)
            ax1.plot([535,535],[0,max(media_verdes[-1])/2], '--', color='k', linewidth=0.2)
            ax1.annotate('Min. Power',xy=(0.55, 0.6), xycoords='axes fraction')
            #ax1.fill_between(self.PTC[0].verdes[0][0][:index], media_verdes[-1][:index], color = '#4aff00')
            #ax1.fill_between(self.PTC[0].verdes[0][0][index:], media_verdes[-1][index:], color = '#8fff00')
            #ax1.fill_between(self.PTC[0].vermelho[0][0], media_vermelho[-1], color = '#ff0000')
            ax1.set_ylabel('Intensity [a.u]', size=7)
            ax1.set_xlabel('Wavelength [nm]', size=7)
            ax1.xaxis.set_label_position('top')
            ax1.xaxis.tick_top()
            ax1.tick_params(direction='in',which='both', labelsize=7)
            ax1.set_yticks(np.arange(0, 75, 25))

        fig.set_size_inches(3.54,1.97)
        fig.savefig('spectrum.pdf', dpi = 600)

    def add_subplot_axes(self, ax,rect,axisbg='w'):
        fig = plt.gcf()
        box = ax.get_position()
        width = box.width
        height = box.height
        inax_position  = ax.transAxes.transform(rect[0:2])
        transFigure = fig.transFigure.inverted()
        infig_position = transFigure.transform(inax_position)
        x = infig_position[0]
        y = infig_position[1]
        width *= rect[2]
        height *= rect[3]  # <= Typo was here
        subax = fig.add_axes([x,y,width,height])
        x_labelsize = subax.get_xticklabels()[0].get_size()
        y_labelsize = subax.get_yticklabels()[0].get_size()
        x_labelsize *= rect[2]**0.5
        y_labelsize *= rect[3]**0.5
        subax.xaxis.set_tick_params(labelsize=x_labelsize)
        subax.yaxis.set_tick_params(labelsize=y_labelsize)
        return subax

class Pot_Estat_Fit():
    def __init__(self,PTC):
        self.PTC = PTC
        self.calculate()
        #self.calculate_FIR()

    def calculate(self):
        qtd_particles = len(self.PTC )
        #Media das integrais:
        self.int_media_525 = []
        self.int_media_547 = []
        self.int_media_660 = []
        for i in range(self.PTC[0].qtd_spectra):
            self.int_media_525.append(0)
            self.int_media_547.append(0)
            self.int_media_660.append(0)
            for p in self.PTC:
                self.int_media_525[i] += p.integral_525[i]
                self.int_media_547[i] += p.integral_547[i]
                self.int_media_660[i] += p.integral_660[i]

        self.int_media_525 = [x/qtd_particles for x in self.int_media_525]
        self.int_media_547 = [x/qtd_particles for x in self.int_media_547]
        self.int_media_660 = [x/qtd_particles for x in self.int_media_660]

        self.int_std_525 = np.zeros(self.PTC[0].qtd_spectra)
        self.int_std_547 = np.zeros(self.PTC[0].qtd_spectra)
        self.int_std_660 = np.zeros(self.PTC[0].qtd_spectra)

        for i in range(self.PTC[0].qtd_spectra):
            for j in range(len(self.PTC)):
                self.int_std_525[i] = (self.PTC[j].integral_525[i] - self.int_media_525[i])**2
                self.int_std_547[i] = (self.PTC[j].integral_547[i] - self.int_media_547[i])**2
                self.int_std_660[i] = (self.PTC[j].integral_660[i] - self.int_media_660[i])**2
            self.int_std_525[i] = np.sqrt(self.int_std_525[i]/(len(self.PTC)-1))
            self.int_std_547[i] = np.sqrt(self.int_std_547[i]/(len(self.PTC)-1))
            self.int_std_660[i] = np.sqrt(self.int_std_660[i]/(len(self.PTC)-1))

            self.int_std_525[i] = self.int_std_525[i]/np.sqrt(len(self.PTC))
            self.int_std_547[i] = self.int_std_547[i]/np.sqrt(len(self.PTC))
            self.int_std_660[i] = self.int_std_660[i]/np.sqrt(len(self.PTC))

        #LINEAR FIT
        log_int_525 = [np.log10(a) for a in self.int_media_525]
        log_int_547 = [np.log10(a) for a in self.int_media_547]
        log_int_660 = [np.log10(a) for a in self.int_media_660]
        self.pot_POT = self.PTC[0].pot
        sigma = 350*10**(-7)#cm Gaussian beam standard deviation, RESPONSE TO REFEREE,
        self.pot = [x/(np.pi*(1.8*sigma)**2) for x in self.pot_POT] #power density RESPONSE TO REFEREE, 12º OCTOBER 2020
        log_pot = [np.log10(a) for a in self.pot]

        self.fit_525 = np.polyfit(log_pot, log_int_525, 1, full=True) #([a[0]*x + a[1]],residual error, rank, singular values, rcond)
        self.fit_547 = np.polyfit(log_pot, log_int_547, 1, full=True)
        self.fit_660 = np.polyfit(log_pot, log_int_660, 1, full=True)

        self.y_fit_525 = [10**(x*self.fit_525[0][0]+ self.fit_525[0][1]) for x in log_pot] #Int = 10^Y = 10^(AX + B) = 10^(A*log(pot)+B)
        self.y_fit_547 = [10**(x*self.fit_547[0][0]+ self.fit_547[0][1]) for x in log_pot]
        self.y_fit_660 = [10**(x*self.fit_660[0][0]+ self.fit_660[0][1]) for x in log_pot]

        #self.plot_indiv_line()
        self.plot_estat_pot()

    def calculate_FIR(self):
        FIR_mean = np.zeros(self.PTC[0].qtd_spectra)
        FIR_all = []
        std_fir = np.zeros(self.PTC[0].qtd_spectra)

        for j in range(self.PTC[0].qtd_spectra):
            sum = 0
            for i in range(len(self.PTC)):
                sum += (self.PTC[i].integral_525[j]/self.PTC[i].integral_547[j])
            FIR_mean[j] = sum/(len(self.PTC))

        for i in range(len(self.PTC)):
            FIR_all.append([])
            for j in range(self.PTC[0].qtd_spectra):
                FIR_all[i].append((self.PTC[i].integral_525[j]/self.PTC[i].integral_547[j]))

        for j in range(self.PTC[0].qtd_spectra):
            sum = 0
            for i in range(len(self.PTC)):
                sum += ((self.PTC[i].integral_525[j]/self.PTC[i].integral_547[j]) - FIR_mean[j])**2
            std_fir[j] = np.sqrt(sum/(len(self.PTC)-1))/np.sqrt(len(self.PTC))

        ln_FIR = [np.log(x) for x in FIR_mean]
        ln_FIR_std = [std_fir[i]/FIR_mean[i] for i in range(len(FIR_mean))]

        #Linear fit to calculate relatively sensibility
        # fit_fir = np.polyfit(self.pot, ln_FIR, 1, full=True)
        # fit_curve = [x*fit_fir[0][0]+fit_fir[0][1] for x in self.pot]

        fig, ax = plt.subplots(constrained_layout=True)

        ax.errorbar(self.pot, FIR_mean, std_fir, color='000000', linestyle='none', marker='o', ms=2)
        for i in range(len(self.PTC)):
            ax.scatter(self.pot, FIR_all[i])
        #ax.plot(self.pot, fit_curve, color='000000')

        ax.set_xlabel(r'Excitation Power [W]$', size=10)
        ax.set_ylabel('FIR', size=10)
        ax.tick_params(direction='in',which='both')

        ax.set_xlim([min(self.pot)*0.5,max(self.pot)*1.05])
        fig.set_size_inches(3.54,1.97)
        fig.savefig("fir_power_dependence.png", dpi=600)

    def plot_estat_pot(self):
        #PLOTS
        #self.pot = [x*10**8 for x in self.pot] ##12/10/2020
        fig, ax = plt.subplots(constrained_layout=True)
        ax.errorbar(self.pot, self.int_media_525, color = '#4aff00', marker = 'o', linestyle='None', ms=4)
        ax.plot(self.pot, self.y_fit_525, color = '#4aff00', linestyle='dotted', label = "$^4H_{11/2} → ^4I_{15/2}$ (525 nm)")
        # ax.annotate('$^4H_{11/2}$',
        #         xy=(0.23, 0.05), xycoords='axes fraction',
        #         xytext=(0.03, 0.05),
        #         arrowprops={'arrowstyle': '->'}, va='center')
        ax.errorbar(self.pot, self.int_media_547, color = '#8fff00', marker = 's', linestyle='None', ms=4)
        ax.plot(self.pot, self.y_fit_547, color = '#8fff00', linestyle='-.', label = "$^4S_{3/2} → ^4I_{15/2}$ (547 nm)")
        # ax.annotate('$^4S_{3/2}$',
        #         xy=(0.23, 0.28), xycoords='axes fraction',
        #         xytext=(0.05, 0.35),
        #         arrowprops={'arrowstyle': '->'}, va='center')
        ax.errorbar(self.pot, self.int_media_660, color = '#ff0000', marker = 'v', linestyle='None', ms=4)
        ax.plot(self.pot, self.y_fit_660, color = '#ff0000', linestyle='--', label = "$^4F_{9/2} → ^4I_{15/2}$ (660 nm)")
        # ax.annotate('$^4F_{9/2}$',
        #         xy=(0.22, 0.20), xycoords='axes fraction',
        #         xytext=(0.02, 0.2),
        #         arrowprops={'arrowstyle': '->'}, va='center')

        ax.set_xlabel('Excitation Power Density [$\mathrm{W/cm^2}$]', fontsize=10)
        ax.set_ylabel('Integrated Intensity [arb.units]', fontsize=10)
        ax.set_xscale('log')
        ax.set_yscale('log')
        #ax.legend(fontsize=7, loc='upper left')
        #ax.grid(which='both')

        ax.tick_params(direction='in',which='both',labelsize=7)
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
        #ax.set_xlim(1*10**(-8))
        fig.set_size_inches(3.54,2.52) #7:5 scale in single column
        fig.savefig('power_dependence.pdf', dpi=600)

    def plot_indiv_line(self):
        fig = [0,0,0]
        axs = [0,0,0]
        marker=['o', 's', '^', '*', 'd']
        color = ['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628']

        fig[0], axs[0] = plt.subplots()
        for j in range(len(self.PTC)):
            axs[0].plot(self.PTC[j].pot, self.PTC[j].integral_525, color = color[j], marker = marker[j], linestyle='None',
                        label='P'+str(int(self.PTC[j].number)-1).zfill(1)+ r': $\alpha = $'+'{:.1f}(1)'.format(self.PTC[j].fit_525[0][0])+" $[W^{-1}]$")
            axs[0].plot(self.PTC[j].pot, self.PTC[j].y_fit_525, color = color[j], linestyle='--')
            #bbox_props = dict(boxstyle="round,pad=0.5", fc="#4aff00", ec="k", lw=1)
            #axs[0].annotate('525 nm', xy=(0.55, 0.8), xycoords='axes fraction',bbox=bbox_props)

        fig[1], axs[1] = plt.subplots()
        for j in range(len(self.PTC)):
            axs[1].plot(self.PTC[j].pot, self.PTC[j].integral_547, color = color[j], marker = marker[j], linestyle='None',
                        label='P'+str(int(self.PTC[j].number)-1).zfill(1)+ r': $\alpha = $'+'{:.1f}(1)'.format(self.PTC[j].fit_547[0][0])+" $[W^{-1}]$")
            axs[1].plot(self.PTC[j].pot, self.PTC[j].y_fit_547, color = color[j], linestyle='--')
            #bbox_props = dict(boxstyle="round,pad=0.5", fc="#8fff00", ec="k", lw=1)
            #axs[1].annotate('547 nm', xy=(0.55, 0.8), xycoords='axes fraction',bbox=bbox_props)

        fig[2], axs[2] = plt.subplots()
        for j in range(len(self.PTC)):
            axs[2].plot(self.PTC[j].pot, self.PTC[j].integral_660, color = color[j], marker = marker[j], linestyle='None',
                        label='P'+str(int(self.PTC[j].number)-1).zfill(1)+ r': $\alpha = $'+'{:.1f}(1)'.format(self.PTC[j].fit_660[0][0])+" $[W^{-1}]$")
            axs[2].plot(self.PTC[j].pot, self.PTC[j].y_fit_660, color = color[j], linestyle='--')
            #bbox_props = dict(boxstyle="round,pad=0.5", fc="#ff0000", ec="k", lw=1)
            #axs[2].annotate('660 nm', xy=(0.55, 0.8), xycoords='axes fraction',bbox=bbox_props)

        for i in range(len(axs)):
            x=['525 nm', '547 nm', '660 nm']
            axs[i].set_xlabel('Power Excitation [W]', fontsize='x-large')
            axs[i].set_ylabel('Signal Integral [a.u]', fontsize='x-large')
            axs[i].set_xscale('log')
            axs[i].set_yscale('log')
            axs[i].legend(title=x[i], labelspacing=0.2, loc='upper-right')
            axs[i].grid(which='both')
            axs[i].tick_params(direction='in',which='both', labelsize=12)
            axs[i].set_xlim(7*10**(-9))

        for i in range(len(fig)):
            x=[525,544,660]
            fig[i].set_size_inches(6,4.5)
            fig[i].savefig('part_juntas_'+str(x[i])+'.png', dpi=600)

PTC = []
PTC.append(Pot_Particula('02',5,norm=False,plot=False))
PTC.append(Pot_Particula('03',5,norm=False,plot=False))
PTC.append(Pot_Particula('04',2,norm=False,plot=False))
PTC.append(Pot_Particula('05',2,norm=False,plot=False))
PTC.append(Pot_Particula('06',2,norm=False,plot=False))

#Est_Espectra = Pot_Estat_Espectra(PTC)
Est_Potencia = Pot_Estat_Fit(PTC)

# %%
# fig, ax = plt.subplots()
# for i in range(len(PTC)):
#     ax.plot(PTC[0].pot, PTC[i].integral_547, label="P"+str(i))
# ax.legend()
# ax.set_xlabel("Potência [W]")
# ax.set_ylabel("Intensidade [a.u]")
# ax.set_title("Emissão 547 nm", size='x-large')
# fig.savefig('emissao_547_pot.png')
#
# fig, ax = plt.subplots()
# for i in range(len(PTC)):
#     ax.plot(PTC[0].pot, PTC[i].integral_525, label="P"+str(i))
# ax.legend()
# ax.set_xlabel("Potência [W]")
# ax.set_ylabel("Intensidade [a.u]")
# ax.set_title("Emissão 525 nm", size='x-large')
# fig.savefig('emissao_525_pot.png')
