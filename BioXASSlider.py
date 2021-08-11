from flask import Flask, render_template,request, send_file
import plotly
import plotly.graph_objs as go

#import pandas as pd
import numpy as np
import json
import time
import os, sys;
#sys.path.append(r'/home/ubuntu/xrt-1.3.4')  # analysis:ignore
#sys.path.append(r'c:\github\xrt')  # analysis:ignore
sys.path.append(r'D:\xrt-1.3.5')  # analysis:ignore
import xrt.backends.raycing.materials as rmats
import xrt.backends.raycing.sources as rsources
import copy

energyRange = np.linspace(5000, 30000, 501)

CVD = rmats.Material(
    elements=r"C",
    kind=r"plate",
    rho=3.52,
    name=r"Diamond plate")

Rh = rmats.Material(
    elements=r"Rh",
    kind=r"mirror",
    rho=12.41,
    name=r"Rhodium bulk")

Be = rmats.Material(
    elements=r"Be",
    kind=r"plate",
    rho=1.85,
    name=r"Beryllium window")

Al = rmats.Material(
    elements=r"Al",
    kind=r"plate",
    rho=2.7,
    name=r"Aluminum foil")

Kapton = rmats.Material(
    elements=('C', 'H', 'N', 'O'),
    quantities=(35, 28, 2, 7),
    kind=r"plate",
    rho=1.42,
    name=r"Kapton film")

Air = rmats.Material(
    elements=('N', 'O'),
    quantities=(0.8, 0.2),
    kind=r"plate",
    rho=1.225e-3,
    name=r"Air")

N2 = rmats.Material(
    elements=('N',),
    quantities=(2,),
    kind=r"plate",
    rho=1.25e-3,
    name=r"N2")

Ar = rmats.Material(
    elements=('Ar',),
    kind=r"plate",
    rho=1.784e-3,
    name=r"Ar")

He = rmats.Material(
    elements=('He'),
    kind=r"plate",
    rho=0.179e-3,
    name=r"He")

Si220 = rmats.CrystalSi(
    hkl=[2, 2, 0],
    name=r"Si220")

CVDcoating = rmats.Material(
    elements=r"C",
    kind=r"mirror",
    rho=3.52,
    name=r"CVD coating")

Si = rmats.Material(
    elements=r"Si",
    kind=r"mirror",
    rho=2.33,
    name=r"Si bulk")

RhOnSi = rmats.Coated(
    coating=Rh,
    cThickness=930,  # Angstroem
    surfaceRoughness=3.5,  # Angstroem
    substrate=Si,
    substRoughness=0.5,
    name=r"Rhodium on Silicon")

CVDonSi = rmats.Coated(
    coating=CVDcoating,
    cThickness=300,
    surfaceRoughness=5,
    substrate=Si,
    substRoughness=5,
    name=r"CVD on Silicon")

Si220harm = rmats.CrystalHarmonics(
    Nmax=2,
    name=r"Si220 with harmonics",
    hkl=[2, 2, 0],
    a=5.41949,
    tK=297.15)

thetaRange = np.linspace(-6/12000, 6/12000, 101)
psiRange = np.linspace(-1.05/12000, 1.05/12000, 101)
Wiggler = rsources.Wiggler(
        bl=None,
        name=r"wiggler",
        center=[0, 0, 0],
        yaw=-0.005,
        nrays=500000,
        eE=2.9,
        eI=0.25,
        eEspread=0.001,
        eEpsilonX=18.1,
        eEpsilonZ=0.0362,
        betaX=9.1,
        betaZ=2.8,
        xPrimeMax=6.1,
        zPrimeMax=0.175,
        eMin=energyRange[0],
        eMax=energyRange[-1],
        distE='BW',
        K=35,
        period=150,
        n=11)

class TransmissionMatrix(object):

    def __init__(self):
        nElements = 11
        self.h1Matrix = np.ones((nElements, len(energyRange)))
        self.h2Matrix = np.ones((nElements, len(energyRange)))
        self.energyRange = energyRange
        er2 = energyRange*2
        try:
            self.h1Matrix[0, :] = np.loadtxt('wigglerFluxBW.dat') * 0.1
        except:
            I0, l1, l2, l3 = Wiggler.intensities_on_mesh(
                    energyRange, thetaRange, psiRange)
            self.h1Matrix[0, :] = self.flux_through_aperture(
                    thetaRange, psiRange, I0 * 0.1)
            np.savetxt('wigglerFluxBW.dat', self.h1Matrix[0, :])

        self.tCFF1 = [0, 50]
        self.tCFF2 = [0, 700]

        self.CVDabs = CVD.get_absorption_coefficient(energyRange)
        self.h1Matrix[1, :] = np.abs(np.exp(-(self.tCFF1[0])*1e-4*self.CVDabs))
        self.h1Matrix[2, :] = np.abs(np.exp(-(self.tCFF2[0])*1e-4*self.CVDabs))

        mirrPitch = 0.16
        mirrAbs = abs(RhOnSi.get_amplitude(energyRange, np.sin(np.radians(mirrPitch)))[0])**4
        self.h1Matrix[3, :] = mirrAbs

        dbhrPitch = 0.20
        if dbhrPitch > 0:
            dbhrAbs = abs(CVDonSi.get_amplitude(energyRange, np.sin(np.radians(dbhrPitch)))[0])**4
            self.h1Matrix[4, :] = dbhrAbs

        self.h1Matrix[5, :] = np.abs(np.exp(-50e-4 * Be.get_absorption_coefficient(energyRange)))
        self.h1Matrix[6, :] = np.abs(np.exp(-50e-4 * 4 * Kapton.get_absorption_coefficient(energyRange)))

        self.AirAbs = Air.get_absorption_coefficient(energyRange)
        self.AirAbs2 = Air.get_absorption_coefficient(er2)
        self.HeAbs = He.get_absorption_coefficient(energyRange)
        self.HeAbs2 = He.get_absorption_coefficient(er2)
        self.h1Matrix[7, :] = np.abs(np.exp(-6 * self.AirAbs))
        self.h1Matrix[8, :] = np.abs(np.exp(-70 * self.AirAbs))
        self.AlAbs = Al.get_absorption_coefficient(energyRange)
        self.AlAbs2 = Al.get_absorption_coefficient(er2)
        self.N2Abs = N2.get_absorption_coefficient(energyRange)
        self.N2Abs2 = N2.get_absorption_coefficient(er2)
        self.ArAbs = Ar.get_absorption_coefficient(energyRange)
        self.ArAbs2 = Ar.get_absorption_coefficient(er2)

        try:
            self.h2Matrix[0, :] = np.loadtxt('wigglerFluxBW2.dat') * 0.1
        except:
            I0, l1, l2, l3 = Wiggler.intensities_on_mesh(
                    er2, thetaRange, psiRange)
            self.h2Matrix[0, :] = self.flux_through_aperture(
                    thetaRange, psiRange, I0 * 0.1)
            np.savetxt('wigglerFluxBW2.dat', self.h2Matrix[0, :])

        self.CVDabs2 = CVD.get_absorption_coefficient(er2)
        self.h2Matrix[1, :] = np.abs(np.exp(-(self.tCFF1[0])*1e-4*self.CVDabs2))
        self.h2Matrix[2, :] = np.abs(np.exp(-(self.tCFF2[0])*1e-4*self.CVDabs2))

        mirrAbs = abs(RhOnSi.get_amplitude(er2, np.sin(np.radians(mirrPitch)))[0])**4
        self.h2Matrix[3, :] = mirrAbs

        if dbhrPitch > 0:
            dbhrAbs = abs(CVDonSi.get_amplitude(er2, np.sin(np.radians(dbhrPitch)))[0])**4
            self.h2Matrix[4, :] = dbhrAbs

        self.h2Matrix[5, :] = np.abs(np.exp(-50e-4 * Be.get_absorption_coefficient(er2)))
        self.h2Matrix[6, :] = np.abs(np.exp(-50e-4 * 4 * Kapton.get_absorption_coefficient(er2)))
        self.h2Matrix[7, :] = np.abs(np.exp(-6 * self.AirAbs2))
        self.h2Matrix[8, :] = np.abs(np.exp(-70 * self.AirAbs2))
#        self.h2Matrix[10, :] = 1. - np.abs(np.exp(-12 * self.N2signal2))

    def flux_through_aperture(self, theta, psi, I0):
        dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]
        flux = I0.sum(axis=(1, 2)) * dtheta * dpsi
        return flux

    def _recalc_CFF1(self, nr):
        print("NR", nr, file=sys.stderr)
        intId = int(nr)
        self.h1Matrix[1, :] = np.abs(np.exp(-self.tCFF1[intId]*1e-4*self.CVDabs))
        self.h2Matrix[1, :] = np.abs(np.exp(-self.tCFF1[intId]*1e-4*self.CVDabs2))

    def _recalc_CFF2(self, nr):
        intId = int(nr)
        self.h1Matrix[2, :] = np.abs(np.exp(-self.tCFF2[intId]*1e-4*self.CVDabs))
        self.h2Matrix[2, :] = np.abs(np.exp(-self.tCFF2[intId]*1e-4*self.CVDabs2))

    def _recalc_DBHR_Gas(self, nr):
        intId = int(nr)
        if intId == 0:
            absTable = np.zeros(len(self.energyRange))
            absTable2 = np.zeros(len(self.energyRange))
        elif intId == 1:
            absTable = self.HeAbs
            absTable2 = self.HeAbs2
        else:
            absTable = self.AirAbs
            absTable2 = self.AirAbs2
        self.h1Matrix[8, :] = np.abs(np.exp(-70 * absTable))
        self.h2Matrix[8, :] = np.abs(np.exp(-70 * absTable2))

#    def _recalc_DBHR_IO(self, nr):
#        intId = self.DBHRIOButtonGroup.id(nr)
#        if intId == 1:
#            pitch = 0
#            self.dbhrSlider.setEnabled(False)
#        else:
#            pitch = self.dbhrSlider.value()
#            self.dbhrSlider.setEnabled(True)
#        self._recalc_dbhr(pitch)

    def _recalc_XIA(self, nr):
#        t = 0
#        for button in self.XIAButtonGroup.buttons():
#            if button.isChecked():
#                t += float(button.text())
        t = float(nr)*250.
        self.h1Matrix[9, :] = np.abs(np.exp(-t*1e-4*self.AlAbs))
        self.h2Matrix[9, :] = np.abs(np.exp(-t*1e-4*self.AlAbs2))        
#        self._update_canvas()

    def _recalc_M1(self, pitch):
#        self.sliderQE1.setText(str(pitch))
        sinpitch = np.sin(np.radians(float(pitch)))
        self.h1Matrix[3, :] = abs(RhOnSi.get_amplitude(
                energyRange, sinpitch)[0])**4
        self.h2Matrix[3, :] = abs(RhOnSi.get_amplitude(
                energyRange*2, sinpitch)[0])**4

    def _recalc_DBHR(self, pitch, dbhr_io):
        if int(dbhr_io) == 0:
#            self.sliderQE2.setText(str(pitch))
#            self.dbhrPitch = pitch
            sinpitch = np.sin(np.radians(float(pitch)))
            self.h1Matrix[4, :] = abs(CVDonSi.get_amplitude(
                    energyRange, sinpitch)[0])**2 #4 one mirror temporary
            self.h2Matrix[4, :] = abs(CVDonSi.get_amplitude(
                    energyRange*2, sinpitch)[0])**2 #4 one mirror temporary
        else:
            self.h1Matrix[4, :] = np.ones(len(energyRange))
            self.h2Matrix[4, :] = np.ones(len(energyRange))

    def _recalc_IC_Gas(self, nr):
        intId = int(nr)
        if intId == 0:
#            absTable = np.zeros(len(self.energyRange))
#            absTable2 = np.zeros(len(self.energyRange))
            IEeV = 1
#        elif intId == 1:
#            absTable = self.HeAbs
#            absTable2 = self.HeAbs2
#            IEeV = 41
        elif intId == 1:  # 2:
            absTable = self.N2Abs
            absTable2 = self.N2Abs2
            IEeV = 36
        elif intId == 2: #3:
            absTable = self.AirAbs
            absTable2 = self.AirAbs2
            IEeV = 34.4
#        else:
#            absTable = self.ArAbs
#            absTable2 = self.ArAbs2
#            IEeV = 26
        if IEeV > 1:
            self.h1Matrix[10, :] = self.energyRange * (1. - np.abs(np.exp(-12 * absTable))) * 1.6e-19 / IEeV
            self.h2Matrix[10, :] = self.energyRange * (1. - np.abs(np.exp(-12 * absTable2)))* 1.6e-19 / IEeV
        else:
            self.h1Matrix[10, :] = np.ones(len(self.energyRange))
            self.h2Matrix[10, :] = np.ones(len(self.energyRange))
#        self._update_canvas()

#    def _update_canvas(self):
#        flux1 = np.prod(self.h1Matrix, axis=0)
#        self.fluxPlot.setData(x=self.energyRange,
#                              y=flux1)
#        self.h2plot.setData(x=self.energyRange,
#                            y=np.prod(self.h2Matrix, axis=0)/flux1)

    def update_matrices(self, params):
#        print(params)
        self._recalc_CFF1(params['cff1'][0])
        self._recalc_CFF2(params['cff2'][0])
        self._recalc_DBHR_Gas(params['dbhrgas'][0])
#        self._recalc_DBHR_IO()
        self._recalc_XIA(params['xia'][0])
        self._recalc_M1(params['m1p'][0])
        self._recalc_DBHR(params['dbhrp'][0], params['dbhrio'][0])
        self._recalc_IC_Gas(params['i0'][0])

default_tm = TransmissionMatrix()

kwargs_BL = {'cff1': '0', 
          'cff2': '0', 
          'dbhrgas': '2', 
          'dbhrio': '1', 
          'i0': '0', 
          'm1p': '0.15', 
          'dbhrp': '0.2', 
          'xia': '0'}

kwargs_calc = {'geometry': "Bragg",
          'material': "Si",
          'energy': 9000,
          'hkl_h': 1,
          'hkl_k': 1,
          'hkl_l': 1,
          'axis_min': -100,
          'axis_max': 100,
          'thickness': 1.,
          'asymmetry': 0}

def calc_vectors(theta, alphaDeg, geom):
    alpha = np.radians(alphaDeg)
    s0 = (np.zeros_like(theta), np.cos(theta+alpha), -np.sin(theta+alpha))
    sh = (np.zeros_like(theta), np.cos(theta-alpha), np.sin(theta-alpha))
    if geom.startswith('Bragg'):
        n = (0, 0, 1)  # outward surface normal
    else:
        n = (0, -1, 0)  # outward surface normal
    hn = (0, np.sin(alpha), np.cos(alpha))  # outward Bragg normal
    gamma0 = sum(i*j for i, j in zip(n, s0))
    gammah = sum(i*j for i, j in zip(n, sh))
    hns0 = sum(i*j for i, j in zip(hn, s0))
    return gamma0, gammah, hns0

app = Flask(__name__)

#def calc_vectors(theta, alphaDeg, geom):
#    alpha = np.radians(alphaDeg)
#    s0 = (np.zeros_like(theta), np.cos(theta+alpha), -np.sin(theta+alpha))
#    sh = (np.zeros_like(theta), np.cos(theta-alpha), np.sin(theta-alpha))
#    if geom.startswith('Bragg'):
#        n = (0, 0, 1)  # outward surface normal
#    else:
#        n = (0, -1, 0)  # outward surface normal
#    hn = (0, np.sin(alpha), np.cos(alpha))  # outward Bragg normal
#    gamma0 = sum(i*j for i, j in zip(n, s0))
#    gammah = sum(i*j for i, j in zip(n, sh))
#    hns0 = sum(i*j for i, j in zip(hn, s0))
#    return gamma0, gammah, hns0
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calc')
def index_calc():
    return render_template('MPtableinput.html', plot=create_plot_calc(**kwargs_calc))

def create_plot_calc(**kwargs_in):
    hkl = [int(kwargs_in['hkl_h']),
           int(kwargs_in['hkl_k']),
           int(kwargs_in['hkl_l'])]
    geom = str(kwargs_in['geometry'])
    mat = str(kwargs_in['material'])
    if mat == 'Si':
        crystal = rmats.CrystalSi(hkl=hkl,
                               geom=geom + ' reflected',
                               t=float(kwargs_in['thickness']))
    else:  # Ge
        crystal = rmats.CrystalDiamond(hkl=hkl, d=5.657, elements='Ge',
                               geom=geom + ' reflected',
                               t=float(kwargs_in['thickness']))
    E = float(kwargs_in['energy'])
    dtheta = np.linspace(float(kwargs_in['axis_min']), 
                         float(kwargs_in['axis_max']), 501)
    theta = crystal.get_Bragg_angle(E) + dtheta*1e-6
    asymmDeg = float(kwargs_in['asymmetry'])  # Degrees

    g0, gh, hs0 = calc_vectors(theta, asymmDeg, geom)
    curS, curP = crystal.get_amplitude(E, g0, gh, hs0)
    data = [go.Scatter(x=dtheta, y=abs(curS)**2, mode='lines', name='s', line=dict(color='rgb(220, 0, 0)', width=2)),
            go.Scatter(x=dtheta, y=abs(curP)**2, mode='lines', name='p', line=dict(color='rgb(0, 0, 220)', width=2, dash='dash'))]
#    data = [go.Line(x=dtheta, y=abs(curS)**2)]
#    print(kwargs_in)
    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/upd_calc', methods=['GET', 'POST'])
def change_param_calc():
#    print(request.args)
#    param_name = request.args['pname']
#    param_value = request.args['pvalue']
#    print(request.args)
#    kwargs = np.copy(default_kwargs)
    localargs = dict(request.args)
#    print("LA", localargs, file=sys.stderr)
#    print(localargs)
#    kwargs[param_name] = param_value
#    print(kwargs)
    graphJSON = create_plot_calc(**localargs)
    return graphJSON

@app.route('/BioXAS')
def index_slider():
    return render_template('BioXAStableinput.html', plot=create_plot(**kwargs_BL))

def create_plot(**kwargs_in):
    print("KW", kwargs_in, file=sys.stderr)
    time0 = time.time()
    locTM = copy.deepcopy(default_tm)
    locTM.update_matrices(kwargs_in)
    
    flux1 = np.prod(locTM.h1Matrix, axis=0)
#    flux2 = np.prod(locTM.h2Matrix, axis=0)/flux1
    time1 = time.time()
#    print("matrix updated in", time1-time0, "s")

    axname = 'flux' if int(kwargs_in['i0'][0]) == 0 else 'counts/s'
    data = [go.Scatter(x=locTM.energyRange, y=flux1, mode='lines', name=axname,
                       line=dict(color='rgb(0, 0, 220)', width=2))] #,
#            go.Scatter(x=locTM.energyRange, y=flux2, mode='lines', name='2nd harmonic', line=dict(color='rgb(220, 0, 0)', width=2))]

    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/upd_bioxas', methods=['GET', 'POST'])
def change_param():
    localargs = dict(request.args)
    graphJSON = create_plot(**localargs)
    return graphJSON


if __name__ == '__main__':
    app.run(host="10.45.0.114", port="5080", debug=True)
#    app.run(host="192.168.2.196", port="5080")
