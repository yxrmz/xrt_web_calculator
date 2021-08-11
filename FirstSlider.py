from flask import Flask, render_template,request
import plotly
import plotly.graph_objs as go

import pandas as pd
import numpy as np
import json

import os, sys;
sys.path.append(r'G:/xrt-1.3.4')  # analysis:ignore
import xrt.backends.raycing.materials as rm

kwargs = {'geometry': "Bragg",
          'material': "Si",
          'energy': 9000,
          'hkl_h': 1,
          'hkl_k': 1,
          'hkl_l': 1,
          'thickness': 1.,
          'asymmetry': 0}

app = Flask(__name__)

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


@app.route('/')
def index_slider():
    return render_template('tableinput.html', plot=create_plot(**kwargs))

def create_plot(**kwargs_in):
    hkl = [int(kwargs_in['hkl_h']),
           int(kwargs_in['hkl_k']),
           int(kwargs_in['hkl_l'])]
    geom = str(kwargs_in['geometry'])
    crystal = rm.CrystalSi(hkl=hkl,
                           geom=geom + ' reflected',
                           t=float(kwargs_in['thickness']))
    E = float(kwargs_in['energy'])
    dtheta = np.linspace(-100, 100, 501)
    theta = crystal.get_Bragg_angle(E) + dtheta*1e-6
    asymmDeg = float(kwargs_in['asymmetry'])  # Degrees

    g0, gh, hs0 = calc_vectors(theta, asymmDeg, geom)
    curS, curP = crystal.get_amplitude(E, g0, gh, hs0)
    data = [go.Line(x=dtheta, y=abs(curS)**2)]
#    print(kwargs_in)
    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

#@app.route('/recalc', methods=['GET', 'POST'])
#def change_features():
#    sfreq = request.args['sfreq']
#    print(sfreq)
#    kwargs = {'freq': sfreq}
#    graphJSON= create_plot(**kwargs)
#
#    return graphJSON

@app.route('/recalc', methods=['GET', 'POST'])
def change_param():
    param_name = request.args['pname']
    param_value = request.args['pvalue']
#    print(request.args)
#    kwargs = np.copy(default_kwargs)
    global kwargs
    kwargs[param_name] = param_value
#    print(kwargs)
    graphJSON = create_plot(**kwargs)
    return graphJSON


if __name__ == '__main__':
    app.run(host="192.168.2.119", port="5080")
