# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 21:47:01 2020

@author: roman
"""

from flask import Flask

app = Flask(__name__)


@app.route('/')
def index():
    return 'Hello World!'


if __name__ == '__main__':
    app.run()