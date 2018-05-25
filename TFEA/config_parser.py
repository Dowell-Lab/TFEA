__author__ = 'Jonathan Rubin'

import configparser

def run(srcdirectory,configfile):
    config = configparser.ConfigParser(interpolation = configparser.ExtendedInterpolation())
    # print config.BasicInterpolation()
    config.read(configfile)
    outfile = open(srcdirectory+'config.py','w')
    for key in config:
        for item in config[key]:
            outfile.write(item.upper()+'='+config[key][item]+'\n')
