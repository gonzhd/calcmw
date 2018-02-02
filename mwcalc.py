# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import json
import re
import numpy as np

''' 
Atomic Number = 1
Atomic Symbol = H
Mass Number = 1
Relative Atomic Mass = 1.00782503223(9)
Isotopic Composition = 0.999885(70)
Standard Atomic Weight = [1.00784,1.00811]
'''

def load_db(fname='atomic_masses.txt'):
    f = open(fname)
    groups = f.read().split('\n\n')
    
    db = {}
    saw = {}
    for group in groups:
        dic = {}
        lines = group.split('\n')
        for line in lines:
            k,v = line.split(' = ')
            dic[k] = v
        atnum = dic['Atomic Number']
        atsym = dic['Atomic Symbol']
        print(atnum)
        db[atnum] = db.get(atnum, []) + [dic]
        if saw.get(atsym, False): 
            continue
        else:
            sawstr = dic['Standard Atomic Weight']
            if sawstr == ' ': 
                rng = None
            elif '[' in sawstr:
                vals = sawstr[1:-1].split(',')
                if len(vals) > 1:
                    rng = [float(vals[0]),float(vals[1])]
                else:
                    rng = [float(vals[0]),float(vals[0])]
            else:
                val = float(sawstr.split('(')[0])
                rng = [val,val]
            saw[atsym] = rng
    return db, saw
    
def make_atomic_weight_db(db):
    atdb = {}
    for k,v in db.items():
        print (k, v)
        atsym = v[0]['Atomic Symbol']
        atdb[atsym] = calc_atomic_weight(v)
    return atdb
        
def get_comp(mf):
    return re.findall(r'([A-Z][a-z]*)(\d*)',mf)

def bounded_mw(mf,saw):
    bmw = np.array([0.,0.])
    atcomp = get_comp(mf)
    for at in atcomp:
        rng = np.array(saw[at[0]])
        num = int(at[1] if at[1] else 1)
        bmw += num * rng
        print (num, rng, bmw)
    return bmw

def extract_values(value):
    ''' extract significant figures from value in the form x.xxxx(uu) '''
    # get numeric values
    print('val',value)
    if '(' in value:
        val, u = re.findall(r'(\d+[.,]*?\d*)\((\d+)\)',value)[0]
        # count significant figures
        sigfig = len(val.split('.')[1])
        # multiply uncertainty value with significant figures
        u = float(u) * 10**(-sigfig)
        val = float(val)
    else:
        val = float(value)
        u = 0
    return val, u

def calc_atomic_weight(atdic):
    ''' atomic weight is calculated as the weighted average of relative
        atomic mass for all isotopes '''
    
    atweight = 0
    totunc = 0
    for iso in atdic:
        if iso['Isotopic Composition'] == '' :
            continue
        # extract values from Isotopic Composition as pair (value, uncertainty)
        isocomp = extract_values(iso['Isotopic Composition'])
        relmass = extract_values(iso['Relative Atomic Mass'])
        # relative contribution to atomic weight
        atw = isocomp[0] * relmass[0]
        # error propagation
        unc = (isocomp[1]/isocomp[0])**2 + (relmass[1]/relmass[0])**2
        print(isocomp,relmass,atw,unc)
        atweight += atw
        totunc += unc * atw**2
    totunc = np.sqrt(totunc)
    return atweight, totunc
        
def calc_mol_weight(mf, atdb):
    ''' calulate molecular weight and corresponding uncertainty'''
    
    comp = get_comp(mf)
    mw = 0
    unc = 0
    for c in comp:
        print(c)
        sym = c[0]
        num = int(c[1]) if c[1] else 1
        mw += num * atdb[sym][0]
        unc += (num * atdb[sym][1])**2
    return mw, np.sqrt(unc)


if __name__ == "__main__":
    db, saw = load_db()
    atdb = make_atomic_weight_db(db)
    print("The molecular weight of water is: %f +/- %f"%calc_mol_weight('H2O',atdb))
    