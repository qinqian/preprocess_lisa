#!/usr/bin/env python

import os
import configparser

c = configparser.ConfigParser()

c.read('lisa.conf')

for s in c.sections():
    for k in c[s].keys():
        if k == 'bwa_index':
            continue
        if os.path.exists(c.get(s, k)):
            print(c.get(s, k))
        else:
            print(c.get(s, k))
            raise Exception
