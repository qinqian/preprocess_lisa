import os
import configparser
import argparse

p = argparse.ArgumentParser()
p.add_argument('--folder', help='input folder')

args = p.parse_args()

conf = configparser.ConfigParser()
conf.read('lisa.conf')

for key in conf.sections():
    for i in conf[key].keys():
        conf.set(key, i, os.path.join(args.folder, os.path.basename(conf.get(key, i))))

with open('lisa.cfg', 'w') as configfile:
    conf.write(configfile)

