#!/bin/bash



awk 'BEGIN{IF=" "}{if (int($1)!=1604 && int($1)!=1061) {print $2}}' nodnase_10sample.du

