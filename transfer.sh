#!/bin/bash

for i in `cut -f 1 ../marge2/data.txt`; do
  rsync -auviz --progress '-e ssh -p 33001' qqin@cistrome.org:${i} .
done

