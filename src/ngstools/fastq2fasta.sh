#!/bin/bash
cat $1 | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $2
