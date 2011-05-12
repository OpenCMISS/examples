#!/bin/bash

#rename the general field to source so that cmgui doesn't complain

sed -i -e "s%4) general%4) source%" Poisson.part0.exnode
sed -i -e "s%5) general%5) source%" Poisson.part0.exelem
