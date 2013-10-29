#!/bin/bash

cd $(dirname ${0})/..
# remove all config.h includes
find . -name "*.hh" | xargs sed -i "s;#include\ \<config\.h\>;;g"
find . -name "*.cc" | xargs sed -i "s;#include\ \<config\.h\>;;g"
# insert as 1st line
find . -name "*.cc" | xargs sed -i "1i #include\ <config\.h>"

