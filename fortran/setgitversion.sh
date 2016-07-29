#!/bin/bash

# Get git version
if git --version 2>&1 1> /dev/null
then
    vers=`git log --pretty="format:%h" -n1 HEAD`
    versstring=`printf "%17s" $vers`
else
    versstring=`printf "%17s" "N\/A"`
fi

sed -e "s/dummy/${versstring}/;" gitversion_template.f90 > gitversion.f90

