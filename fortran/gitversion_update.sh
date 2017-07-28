#!/bin/bash

# Get git version
if git --version 2>&1 1> /dev/null
then
    vers=`git log --pretty="format:%h" -n1 HEAD`
    versionstring=`printf "%17s" $vers`
else
    versionstring=`printf "%17s" "N\/A"`
fi
sed "s/tmp/${versionstring}/" gitversion_template > gitversion
