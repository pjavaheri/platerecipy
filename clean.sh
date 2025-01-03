#!/bin/bash
rm -rf ./src/*.egg-info
cd doc && latexmk -c && rm manual.bbl manual.run.xml