#!/bin/bash
cp setInput.m ./gainCurve/
cd gainCurve
matlab -nosplash -nodisplay -r "setInput;exit"
