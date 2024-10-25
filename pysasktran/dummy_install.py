# Runs some sasktranif things to create registry entries

print("Starting dummy install")

import sasktran

print("Sasktran imported")
import sasktranif.sasktranif as skif

opt_prop = skif.ISKOpticalProperty('MIEAEROSOL_H2SO4')