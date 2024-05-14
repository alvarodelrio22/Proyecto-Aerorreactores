import numpy as np, csv
from Components.ComponentMap import compressor, turbine
from Components.DesignVariables import m_LPC_design

compressor(1,1.08,"beta","N","pi","LPC")
compressor(m_LPC_design,0.9949538584005653,"m","N","pi","LPC")

