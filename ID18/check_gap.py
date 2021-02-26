#
# Dear Marco,
#
# We expect a K = 2.2 at 4.8 mm. This is more or less what we obtained with the ID15 CPMU18, at lower gap.
#
# You should use 130 pm rad horizontal and 10 pm rad vertical emittance.
#
# Regards
#
# Gaël Le Bec

import numpy
import scipy.constants as codata
# following /segfs/tango/jsrund  A(0)=5  A(1)=1  A(>1)=0  for ID15 U18
A0 = 5.0 # 6.0 # 5.0
A1 = 1.0 # 0.98 # 1.0

# what is the K for gap=4.8?

period = 18.0
period_in_m = period / 1000

for gap in [6.0, 4.8]:
    #after email Juan 2020-11-12
    B =  A1 * numpy.exp(- 1 * numpy.pi * (gap - A0) / period )

    K = B /(2 * numpy.pi * codata.m_e * codata.c / (codata.e * period_in_m))
    print("gap: %g mm, B: %g T, K: %g" % (gap, B, K) )
