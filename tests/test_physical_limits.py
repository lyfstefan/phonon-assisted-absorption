import unittest
import numpy as np
from phonon_absorption import Channel,model_i,model_ii,sideband_weights


class Limits(unittest.TestCase):
    def test_uncoupled_lorentzian(self):
        x=np.linspace(.9,1.4,137)
        c=Channel(1.2,0,.7)
        exact=.7*.02/((x-1.2)**2+.02**2)
        for fn in [model_i,model_ii]:
            np.testing.assert_allclose(fn(x,c,.03,.02,.02),exact,rtol=1e-14)

    def test_zero_temperature_weights_and_mean(self):
        m,w=sideband_weights(2.3,.032,0)
        self.assertEqual(w[m<0].sum(),0)
        self.assertAlmostEqual(float(w.sum()),1,places=13)
        self.assertAlmostEqual(float((m*w).sum()),2.3,places=12)

    def test_parameters(self):
        c=Channel(1.2,1,1)
        for omega,line,temp in [(0,.02,.02),(.03,0,.02),(.03,.02,-1)]:
            with self.assertRaises(ValueError):model_i([1.1],c,omega,line,temp)
        with self.assertRaises(ValueError):Channel(1.2,-1,1)


if __name__=="__main__":unittest.main()
