import unittest
from PyQuante.HFSolver import HFSolver
import logging
class UHFSolver(HFSolver):
    def setup(self,**opts):
        self.nalpha,self.nbeta = self.molecule.get_alphabeta()
        HFSolver.setup(self,**opts)
        self.method = "UHF"
        return

    def print_setup_info(self):
        HFSolver.print_setup_info(self)
        logging.debug("Nalpha = %d" % self.nalpha)
        logging.debug("Nbeta = %d" % self.nbeta)
        return

    def setup_guess(self,**opts):
        HFSolver.setup_guess(self,**opts)
        self.orbea = self.orbeb = self.orbe
        self.orbsa = self.orbsb = self.orbs
        return

    def setup_averaging(self,**opts):
        # Taken out because DIIS just doesn't work for UHF
        #  without some additional hacking
        #from PyQuante.Convergence import DIIS
        #self.avga = DIIS(self.S)
        #self.avgb = DIIS(self.S)
        self.avgfact = opts.get('avgfact',0.5)
        return        

    def update_density(self):
        from PyQuante.LA2 import mkdens
        from PyQuante.fermi_dirac import mkdens_fermi
        if self.etemp:
            self.Da,self.entropya = mkdens_fermi(2*self.nalpha,self.orbea,
                                                 self.orbsa,self.etemp)
            self.Db,self.entropyb = mkdens_fermi(2*self.nbeta,self.orbeb,
                                                 self.orbsb,self.etemp)
            self.entropy = (self.entropya+self.entropyb)/2.0
        else:
            self.Da = mkdens(self.orbsa,0,self.nalpha)
            self.Db = mkdens(self.orbsb,0,self.nbeta)
            self.entropy=0
        if self.do_averaging:
            if self.iter > 1:
                self.Da = self.avgfact*self.Da + (1.0-self.avgfact)*self.Da_old
                self.Db = self.avgfact*self.Db + (1.0-self.avgfact)*self.Db_old
            self.Da_old = self.Da
            self.Db_old = self.Db
        self.Dab = self.Da + self.Db
        return

    def update_J(self):
        from PyQuante.Ints import getJ
        from PyQuante.LA2 import trace2
        self.Ja = getJ(self.Ints,self.Da)
        self.Jb = getJ(self.Ints,self.Db)
        return

    def update_K(self):
        from PyQuante.Ints import getK
        from PyQuante.LA2 import trace2
        self.Ka = getK(self.Ints,self.Da)
        self.Kb = getK(self.Ints,self.Db)
        return

    def update_fock(self):
        self.update_J()
        self.update_K()
        self.Fa = self.h + self.Ja + self.Jb - self.Ka
        self.Fb = self.h + self.Ja + self.Jb - self.Kb
        # DIIS doesn't work well for UHF:
        #if self.do_averaging:
        ##    self.Fa = self.avga.getF(self.Fa,self.Da)
        ##    self.Fb = self.avgb.getF(self.Fb,self.Db)
        return

    def solve_fock(self):
        from PyQuante.LA2 import geigh
        self.orbea,self.orbsa = geigh(self.Fa,self.S)
        self.orbeb,self.orbsb = geigh(self.Fb,self.S)
        return

    def calculate_energy(self):
        from PyQuante.LA2 import trace2
        self.Eone = trace2(self.Dab,self.h)
        self.Ej = 0.5*trace2(self.Dab,self.Ja+self.Jb)
        self.Exc = -0.5*(trace2(self.Da,self.Ka)+trace2(self.Db,self.Kb))
        self.energy = self.Eone + self.Ej + self.Exc + self.Enuke + self.entropy
        return        

class UnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)

    def testLiUHF(self):
        li_uhf = UHFSolver(self.li)
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431364,4)

    def testLiUHFFT(self):
        li_uhf = UHFSolver(self.li)
        li_uhf.iterate(etemp=1e4)
        return

def test():
    suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
if __name__ == '__main__': test()
