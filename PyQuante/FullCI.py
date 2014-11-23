"""
Full CI based on determinants rather than on CSFs.

The approach is the one introduced by Olsen J Chem Phys 89 2185 (1988)

It is also described in the book Molecular electronic structure
theory, by Helgaker, Olsen, Jorgensen.  There it is called 'Minimal
operator count (MOC) method'

written by Simon Sala

Notation:
Book of Helgaker.

"""

import itertools
import numpy as np
import scipy.sparse as spspa
import scipy.sparse.linalg as spspalin
import scipy.linalg as splin
from scipy.special import binom


try: from PyQuante.cints import ijkl2intindex
except: 
    print "cints import failed in CI.py"	
    from PyQuante.pyints import ijkl2intindex

from PyQuante.CI import TransformInts


def single_excitations(n):
    singles = []
    for p in xrange(n):
        for q in xrange(n):
            singles.append((p,q))
    return singles

def double_excitations(n):
    doubles = []
    for p in xrange(n):
        for q in xrange(n):
            for r in xrange(n):
                for s in xrange(n):
                    doubles.append((p,q,r,s))
    return doubles

def transform_one_ints(h,orbs):
    """ Transform the one-electron Hamilton matrix from basis function
    representation to MO basis, 

    orbs is the coefficient matrix with rows indexing orbitals and
    colums indexing the basis function coefficients.

    See 
    http://vergil.chemistry.gatech.edu/resources/programming/mp2-transform-project.pdf
    for details.

    For very large basis sizes, this might need to be calculated on
    the fly.
    """
    return np.dot(orbs.T, np.dot(h,orbs))

def e_pq_on_string(p,q,string):
    """ 
    apply the excitation operator a^+_p a_q on a string 
    
    This gives new string and a phase factor.

    It must have been checked that q is in string and p is not!
    """

    if q not in string:
        """ annihilate vacuum """
        return 0,0
    if p in string and p!=q:
        """ try to create already occupied orbital which was
        not destroyed """
        return 0,0

    # action of E_pq on string j gives new string e_pq_string:
    e_pq_string = list(string)
    # determine phase factor
    phase_q = (-1)**e_pq_string.index(q)
    # apply annihilator q
    e_pq_string.remove(q)
    # apply creator p
    e_pq_string.append(p)
    e_pq_string.sort()
    phase_p = (-1)**e_pq_string.index(p)
    
    return phase_p*phase_q, e_pq_string


class FCISolver(object):
    """ Interface to the scipy.sparse.linalg.eigs eigenvalue solver"""

    def __init__(self, h, ERI, enuke, orbs, n_elec, multiplicity, m_s, k=4, sigma_eigs=None, which='SA', v0=None, maxiter=None, tol=0, return_eigenvectors=True ):
        """
        Parameters: 
        
        h           :       one-electron integrals over basis functions
        
        ERI         :       electron repulsion integrals over basis functions

        enuke       :       The nuclear attraction energy, Molecule.get_enuke()
        
        orbs        :       coefficient matrix from HF calculation giving the 
                            orbitals in rows and the bfs coeffs in columns
       
        n_elec      :       total number of electron

        multiplicity:       2*S+1
        
        m_s         :       M_s component of total spin.

        
        keyword parameters passed to eigs solver (see scipy docs):
        
        k           :       number of eigenvalues computed
        
        sigma_eigs  :       number to which eigenvalues are close 
                            (should be set to increase performance)
        
        which       :       if set to 'SR' calculate k smalles real part egenvalues
        
        v0          :       initial vector to start from
                            perhaps HF vector (1,0,0...)
                            
        maxiter     :       maximum number of Arnoldi updates allowed
        
        tol         :       tolerance in calculation of eigenvalues
                            0 means machine precision
                            
        return_eigenvectors: return eigenvector in addition to eigenvalues
                             if set to True
        
        """

        self.enuke = enuke
        self.k = k
        self.sigma_eigs = sigma_eigs
        self.which = which
        self.v0 = None
        self.maxiter=maxiter
        self.tol=tol
        self.return_eigenvectors = return_eigenvectors

        # number of alpha electrons
        self.n_alpha = 0.5*n_elec + m_s
        # number of beta electrons
        self.n_beta = 0.5*n_elec - m_s
        # number of orbitals
        self.n_orbs = orbs.shape[0]
        # number of alpha strings
        self.len_alpha = int(binom(self.n_orbs,self.n_alpha))
        # number of beta strings
        self.len_beta = int(binom(self.n_orbs,self.n_beta))

        assert self.n_alpha +self.n_beta == n_elec

        # Instantiate Sigma class
        self.SigmaInst = Sigma(np.eye(self.len_alpha, self.len_beta), h, ERI, orbs, n_elec, multiplicity, m_s)
        
        # shape of the H matrix
        self.H_mat_shape = (self.len_alpha*self.len_beta , self.len_alpha*self.len_beta)
        # shape of the coefficient matrix in Sigma class
        self.c_mat_shape = (self.len_alpha , self.len_beta)
        # shape of the corresponding vector passed to eigs
        self.c_vec_shape = self.len_alpha*self.len_beta

        # Linear operator passed to eigensolver
        self.LinOp = spspalin.LinearOperator(self.H_mat_shape, self.matvec, dtype=np.float64)
        

    def matvec(self, vec):
        """ The reshaped matrix vector step needed for the iterations
        in eigs solver. 

        The steps are:
        1. reshape vec to matrix
        2. get sigma
        3. reshape back and return
        """
        
        vec_mat = vec.reshape( self.c_mat_shape )
        self.SigmaInst.c_mat = vec_mat
        new_vec_mat = self.SigmaInst.get_sigma()
        return new_vec_mat.reshape(self.c_vec_shape)
        
        
    def iterate(self):
        eva, eve = spspalin.eigsh(self.LinOp,k=self.k, sigma = self.sigma_eigs, which = self.which, v0 = self.v0, maxiter= self.maxiter, tol=self.tol, return_eigenvectors = self.return_eigenvectors)
        print "diagonalization sucessful"

        self.eva, self.eve = self.sort_and_add_enuke(eva,eve)
        
        return self.eva, self.eve
        

    def sort_and_add_enuke(self, eva, eve):
        """ sort the eva end eve and add the nuclear attraction energy
        to eva. """

        # sort
        indx = eva.argsort()
        eva = eva[indx] 
        eve = eve[:,indx]

        # add enuke
        eva += self.enuke

        return eva, eve
        
        
        

class FCIExactSolver(object):
    """
    In contrast to FCISolver, this method build the full CI
    Hamiltonian matrix explicitly, and then diagonalizes it exactly.
    It is only suitable for small CI spaces and is more intendend for
    debugging purposes.
    """

    def __init__(self, h, ERI, enuke, orbs, n_elec, multiplicity, m_s):
        """
        Parameters: 
        
        h           :       one-electron integrals over basis functions
        
        ERI         :       electron repulsion integrals over basis functions
        
        orbs        :       coefficient matrix from HF calculation giving the 
                            orbitals in rows and the bfs coeffs in columns
       
        n_elec      :       total number of electron

        multiplicity:       2*S+1
        
        m_s         :       M_s component of total spin.

        """
        
        # Instantiate FCISolver class to access necessarry structures.
        self.FCISolverInst = FCISolver(h, ERI, enuke, orbs, n_elec, multiplicity, m_s)

        
    def get_H_mat(self):
        """ build the Hamiltonian matrix in the I_c = I_alpha I_beta space.
        The principle is as follows:
        
        With the Sigma class we have a (hopefully efficient) method to
        calculate priducts of the Hamiltonian matrix (tensor) with a
        coefficient vector (matrix). 
        
        The, e.g., 1st column of a matrix A is obtained by the
        multiplication of A with the vector (1,0,0...,0).

        This principle is applied for each of the len_alpha*len_beta
        components of the coefficient vector.

        The reshaping of the coeffitient vector to an coefficient
        matrix is handled by the matvec method of FCISolver class.
        """
        
        self.H_mat = np.zeros((self.FCISolverInst.len_alpha*self.FCISolverInst.len_beta,self.FCISolverInst.len_alpha*self.FCISolverInst.len_beta))
        
        for col in xrange(self.FCISolverInst.len_alpha*self.FCISolverInst.len_beta):
            """ loop over c_mat vector """
            vec = np.zeros((self.FCISolverInst.len_alpha*self.FCISolverInst.len_beta))
            vec[col] = 1.
            self.H_mat[:,col] = self.FCISolverInst.matvec(vec)
        
        print "build of H_mat successful."

        
    def diagonalize(self):
        """ diagonalize the Hamiltonian matrix """
        try: self.H_mat
        except: self.get_H_mat()
        
        eva, eve = splin.eigh(self.H_mat)
        
        self.eva, self.eve = self.FCISolverInst.sort_and_add_enuke(eva,eve)

        print "diagonalization successful"""
        return self.eva, self.eve
    


class Graph(object):
    """ graph object determining vertex weights and arc weights in
    reverse lexical ordering. 

    see Helgaker section 11.8.2

    Different to Helgaker:
    Attention: orbital numbering starts with 0!
    Attention: address starts also with 0!

    """
    
    def __init__(self, n_orbs, n_electron):
        """
        n_orbs      :        number of orbitals
        n_electron  :        number of electrons
        """
        self.n_orbs = n_orbs
        self.n_electron = n_electron
        self.get_vertex_weights()
        self.get_arc_weights()
        self.get_occupations()

        assert int(binom(self.n_orbs,self.n_electron)) == len(self.occupations)
        
    def get_vertex_weights(self):
        """ get the vertex weights

        vertices are indexed as a two-dimensional n_orbs+1 x
        n_electron+1 array:
        
        rows: orbitals
        columns: number of electrons
        """
        self.vert_weights = np.zeros((self.n_orbs+1,self.n_electron+1), dtype=np.int32)
        self.vert_weights[0,0] = 1
        for row in xrange(1,self.n_orbs+1):
            for column in xrange(self.n_electron+1):
                
                if column > row:
                    """ upper triangle is left out """
                    continue
                if row > column+ self.n_orbs - self.n_electron:
                    continue
                
                if column==0:
                    """check if vertex is allowed"""
                    self.vert_weights[row,column] = self.vert_weights[row-1,column] 
                else:
                    self.vert_weights[row,column] = self.vert_weights[row-1,column] + self.vert_weights[row-1,column-1]

            

    def get_arc_weights(self):
        """ get the arc weights 
        
        arc weigths for vertical arcs. Represented as (n,N) array
        
        """
        self.arc_weights = np.zeros((self.n_orbs, self.n_electron), dtype=np.int32)
        for row in xrange(self.n_orbs):
            for column in xrange(self.n_electron):
                if column > row:
                    """ upper triangle is left out """
                    continue
                if row > column+ self.n_orbs - self.n_electron:
                    """ lower part """
                    continue
                
                self.arc_weights[row,column] = self.vert_weights[row,column+1]
       
        
    def address(self, occupation):
        """ get the address of a string given its occupation as, e.g.,
        (0,2,3) means string a^+_0 a^+_2 a^+_3
        
        Attention: orbital numbering starts with 0!
        Attention: address starts also with 0!

        occupation      :     SORTED list of creation operators (integers)
        """
        address = 0
        for index in xrange(self.n_electron):
            address += self.arc_weights[occupation[index],index]
        return address
            
            
    def get_occupations(self):
        """ return a list of occupations (list of lists) in reverse
        lexical order 

        Strategy:
        create all occupations and the sort by address.
        """
        occs = list(itertools.combinations(range(self.n_orbs), self.n_electron))
        occs = sorted(occs, key=lambda occ:self.address(occ))
        self.occupations = occs

        

class Sigma(object):
    """ class to compute sigma matrix with a given index matrix c_mat """
    
    def __init__(self, c_mat, h, ERI, orbs, n_elec, multiplicity, m_s):
        """
        
        c_mat       :       The coefficient matrix indexed with alpha and
                            beta string addressing

        h           :       one-electron integrals over basis functions
        
        ERI         :       electron repulsion integrals over basis functions

        
        orbs        :       coefficient matrix from HF calculation giving the 
                            orbitals in rows and the bfs coeffs in columns
       
        n_elec      :       total number of electron

        multiplicity:       2*S+1
        
        m_s         :       M_s component of total spin.
        
        """

        self.c_mat = c_mat
        self.h = h
        self.ERI = ERI
        self.orbs = orbs
        self.n_elec = n_elec
        self.multiplicity = multiplicity
        self.m_s = m_s
        

        # see (11.6.1)
        self.n_alpha = int(0.5*self.n_elec + self.m_s)
        self.n_beta = int(0.5*self.n_elec - self.m_s)
        assert self.n_alpha+self.n_beta == self.n_elec

        self.n_orbs = self.orbs.shape[0]
        self.singles = single_excitations(self.n_orbs)
        self.doubles = double_excitations(self.n_orbs)
        
        self.AlphaStrings = Graph(self.n_orbs, self.n_alpha)
        self.BetaStrings = Graph(self.n_orbs, self.n_beta)

        
    def get_sigma(self):
        """ add all sigma components togetther:

        one electron terms: sigma_alpha, sigma_beta (11.8.20)

        two electron terms: sigma_alpha_alpha, sigma_beta_beta,
        sigma_alpha_beta (11.8.29)
        
        """
        
        # if self.m_s == 0:
        if False:
            """ In case of 0 spin projection, i.e. as many spin up
            (alpha) as spin down(beta) electrons, sigma_beta can be
            obtained from sigma_alpha and sigma_alpha_beta can be
            obtained from sigma_alpha_alpha. See 
            
            http://vergil.chemistry.gatech.edu/notes/ci/ci.pdf
            
            and 

            J Chem Phys 89 2185 (1988)
            
            for details.

            ATTENTION. Disable this option for the moment. It gives me
            incorrect results. Gives correct results if second line is
            
            sigma_beta = ...* sigma_alpha
            
            i.e. no transpose of sigma_alpha.

            As long as I don't understand this, it should be disabled!

            """
            sigma_alpha = self.get_sigma_alpha() + self.get_sigma_alpha_alpha()
            # (-1)**S * sigma_alpha.T (Eq.15 in J Chem Phys paper)

            sigma_beta = (-1)**(0.5*(self.multiplicity - 1)) * sigma_alpha.T
            
            return sigma_alpha + sigma_beta + self.get_sigma_alpha_beta()
            
        else:
            """ M_s != 0 """
            return self.get_sigma_alpha() + self.get_sigma_beta() + self.get_sigma_alpha_alpha() + self.get_sigma_beta_beta() + self.get_sigma_alpha_beta()


    # one electron part


    def get_sigma_alpha(self):
        """ (11.8.26) """
        try: self.k_alpha
        except: self.get_k_alpha()
        # dot product of sparse matrix k_alpha and dense matrix c_mat
        return self.k_alpha*self.c_mat
    
    def get_sigma_beta(self):
        """ (11.8.27), beware the transpose!!!"""
        try: self.k_beta
        except: self.get_k_beta()
        # dot product of sparse matrix k_beta and dense matrix c_mat
        return self.c_mat*self.k_beta.transpose()

    def get_k_alpha(self):
        """ (11.8.24) """
        self.k_alpha = self.get_k_gamma("alpha")
        print " build of k_alpha successful"

    def get_k_beta(self):
        """ (11.8.25) """
        self.k_beta = self.get_k_gamma("beta")
        print " build of k_beta successful"

    def get_k_gamma(self, alpha_beta):
        """ get k_sigma matrices (11.8.24) 
        gamma is alpha or beta.

        It is sparse! Hence is it constructed solely from the non-zero
        elements. 
        """
        try: self.k_mat
        except: self.get_k_mat()
        
        if alpha_beta == "alpha":
            Strings = self.AlphaStrings
        elif alpha_beta == "beta":
            Strings = self.BetaStrings
        else:
            raise ValueError, 'argument alpha_beta must be alpha or beta'

        row_index = []
        column_index = []
        data = []
        
        for i_string in Strings.occupations:
            for j_string in Strings.occupations:
                
                row = Strings.address(i_string)
                column = Strings.address(j_string)
                elem = 0

                if row == column:
                    """ strings are equal """
                    elem = 0
                    for occ in i_string:
                        elem += self.k_mat[occ,occ]

                    row_index.append(row)
                    column_index.append(column)
                    data.append(elem)
                    continue
                
                for p,q in self.singles:
                    """ loop over excitations """
    
                    # apply excitation operator on string:
                    phase, e_pq_j = e_pq_on_string(p,q,j_string)
                    if e_pq_j == 0:
                        """ tried to annihilate vaccum or to create
                        doubly """
                        continue
                
                    if row != Strings.address(e_pq_j):
                        """ strings differed by more than the pair p q """
                        continue
                    else:
                        elem += phase*self.k_mat[p,q]
                        # row_index.append(row)
                        # column_index.append(column)
                        # data.append(elem)
                        # # there will not be another pq that satisfies,
                        # # exit pq loop
                        # break

                if abs(elem) > 1e-14:
                    row_index.append(row)
                    column_index.append(column)
                    data.append(elem)

                    
        return spspa.csr_matrix( (np.array(data),(np.array(row_index),np.array(column_index))), shape=(len(Strings.occupations),len(Strings.occupations)) )


    def get_k_mat(self):
        """ build k_pq from (11.8.8) """
        try: self.h_mat
        except: self.h_mat = transform_one_ints(self.h,self.orbs)
        try: self.MOInts
        except: self.MOInts = TransformInts(self.ERI,self.orbs)
        # except: self.MOInts = TransformInts_test(self.ERI,self.orbs)

        self.k_mat = np.zeros((self.n_orbs,self.n_orbs))
        for p,q in self.singles:
            for r in xrange(self.n_orbs):
                self.k_mat[p,q] -= 0.5*self.MOInts[ijkl2intindex(p,r,r,q)]
                # self.k_mat[p,q] -= 0.5*self.MOInts[p,r,r,q]
            self.k_mat[p,q] += self.h_mat[p,q]

        print "build of k_mat successful"    
        

    # two-electron part

    
    def get_sigma_alpha_alpha(self):
        """ (11.8.35) """
        try: self.G_alpha
        except: self.get_G_alpha()
        # dot product of sparse matrix k_alpha and dense matrix c_mat
        return self.G_alpha*self.c_mat
    
    def get_sigma_beta_beta(self):
        """ (11.8.36), beware the transpose!!!"""
        try: self.G_beta
        except: self.get_G_beta()
        # dot product of sparse matrix k_beta and dense matrix c_mat
        return self.c_mat* self.G_beta.transpose()

    def get_G_alpha(self):
        """ (11.8.33) """
        self.G_alpha = self.get_G_gamma("alpha")
        print "build of G_alpha successful"

    def get_G_beta(self):
        """ (11.8.34) """
        self.G_beta = self.get_G_gamma("beta")
        print "build of G_beta successful"

    def get_G_gamma(self, alpha_beta):
        """ get G_sigma matrices (11.8.33/34) 
        gamma is alpha or beta.

        It is sparse! Hence is it constructed solely from the non-zero
        elements. 
        """

        try: self.MOInts
        except: self.MOInts = TransformInts(self.ERI,self.orbs)
        # except: self.MOInts = TransformInts_test(self.ERI,self.orbs)

            
        if alpha_beta == "alpha":
            Strings = self.AlphaStrings
        elif alpha_beta == "beta":
            Strings = self.BetaStrings
        else:
            raise ValueError, 'argument alpha_beta must be alpha or beta'

        row_index = []
        column_index = []
        data = []
        
        for i_string in Strings.occupations:
            for j_string in Strings.occupations:
                
                row = Strings.address(i_string)
                column = Strings.address(j_string)
                
                elem = 0

                for p,q,r,s in self.doubles:
                    """ loop over excitations """
                    
                    # apply excitation E_rs operator on string:
                    phase_rs, e_rs_j = e_pq_on_string(r,s,j_string)
                    if e_rs_j == 0:
                        """ tried to annihilate vaccum or to create
                        doubly """
                        continue

                    # apply excitation E_pq operator on string:
                    phase_pq, e_pqrs_j = e_pq_on_string(p,q,e_rs_j)
                    if e_pqrs_j == 0:
                        """ tried to annihilate vaccum or to create
                        doubly """
                        continue
                    
                    if row != Strings.address(e_pqrs_j):
                        """ strings differed by more than two pairs p q r s """
                        continue
                    else:
                        elem += 0.5 *phase_pq *phase_rs *self.MOInts[ijkl2intindex(p,q,r,s)]
                        # elem += 0.5 *phase_pq *phase_rs *self.MOInts[p,q,r,s]

                        ### Need to think when can exit the loop. For sure if p!=q!=r!=s
                        # if p!=q and q!=r and r!=s:
                        #     """ there will not be another pqrs that
                        #     satisfies, exit pqrs loop """
                        #     row_index.append(row)
                        #     column_index.append(column)
                        #     data.append(elem)
                        #     break

                if abs(elem) > 1e-14:
                    row_index.append(row)
                    column_index.append(column)
                    data.append(elem)

        return spspa.csr_matrix( (np.array(data),(np.array(row_index),np.array(column_index))), shape=(len(Strings.occupations),len(Strings.occupations)) )
    


    def get_sigma_alpha_beta(self):
        """ (11.8.39) """
        
        sigma_alpha_beta = np.zeros( (len(self.AlphaStrings.occupations),len(self.BetaStrings.occupations)) )
        
        for p,q in self.singles:
            """ Matrix summation """
            sigma_alpha_beta += self.get_sigma_alpha_beta_pq(p,q)
        
        return sigma_alpha_beta
        
    
    def get_sigma_alpha_beta_pq(self,p,q):
        """ (11.8.43) 

        Dot product of dense matrix with sparse matrix.
        
        temp is the matrix multiplication of <I_a|E_pq|J_a>* C_Ia,Ja
        from Eq. (11.8.41)

        """
        try: self.D_alpha_ia_jb_pq_list
        except: self.get_D_alpha_ia_jb_pq_list()

        try: self.G_beta_ib_jb_pq_list
        except : self.get_G_beta_ib_jb_pq_list()

        temp = self.D_alpha_ia_jb_pq_list[p][q] * self.c_mat
        return temp* self.G_beta_ib_jb_pq_list[p][q].transpose()



    def get_D_alpha_ia_jb_pq_list(self):
        """ 
        create a list of lists containing the sparse matrices
        D_alpha_ia_jb_pq(p,q) from (11.8.41) without the product with
        C_JaJb
        """
        self.D_alpha_ia_jb_pq_list = []
        
        for p in xrange(self.n_orbs):
            """ loop manually over excitations """
            row=[]
            for q in xrange(self.n_orbs):
                row.append(self.get_D_alpha_ia_jb_pq(p,q))
            self.D_alpha_ia_jb_pq_list.append(row)
        
        print "build of D_alpha_ia_jb_pq_list successful"
        

    def get_D_alpha_ia_jb_pq(self,p,q):
        """ (11.8.41) but matrix product with c_mat is pulled into
        (11.8.39) """

        row_index=[]
        column_index=[]
        data=[]

        for ia_string in self.AlphaStrings.occupations:
            """
        
            set up matrix <I_alpha|e^alpha_pq|J_alpha>


            apply <I_alpha| E_pq^alpha = <I_alpha| a^+_p a_q 
            = a^+_q a_p |I_alpha>.

            If this vanishes, the row I_alpha is zero, since
            D_IJ^alpha [pq] is zero in (11.8.41).

            This is done to get minimal operator count.
            """
            phase, e_qp_ia = e_pq_on_string(q,p,ia_string)
            if e_qp_ia == 0:
                continue

            for ja_string in self.AlphaStrings.occupations:
                
                row = self.AlphaStrings.address(ia_string)
                column = self.BetaStrings.address(ja_string)
                
                if column == self.AlphaStrings.address(e_qp_ia):
                    """ nonzero element in <I_alpha|e^alpha_pq|J_alpha>"""
                    row_index.append(row)
                    column_index.append(column)
                    data.append(phase)
                
        return spspa.csr_matrix( (np.array(data),(np.array(row_index),np.array(column_index))), shape=(len(self.AlphaStrings.occupations),len(self.AlphaStrings.occupations)) )

        
    def get_G_beta_ib_jb_pq_list(self):
        """ 
        create a list of lists containing the sparse matrices 
        G_beta_ib_jb_pq(p,q) from (11.8.42)
        """
        self.G_beta_ib_jb_pq_list = []
        
        for p in xrange(self.n_orbs):
            """ loop manually over excitations """
            row=[]
            for q in xrange(self.n_orbs):
                row.append(self.get_G_beta_ib_jb_pq(p,q))
            self.G_beta_ib_jb_pq_list.append(row)

        print "build of G_beta_ib_jb_pq_list successful"
            

    def get_G_beta_ib_jb_pq(self,p,q):
        """ (11.8.42) """

        row_index=[]
        column_index=[]
        data=[]

        
        for ib_string in self.BetaStrings.occupations:
            """
            """
            for jb_string in self.BetaStrings.occupations:
                
                row = self.BetaStrings.address(ib_string)
                column = self.BetaStrings.address(jb_string)
                elem = 0

                for r,s in self.singles:
                    """ loop over excitations """
    
                    # apply excitation operator on string:
                    phase, e_rs_jb = e_pq_on_string(r,s,jb_string)
                    if e_rs_jb == 0:
                        """ tried to annihilate vaccum or to create
                        doubly """
                        continue
                
                    if row != self.BetaStrings.address(e_rs_jb):
                        """ strings differed by more than the pair p q """
                        continue
                    else:
                        elem += phase*self.MOInts[ijkl2intindex(p,q,r,s)]
                        # elem += phase*self.MOInts[p,q,r,s]

                if abs(elem) > 1e-14:
                    row_index.append(row)
                    column_index.append(column)
                    data.append(elem)
                
        return spspa.csr_matrix( (np.array(data),(np.array(row_index),np.array(column_index))), shape=(len(self.BetaStrings.occupations),len(self.BetaStrings.occupations)) )
                        
        


def test_graph():
    Inst = Graph(5,3)
    print "address 0,1,2 = %i"%Inst.address([0,1,2])
    print "address 0,1,3 = %i"%Inst.address([0,1,3])
    print "address 1,2,3 = %i"%Inst.address([1,2,3])    
    print "address 0,3,4 = %i"%Inst.address([0,3,4])    
    print "address 2,3,4 = %i"%Inst.address([2,3,4])
    print ""
    print "All occupations in reverse lexical ordering:"
    print Inst.occupations


def test_e_pq_on_string():
    print "E_11 on [0,2] = ", e_pq_on_string(1,1,[0,2])
    print "E_22 on [0,2] = ", e_pq_on_string(2,2,[0,2])
    print "E_12 on [0,2] = ", e_pq_on_string(1,2,[0,2])
    print "E_21 on [0,2] = ", e_pq_on_string(2,1,[0,2])
    print "E_02 on [1,2] = ", e_pq_on_string(0,2,[1,2])


def test_fci():
    """ test FCI calculation"""
    from Molecule import Molecule
    from PyQuante import SCF


    nel = 2
    mult = 1
    m_s = 0
    k=10
    
    # h2 = Molecule('h2',[(1,(0,0,0)),(1,(1.4,0,0))], multiplicity = mult)
    h2 = Molecule('h2',[(1,(0,0,0)),(1,(1.4,0,0))])

    Solver = SCF(h2,method = "HF")
    Solver.iterate()
    print "orbital energies ",Solver.solver.orbe
    print "HF energy = ",Solver.energy
    
    # FCIInst = FCISolver(Solver.h, Solver.ERI, Solver.solver.orbs, nel, mult, m_s, k=k, sigma_eigs=None)
    # eva, eve = FCIInst.iterate()    

    FCIInst = FCIExactSolver(Solver.h, Solver.ERI, h2.get_enuke(), Solver.solver.orbs, nel, mult, m_s)
    eva,eve = FCIInst.diagonalize()

    print "eva = ", eva
    # print "eve = ",eve
    print "correlation energy = ", eva[0] - Solver.energy
    print "correlation energy should be (with 6-31g**) -0.03387 a.u."


if __name__ == "__main__":
    # test_graph()
    test_fci()
    # test_e_pq_on_string()
    
