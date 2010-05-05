cimport primitive_gto
import primitive_gto

cdef extern from "contracted-gto.h":
    ctypedef struct cContractedGTO "ContractedGTO":
        primitive_gto.cPrimitiveGTO *primitives
        double *coefs
        int nprims
        double norm
    
    cContractedGTO *contracted_gto_new ()
    void contracted_gto_init(cContractedGTO *cgto, double x, double y, double z, int l, int m, int n, int atid )

    void contracted_gto_from_primitives( cContractedGTO *cgto, primitive_gto.cPrimitiveGTO **pgtos, 
                                        int nbasis)		

    void contracted_gto_add_primitive(cContractedGTO *cgto, primitive_gto.cPrimitiveGTO *pgto, 
                                        double coef)

    void contracted_gto_free(cContractedGTO *cgto)
    void contracted_gto_normalize(cContractedGTO *cgto)
    double contracted_gto_overlap(cContractedGTO *cgto1, cContractedGTO *cgto2)
    double contracted_gto_amp(cContractedGTO *cgto, double x, double y, double z)
#    double contracted_gto_coulomb(cContractedGTO *c1,cContractedGTO *c2,
#                                  cContractedGTO *c3, cContractedGTO *c4)
    double contracted_gto_renorm_prefactor(cContractedGTO *c1,cContractedGTO *c2,
                                           cContractedGTO *c3, cContractedGTO *c4)
    

cdef class ContractedGTO:
    cdef cContractedGTO *this
#cpdef double coulomb( ContractedGTO a, ContractedGTO b, ContractedGTO c,  ContractedGTO d )
