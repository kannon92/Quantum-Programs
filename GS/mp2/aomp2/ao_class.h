#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/vector.h>
#include <libmints/wavefunction.h>


namespace psi { namespace aomp2 {

class ao_class : public Wavefunction
{
private:
    //The matrix containing ( pq | rs) integrals (not antisymmetrized)
    boost::shared_ptr<Matrix> Full_AO_Ints_;
    //The matrix containing the laplace factors for occupied orbitals
    boost::shared_ptr<Matrix> Laplace_occ_;
    //The matrix containing the laplace factors for virtual orbitals
    boost::shared_ptr<Matrix> Laplace_vir_;
    //The occupied psudodensity matrix
    //P_{/mu /nu} = \omega \sum_i^occ CCe^{-\epsilon_i t}
    boost::shared_ptr<Matrix> PsuedoDOcc_;
    //P_{/mu /nu} = \omega \sum_a^virt CCe^{-\epsilon_i t}
    boost::shared_ptr<Matrix> PsuedoDVir_;
    //The AO->MO coefficient matrix.  Get from SCF 
    boost::shared_ptr<Matrix>  AOMO_;
    //The wavefunction to get init stuff
    boost::shared_ptr<Wavefunction> wavefunction_;
    //size of occupied orbitals
    int naocc_;
    //size of virtual orbitals
    int navir_;
    //size of nmo
    int nmo_;
    
    //Compute the psuedo matrices
    void compute_psuedo();
    //Compute the occupied and virtual laplace factors
    void get_laplace_factors();
    double compute_energy(){return 0.0;}
    double compute_mp2_no_approx();
    double compute_mp2_laplace_denom();
    double compute_ao_mp2();


public:
    ao_class(boost::shared_ptr<Wavefunction> WFN, Options &options);
    void common_init();
    void compute_mp2();


};
class conventional_integrals 
{
private:
    boost::shared_ptr<Matrix> integrals_;
    
public:
    conventional_integrals();
    boost::shared_ptr<Matrix> return_integrals()
    {return integrals_;}
};

class df_integrals
{
private:
    boost::shared_ptr<BasisSet> primary_;
    boost::shared_ptr<BasisSet> auxiliary_;
    boost::shared_ptr<Matrix> uvQ_;
public:
    df_integrals();
};
class cholesky_integrals
{
private:
    boost::shared_ptr<Matrix> L_;
public:
    cholesky_integrals();
};

}}

