#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include "ao_class.h"
#include <lib3index/denominator.h>

namespace psi{namespace aomp2 {

ao_class::ao_class()
{
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    wavefunction_ = Process::environment.wavefunction();
    common_init();
}
void ao_class::common_init()
{
    AOMO_ = wavefunction_->Ca();

}
void ao_class::compute_mp2()
{
    get_laplace_factors();
    compute_psuedo();

    //if(options_.get_str("INTEGRAL")=="CHOL")
    //{
        boost::shared_ptr<conventional_integrals> conv_ints(new conventional_integrals());
    //}
    
}
void ao_class::get_laplace_factors()
{
    //Laplace_occ_
    //Laplace_vir_
    boost::shared_ptr<Vector> eps_occ = wavefunction_->epsilon_a_subset("AO", "OCC");
    eps_occ->print();
    boost::shared_ptr<Vector> eps_vir = wavefunction_->epsilon_a_subset("AO", "VIR");
    eps_vir->print();

    boost::shared_ptr<LaplaceDenominator> laplace(new LaplaceDenominator(eps_occ, eps_vir, 1e-14));

    Laplace_occ_ = laplace->denominator_occ();
    Laplace_vir_ = laplace->denominator_vir();

}

void ao_class::compute_psuedo()
{
    //This is only needed for full AO-transformed.  Lazy right now 
}

conventional_integrals::conventional_integrals()
{
    boost::shared_ptr<MintsHelper> Mints(new MintsHelper());
    integrals_ = Mints->ao_eri();
}

}}
