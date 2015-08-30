#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include "ao_class.h"
#include <libthce/laplace.h>
#include <libthce/thcew.h>
#include <libthce/thce.h>

namespace psi{namespace aomp2 {

ao_class::ao_class( Options& options)
{
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    wavefunction_ = Process::environment.wavefunction();
    common_init();
}
void ao_class::common_init()
{
    AOMO_ = wavefunction_->Ca();

}
void ao_class::begin()
{
    get_laplace_factors();
    
}
void ao_class::get_laplace_factors()
{
    //Laplace_occ_
    //Laplace_vir_
    boost::shared_ptr<Vector> eps_occ = wavefunction_->epsilon_a_subset("AO", "OCC");
    eps_occ->print();
    boost::shared_ptr<Vector> eps_vir = wavefunction_->epsilon_a_subset("AO", "VIR");
    eps_vir->print();
    boost::shared_ptr<LaplaceDenom> laplace (new LaplaceDenom::LaplaceDenom(eps_occ, eps_vir, 1e-6, 0.0, 2));
    //laplace->compute();
    //boost::shared_ptr<Tensor> tau_occ = laplace->tau_occ();
    //boost::shared_ptr<Tensor> tau_vir = laplace->tau_vir();
    //int size1 = tau_occ->active_sizes()[1];
    //int size0 = tau_occ->active_sizes()[0];

    //outfile->Printf("\n %d %d", size0, size1);


}
}}
