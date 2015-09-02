#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include "ao_class.h"
#include <lib3index/denominator.h>
#include <libpsio/psio.h>

namespace psi{namespace aomp2 {

ao_class::ao_class(boost::shared_ptr<Wavefunction> wavefunction_, Options& options) 
    : Wavefunction(options, _default_psio_lib_)
{
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
     
    common_init();
}
void ao_class::common_init()
{
    
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Matrix> AOMO = wfn->Ca_subset("AO", "ALL");
    
    AOMO_ = AOMO;

}
void ao_class::compute_mp2()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    get_laplace_factors();
    compute_psuedo();
   
    if(options_.get_str("INTEGRALS")=="CONV")
    {
        boost::shared_ptr<conventional_integrals> conv_ints(new conventional_integrals());
        Full_AO_Ints_ = conv_ints->return_integrals();
    }
    double MP2_full_mo, MP2_full_laplace;
    MP2_full_mo = compute_mp2_no_approx();
    MP2_full_laplace = compute_mp2_laplace_denom();
    outfile->Printf("\n MP2_full_mo:%8.8f \n MP2_full_laplace:%8.8f", MP2_full_mo, MP2_full_laplace);
}
double ao_class::compute_mp2_no_approx()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Vector> eps_occ = wfn->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wfn->epsilon_a_subset("AO", "VIR");
    int nbf = naocc_ + navir_;
    outfile->Printf("\n naocc_:%d navir_:%d", naocc_, navir_);
    double MP2_energy = 0.0;
    for(int i = 0; i < naocc_; i++){
       for(int j = 0; j < naocc_; j++){
          for(int a = 0; a < navir_; a++){
             for(int b = 0; b < navir_; b++){
                       double iajb = Full_AO_Ints_->get(i*nbf + (a + naocc_), j*nbf + (b + naocc_));
                       double ibja = Full_AO_Ints_->get(i*nbf + (b + naocc_), j*nbf + (a + naocc_));
                       double denom = -1.0 / (eps_vir->get(a) + eps_vir->get(b) - eps_occ->get(i) - eps_occ->get(j));
                       MP2_energy+= 2.0*(iajb*iajb)*denom - iajb*ibja*denom;
             }
          }
       }
    }
    return MP2_energy;
}
double ao_class::compute_mp2_laplace_denom()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Vector> eps_occ = wfn->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wfn->epsilon_a_subset("AO", "VIR");
    int nbf = naocc_ + navir_;
    outfile->Printf("\n naocc_:%d navir_:%d", naocc_, navir_);
    double MP2_energy = 0.0;
    int weights = Laplace_occ_->rowspi()[0];
    outfile->Printf("\n weights_occ:%d weights_vir:%d", weights, Laplace_vir_->rowspi()[0]);
    for(int w = 0; w < weights; w++){
        for(int i = 0; i < naocc_; i++){
           for(int j = 0; j < naocc_; j++){
              for(int a = 0; a < navir_; a++){
                 for(int b = 0; b < navir_; b++){
                           double iajb = Full_AO_Ints_->get(i*nbf + (a + naocc_), j*nbf + (b + naocc_));
                           double ibja = Full_AO_Ints_->get(i*nbf + (b + naocc_), j*nbf + (a + naocc_));
                           double denom = -1 * Laplace_occ_->get(w,i) * Laplace_occ_->get(w,j) * Laplace_vir_->get(w,a) * Laplace_vir_->get(w,b);
                           MP2_energy+= 2.0*(iajb*iajb)*denom - iajb*ibja*denom;
                 }
              }
           }
        }
    }
    return MP2_energy;
}
void ao_class::get_laplace_factors()
{
    boost::shared_ptr<Wavefunction> wavefunction_ = Process::environment.wavefunction();
    boost::shared_ptr<Vector> eps_occ = wavefunction_->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wavefunction_->epsilon_a_subset("AO", "VIR");
    naocc_ = eps_occ->dim();
    navir_ = eps_vir->dim();
    double laplace_tolerance = options_.get_double("LAPLACE_TOLERANCE");

    boost::shared_ptr<LaplaceDenominator> laplace(new LaplaceDenominator(eps_occ, eps_vir, laplace_tolerance));

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
    boost::shared_ptr<Wavefunction> wavefunction_ = Process::environment.wavefunction();
    boost::shared_ptr<Matrix> AOMO = wavefunction_->Ca_subset("AO", "ALL");
    integrals_ = Mints->mo_eri(AOMO,AOMO);
}

}}
