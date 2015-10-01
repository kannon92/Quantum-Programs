#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include "ao_class.h"
#include <lib3index/denominator.h>
#include <libpsio/psio.h>
#include <ambit/tensor.h>

namespace psi{namespace aomp2 {

ao_class::ao_class(boost::shared_ptr<Wavefunction> wavefunction_, Options& options) 
    : Wavefunction(options, _default_psio_lib_)
{
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    ambit::initialize();
     
    common_init();

}
void ao_class::common_init()
{
    
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Matrix> AOMO = wfn->Ca_subset("AO", "ALL");
    
    AOMO_ = AOMO;
    nmo_ = wfn->nmo();

}
void ao_class::compute_mp2()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    get_laplace_factors();
    compute_psuedo();
   
    if(options_.get_str("INTEGRALS")=="CONV")
    {
        boost::shared_ptr<conventional_integrals> conv_ints(new conventional_integrals());
        Full_MO_Ints_ = conv_ints->mo_integrals();
        Full_AO_Ints_ = conv_ints->ao_integrals();
    }
    double MP2_full_mo, MP2_full_laplace;
    MP2_full_mo = compute_mp2_no_approx();
    MP2_full_laplace = compute_mp2_laplace_denom();
    double MP2_full_ao_mp2  = compute_mp2_ao_laplace();
    outfile->Printf("\n MP2_full_mo:%8.8f \n MP2_full_laplace:%8.8f MP2_full_ao_mp2:%8.8f", MP2_full_mo, MP2_full_laplace, MP2_full_ao_mp2);

}
double ao_class::compute_mp2_no_approx()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Vector> eps_occ = wfn->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wfn->epsilon_a_subset("AO", "VIR");
    int nbf = naocc_ + navir_;
    outfile->Printf("\n naocc_:%d navir_:%d", naocc_, navir_);
    double MP2_energy = 0.0;
    for(size_t i = 0; i < naocc_; i++){
       for(size_t j = 0; j < naocc_; j++){
          for(size_t a = 0; a < navir_; a++){
             for(size_t b = 0; b < navir_; b++){
                       double iajb = Full_MO_Ints_->get(i*nbf + (a + naocc_), j*nbf + (b + naocc_));
                       double ibja = Full_MO_Ints_->get(i*nbf + (b + naocc_), j*nbf + (a + naocc_));
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
        for(size_t i = 0; i < naocc_; i++){
           for(size_t j = 0; j < naocc_; j++){
              for(size_t a = 0; a < navir_; a++){
                 for(size_t b = 0; b < navir_; b++){
                           double iajb = Full_MO_Ints_->get(i*nbf + (a + naocc_), j*nbf + (b + naocc_));
                           double ibja = Full_MO_Ints_->get(i*nbf + (b + naocc_), j*nbf + (a + naocc_));
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
    weights_ = Laplace_occ_->rowspi()[0];

}

void ao_class::compute_psuedo()
{
    SharedMatrix Xocc(new Matrix("DensityOccupied", weights_, nmo_* nmo_));
    SharedMatrix Yvir(new Matrix("DensityVirtual", weights_, nmo_* nmo_));

    double value_occ, value_vir = 0;
    for(int w = 0; w < weights_; w++){
        for(size_t mu = 0; mu < nmo_; mu++){
            for(size_t nu = 0; nu < nmo_; nu++){
                for(size_t i = 0; i < naocc_; i++){
                    value_occ += AOMO_->get(mu, i) * AOMO_->get(nu, i) * Laplace_occ_->get(w, i);
                    
                }
                Xocc->set(w, mu * nmo_ + nu, value_occ);
                for(size_t a = 0; a < navir_; a++){
                    value_vir += AOMO_->get(mu, naocc_ + a) * AOMO_->get(nu, naocc_ + a) * Laplace_vir_->get(w, a);
                }
                Yvir->set(w, mu * nmo_ + nu, value_vir);
                value_occ = 0.0;
                value_vir = 0.0;
            }
        }
    }
    POcc_ = Xocc->clone();
    PVir_ = Yvir->clone();

}
double ao_class::compute_mp2_ao_laplace()
{
    ambit::Tensor TransAO = ambit::Tensor::build(ambit::kCore, "TransAO", {weights_, nmo_, nmo_, nmo_, nmo_});
    ambit::Tensor POcc    = ambit::Tensor::build(ambit::kCore, "Pocc", {weights_, nmo_, nmo_});
    ambit::Tensor PVir    = ambit::Tensor::build(ambit::kCore, "Pocc", {weights_, nmo_, nmo_});
    ambit::Tensor AOFull = ambit::Tensor::build(ambit::kCore, "TransAO", {nmo_, nmo_, nmo_, nmo_});

    AOFull.iterate([&](const std::vector<size_t>& i,double& value){
        value = Full_AO_Ints_->get(i[0] * nmo_ + i[1], i[2] * nmo_ + i[3]);});

    POcc.iterate([&](const std::vector<size_t>& i,double& value){
        value = POcc_->get(i[0], i[1] * nmo_ + i[2]);});
    PVir.iterate([&](const std::vector<size_t>& i,double& value){
        value = PVir_->get(i[0], i[1] * nmo_ + i[2]);});

    TransAO("w, mu, nu, lam, si") = POcc("w, mu, gam") * PVir("w, nu, del") * AOFull("gam, del, kap, eps") * POcc("w, kap, lam") * PVir("w, eps, si");

    ambit::Tensor E = ambit::Tensor::build(ambit::kCore, "EAlpha", {weights_});
    E("w") = TransAO("w, mu, nu, la, si") * (2 * AOFull("mu, nu, la, si"));
    E("w") -= TransAO("w, mu, nu, la, si") * AOFull("mu, si, la, nu");
    E.print(stdout);

    double mp2_ao = 0.0;

    for(int w = 0; w < weights_; w++)
    { mp2_ao += -1 * E.data()[w];}


return mp2_ao;

}
void SchwartzScreen()
{
   outfile->Printf("Going to do something shortly"); 

    
}


conventional_integrals::conventional_integrals()
{
    boost::shared_ptr<MintsHelper> Mints(new MintsHelper());
    boost::shared_ptr<Wavefunction> wavefunction_ = Process::environment.wavefunction();
    boost::shared_ptr<Matrix> AOMO = wavefunction_->Ca_subset("AO", "ALL");
    mo_integrals_ = Mints->mo_eri(AOMO,AOMO);
    ao_integrals_ = Mints->ao_eri();
}

}}
