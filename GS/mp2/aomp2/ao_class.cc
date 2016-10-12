#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include "ao_class.h"
#include <lib3index/denominator.h>
#include <libpsio/psio.h>
#include <ambit/tensor.h>
#include "ao_helper.h"
//#include <ctf.hpp>
//#include <mpi.h>

namespace psi{namespace aomp2 {

ao_class::ao_class(SharedWavefunction wavefunction, Options& options) 
    : wavefunction_(wavefunction), options_(options)
{
    
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    //ambit::initialize();
     
    common_init();

}
void ao_class::common_init()
{
    
    boost::shared_ptr<Matrix> AOMO = wavefunction_->Ca_subset("AO", "ALL");
    
    AOMO_ = AOMO;
    nmo_ = wavefunction_->nmo();

}
void ao_class::compute_mp2()
{
    get_laplace_factors();
    compute_psuedo();
   
    if(options_.get_str("INTEGRALS")=="CONV")
    {
        boost::shared_ptr<conventional_integrals> conv_ints(new conventional_integrals(wavefunction_));
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
    boost::shared_ptr<Vector> eps_occ = wavefunction_->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wavefunction_->epsilon_a_subset("AO", "VIR");
    int nbf = naocc_ + navir_;
    outfile->Printf("\n naocc_:%d navir_:%d", naocc_, navir_);
    double MP2_energy = 0.0;
    int count = 0;
    for(size_t i = 0; i < naocc_; i++){
       for(size_t j = 0; j < naocc_; j++){
          for(size_t a = 0; a < navir_; a++){
             for(size_t b = 0; b < navir_; b++){
                       double iajb = Full_MO_Ints_->get(i*nbf + (a + naocc_), j*nbf + (b + naocc_));
                       double ibja = Full_MO_Ints_->get(i*nbf + (b + naocc_), j*nbf + (a + naocc_));
                       if(iajb < 1e-10 or ibja < 1e-10)
                            count++;
                       double denom = -1.0 / (eps_vir->get(a) + eps_vir->get(b) - eps_occ->get(i) - eps_occ->get(j));
                       MP2_energy+= 2.0*(iajb*iajb)*denom - iajb*ibja*denom;
             }
          }
       }
    }

    outfile->Printf("\n Number of zero in MP2 %d", count);
    return MP2_energy;
}
double ao_class::compute_mp2_laplace_denom()
{
    boost::shared_ptr<Vector> eps_occ = wavefunction_->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wavefunction_->epsilon_a_subset("AO", "VIR");
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
void ao_class::test_atomic_orbital_class()
{
    boost::shared_ptr<Vector> eps_occ = wavefunction_->epsilon_a_subset("AO", "OCC");
    boost::shared_ptr<Vector> eps_vir = wavefunction_->epsilon_a_subset("AO", "VIR");
    AtomicOrbitalHelper ao_helper(AOMO_, eps_occ, eps_vir, 1e-10);
    boost::shared_ptr<BasisSet> primary = wavefunction_->basisset();
    boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_orbital(primary->molecule(), "DF_BASIS_MP2", options_.get_str("DF_BASIS_MP2"));
    ao_helper.Compute_AO_Screen(primary);
    SharedMatrix AO_Screen = ao_helper.AO_Screen();
    ao_helper.Estimate_TransAO_Screen(primary, auxiliary);
    SharedMatrix TransAO_Screen = ao_helper.TransAO_Screen();
    AO_Screen->print();
    TransAO_Screen->print();
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
    ambit::Tensor TransAO = ambit::Tensor::build(ambit::CoreTensor, "TransAO", {weights_, nmo_, nmo_, nmo_, nmo_});
    ambit::Tensor BAO = ambit::Tensor::build(ambit::CoreTensor, "BTransAO", {weights_, nmo_, nmo_, nmo_, nmo_});
    ambit::Tensor CAO = ambit::Tensor::build(ambit::CoreTensor, "CTransAO", {weights_, nmo_, nmo_, nmo_, nmo_});
    ambit::Tensor POcc    = ambit::Tensor::build(ambit::CoreTensor, "Pocc", {weights_, nmo_, nmo_});
    ambit::Tensor PVir    = ambit::Tensor::build(ambit::CoreTensor, "Pocc", {weights_, nmo_, nmo_});
    ambit::Tensor AOFull = ambit::Tensor::build(ambit::CoreTensor, "TransAO", {nmo_, nmo_, nmo_, nmo_});

    AOFull.iterate([&](const std::vector<size_t>& i,double& value){
        value = Full_AO_Ints_->get(i[0] * nmo_ + i[1], i[2] * nmo_ + i[3]);});

    POcc.iterate([&](const std::vector<size_t>& i,double& value){
        value = POcc_->get(i[0], i[1] * nmo_ + i[2]);});
    PVir.iterate([&](const std::vector<size_t>& i,double& value){
        value = PVir_->get(i[0], i[1] * nmo_ + i[2]);});

    TransAO("w, mu, nu, lam, si") = POcc("w, mu, gam") * PVir("w, nu, del") * AOFull("gam, del, kap, eps") * POcc("w, kap, lam") * PVir("w, eps, si");
    BAO("w, mu, nu, lam, si") = POcc("w, mu, gam") * AOFull("gam, nu, kap, si") * POcc("w, kap, lam");
    CAO("w, mu, nu, kap, si") = PVir("w, nu, gam") * AOFull("mu, gam, kap, lam") * PVir("w, lam, si");

    ambit::Tensor E = ambit::Tensor::build(ambit::CoreTensor, "EAlpha", {weights_});
    ambit::Tensor AOFullV = ambit::Tensor::build(ambit::CoreTensor, "TransAO", {nmo_, nmo_, nmo_, nmo_});
    AOFullV("mu, nu, la, si") = 2.0 * AOFull("mu, nu, la, si");
    AOFullV("mu, nu, la, si") -= AOFull("mu, si, la, nu");
    SharedMatrix AOFullV_screen(new Matrix("AOV", nmo_, nmo_));
    AOFullV.iterate([&](const std::vector<size_t>& i,double& value){
        if(i[0] == i[2] and i[1] == i[3])
        {
        AOFullV_screen->set(i[0], i[1],value);
        }
        });
    //E("w") = TransAO("w, mu, nu, la, si") * (2 * AOFull("mu, nu, la, si"));
    //E("w") -= TransAO("w, mu, nu, la, si") * AOFull("mu, si, la, nu");
    E("w") = TransAO("w, mu, nu, la, si") * (AOFullV("mu, nu, la, si"));
    //E.print(stdout);
    int count = 0;
    TransAO.iterate([&](const std::vector<size_t>& i,double& value){
        if(value < 1e-10) count++;});

    //outfile->Printf("\n Number of zero in AO-MP2 %d", count);

    SharedMatrix A_screen(new Matrix("(uv|uv)", nmo_, nmo_));
    SharedMatrix B_screen(new Matrix("(u_bv|u_bv)", weights_ , nmo_ * nmo_));
    SharedMatrix C_screen(new Matrix("(uv_b|uv_b)", weights_ , nmo_ * nmo_));
    SharedMatrix D_screen(new Matrix("D", weights_ , nmo_ * nmo_));

    
    std::vector<double>& B_vec = BAO.data();
    std::vector<double>& C_vec = CAO.data();
    SharedMatrix BV(new Matrix("B * PV", weights_ , nmo_ * nmo_));
    SharedMatrix CO(new Matrix("C * PO", weights_ , nmo_ * nmo_));


    for(int mu = 0; mu < nmo_; mu++){
        for(int nu = 0; nu < nmo_; nu++){
            A_screen->set(mu, nu, sqrt(fabs(Full_AO_Ints_->get(mu * nmo_ + nu, mu * nmo_ + nu))));
            for(int w = 0; w < weights_; w++){
                B_screen->set(w, mu * nmo_ + nu, sqrt(fabs(B_vec[w * nmo_ * nmo_ * nmo_ * nmo_ + mu * nmo_ * nmo_ * nmo_ + nu * nmo_ * nmo_ + mu * nmo_ + nu]))); 
                C_screen->set(w, mu * nmo_ + nu, sqrt(fabs(C_vec[w * nmo_ * nmo_ * nmo_ * nmo_ + mu * nmo_ * nmo_ * nmo_ + nu * nmo_ * nmo_ + mu * nmo_ + nu]))); 
                double valueB, valueC;
                valueB = 0.0;
                valueC = 0.0;
                for(int sig = 0; sig < nmo_; sig++) 
                {
                    valueB += B_screen->get(w , mu * nmo_ + sig) * fabs(PVir_->get(w, sig * nmo_ + nu));
                    valueC += C_screen->get(w , mu * nmo_ + sig) * fabs(POcc_->get(w, sig * nmo_ + nu));
                }
                BV->set(w ,mu * nmo_ +  nu, valueB);
                CO->set(w ,mu * nmo_ +  nu, valueC);
                double value_min = 0.0;
                if(valueB < valueC)
                    value_min = valueB;
                else 
                    value_min = valueC;
                D_screen->set(w , mu * nmo_ + nu, value_min);

            }
        }
    }
    //if(BV->rms() < CO->rms())
    //    D_screen = BV;
    //else 
    //    D_screen = CO;
    //A_screen->print();
    //D_screen = BV;
    //B_screen->print();
    //D_screen->print();
    //BV->print();
    //C_screen->print();
    //POcc_->set_name("C_occ C_occ^T exp(shit)");
    //POcc_->print();
    //CO->print();
    //PVir_->print();

    //A_screen = AOFullV_screen;
    test_atomic_orbital_class();
    Full_AO_Ints_->print();
    AOFullV_screen->print();
    A_screen->print();
    double screen_energy = 0.0;
    int count_screen = 0;
    int shoulda_screen = 0;
    for(int mu = 0; mu < nmo_; mu++) {
        for(int nu = 0; nu < nmo_; nu++) {
            //if(A_screen->get(mu, nu) > 1e-10)
            {
                for(int rho = 0; rho < nmo_; rho++) {
                    for(int sig = 0; sig < nmo_; sig++) {
                        //if(A_screen->get(rho, sig) > 1e-10)
                        {
                            for(int w = 0; w < weights_; w++) {

                                //if(D_screen->get(w, mu * nmo_ +  nu) * D_screen->get(w , rho * nmo_ + sig) > 1e-10)
                                //if(A_screen->get(mu, nu) * A_screen->get(rho, sig) > 1e-10)
                                {
                                    double value_w = TransAO.data()[w * nmo_ * nmo_ * nmo_ * nmo_ + mu * nmo_ * nmo_ * nmo_ + nu * nmo_ * nmo_ + rho * nmo_ + sig];
                                    double value_ao = AOFullV.data()[mu * nmo_ * nmo_ * nmo_ + nu * nmo_ * nmo_ + rho * nmo_ + sig];
                                    if(fabs(value_w) < 1e-10 or fabs(value_ao) < 1e-10)
                                    {
                                        outfile->Printf("\n mu: %d nu: %d rho: %d sig: %d TransAO: %8.14f value_ao: %8.14f", mu, nu, rho, sig, value_w,
                                        value_ao);
                                        shoulda_screen++;
                                        outfile->Printf("\n A_screen(mu, nu): %8.10f A_screen(rho, sig): %8.10f", A_screen->get(mu,
                                        nu),A_screen->get(rho, sig));
                                    }
                                    else {
                                        screen_energy += -1 * value_w * value_ao;
                                    }

                                    count_screen++;
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    outfile->Printf("\n Skipped %d out of %d actual zeros: %d", -1 * count_screen + weights_ * nmo_ * nmo_ * nmo_ * nmo_, weights_ * nmo_ * nmo_
    * nmo_ * nmo_, shoulda_screen);


    double mp2_ao = 0.0;

    for(int w = 0; w < weights_; w++)
    { mp2_ao += -1 * E.data()[w];}
    outfile->Printf("\n screen_energy: %8.8f mp2_ao: %8.8f", screen_energy, mp2_ao);
    double difference = screen_energy - mp2_ao;
    //if(fabs(screen_energy - mp2_ao) < 1e-5)
    //    throw PSIEXCEPTION("AO-MP2 energies do not agree");

return mp2_ao;

}
void SchwartzScreen()
{
   outfile->Printf("Going to do something shortly"); 
}


conventional_integrals::conventional_integrals(SharedWavefunction ref_wfn)
: wavefunction_(ref_wfn)
{
    MintsHelper Mints(wavefunction_);
    boost::shared_ptr<Matrix> AOMO = wavefunction_->Ca_subset("AO", "ALL");
    mo_integrals_ = Mints.mo_eri(AOMO,AOMO);
    ao_integrals_ = Mints.ao_eri();
}

}}
