#include <libmints/matrix.h>
#include <libmints/vector.h>
#include <libmints/sieve.h>
#include "ao_helper.h"
#include <lib3index/denominator.h>
#include <libfock/jk.h>
namespace psi { namespace aomp2 {

AtomicOrbitalHelper::AtomicOrbitalHelper(SharedMatrix CMO, SharedVector eps_occ, SharedVector eps_vir, double laplace_tolerance)
: CMO_(CMO), eps_rdocc_(eps_occ), eps_virtual_(eps_vir), laplace_tolerance_(laplace_tolerance)
{
    LaplaceDenominator laplace(eps_rdocc_, eps_virtual_, laplace_tolerance_);
    Occupied_Laplace_ = laplace.denominator_occ();
    Virtual_Laplace_  = laplace.denominator_vir();
    weights_          = Occupied_Laplace_->rowspi()[0];
    nrdocc_           = eps_rdocc_->dim();
    nvir_             = eps_virtual_->dim();
    nbf_              = CMO_->rowspi()[0];
}
AtomicOrbitalHelper::~AtomicOrbitalHelper()
{
    outfile->Printf("\n Done with AO helper class");
}
void AtomicOrbitalHelper::Compute_Psuedo_Density()
{
    int nmo_ = nbf_;
    SharedMatrix Xocc(new Matrix("DensityOccupied", weights_, nbf_* nbf_));
    SharedMatrix Yvir(new Matrix("DensityVirtual", weights_, nmo_* nmo_));

    double value_occ, value_vir = 0;
    for(int w = 0; w < weights_; w++){
        for(size_t mu = 0; mu < nbf_; mu++){
            for(size_t nu = 0; nu < nbf_; nu++){
                for(size_t i = 0; i < nrdocc_; i++){
                    value_occ += CMO_->get(mu, i) * CMO_->get(nu, i) * Occupied_Laplace_->get(w, i);

                }
                Xocc->set(w, mu * nmo_ + nu, value_occ);
                for(size_t a = 0; a < nvir_; a++){
                    value_vir += CMO_->get(mu, nrdocc_ + a) * CMO_->get(nu, nrdocc_ + a) * Virtual_Laplace_->get(w, a);
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
void AtomicOrbitalHelper::Compute_AO_Screen(boost::shared_ptr<BasisSet>& primary)
{
    ERISieve sieve(primary, 1e-10);
    std::vector<double> my_function_pair_values = sieve.function_pair_values();
    SharedMatrix AO_Screen(new Matrix("Z", nbf_, nbf_));
    for(int mu = 0; mu < nbf_; mu++)
        for(int nu = 0; nu < nbf_; nu++)
            AO_Screen->set(mu, nu, my_function_pair_values[mu * nbf_ + nu]);

    AO_Screen_ = AO_Screen;
    AO_Screen_->set_name("ScwartzAOInts");
}
void AtomicOrbitalHelper::Estimate_TransAO_Screen(boost::shared_ptr<BasisSet>& primary, boost::shared_ptr<BasisSet>& auxiliary)
{
    Compute_Psuedo_Density();
    boost::shared_ptr<JK> jk(new DFJK(primary, auxiliary));
    jk->initialize();
    jk->compute();
    SharedMatrix AO_Trans_Screen(new Matrix("AOTrans", weights_, nbf_ * nbf_));

    for(int w = 0; w < weights_; w++)
    {
        SharedMatrix COcc(new Matrix("COcc", nbf_, nbf_));
        SharedMatrix CVir(new Matrix("COcc", nbf_, nbf_));
        for(int mu = 0; mu < nbf_; mu++)
            for(int nu = 0; nu < nbf_; nu++)
            {
                COcc->set(mu, nu, POcc_->get(w, mu * nbf_ + nu));
                CVir->set(mu, nu, PVir_->get(w, mu * nbf_ + nu));
            }
                
        SharedVector iaia_w = jk->iaia(COcc, CVir);
        for(int mu = 0; mu < nbf_; mu++)
            for(int nu = 0; nu < nbf_; nu++)
                AO_Trans_Screen->set(w, mu * nbf_ + nu,iaia_w->get(mu * nbf_ + nu));
    }
    TransAO_Screen_ = AO_Trans_Screen;
    TransAO_Screen_->set_name("(u_b {b}^v | u_b {b}^v)");
}

}}
