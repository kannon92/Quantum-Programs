#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "Embed.h"

using namespace boost;

namespace psi{ namespace plugin_fragment {

void Embed::localize_on_fragment(Options& options,boost::shared_ptr<Wavefunction> wfn)
{
    int print = options.get_int("PRINT");
    double localcut = options.get_double("LocalCut");


    // Compute the overlap matrix
    boost::shared_ptr<IntegralFactory> integral_ = wfn->integral();

    boost::shared_ptr<BasisSet> basisset_ = wfn->basisset();
    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S_ao(new Matrix("S_ao",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S_ao);

//    S_ao->print();
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    SharedMatrix AO2SO_ = pet->aotoso();
//    AO2SO_->print();

    // Compute the W matrix for each fragment
    std::vector<int> flist({0});
    std::vector<int> glist;
    boost::shared_ptr<Molecule> frag = mol->extract_subsets(flist,glist);
    int frag_natom = frag->natom();
    //Assuming restricted for now
    outfile->Printf("\n  Fragment contains %d atoms",frag->natom());
    
    // Form a copy of S_ao and zero the rows and columns that are not on this fragment
//    max_a = min_a +
    int nbf = basisset_->nbf();
    int nbfA = 0;
    SharedMatrix S_f(S_ao->clone());
    //This loop is used to get the number of basis functions on a.
    double nEA = 0;
    for (int mu = 0; mu < nbf; mu++) {
        int A = basisset_->function_to_center(mu);
        if (A < frag_natom){
            nbfA += 1;
        }
    }
    double nEB = 0;
    for (int atom = 0; atom < mol->natom(); atom++){
        if(atom < frag_natom){
           nEA+=mol->Z(atom);
        }
        else{ 
           nEB+=mol->Z(atom);
           }

    }
    nEA_ = nEA/2;
    nEB_ = nEB/2;
    nbfA_ = nbfA;
    nbfB_ = nbf - nbfA;
    outfile->Printf( "\nThe number of basis functions on A are %d", nbfA);
    outfile->Printf( "\nThe number of basis functions on B are %d", nbfB_);
    //The overlap matrix for fragment A
    S_f = SharedMatrix(new Matrix("S_f",nbfA,nbfA));
    for (int mu = 0; mu < nbfA; mu++) {
        for (int nu = 0; nu < nbfA; nu++) {
            S_f->set(mu,nu,S_ao->get(mu,nu));
        }
    }
//    S_f->print();
    SharedMatrix S_f_inv(S_f->clone());
    SharedMatrix id(S_f->clone());
    //Inverse of S_f
    S_f_inv->general_invert();
    
    //zero everywhere except for inverse of overlap for fragment A basis functions
    SharedMatrix S_f_inv_large(S_ao->clone());
    S_f_inv_large->zero();
    for (int mu = 0; mu < nbfA; mu++) {
        for (int nu = 0; nu < nbfA; nu++) {
            S_f_inv_large->set(mu,nu,S_f_inv->get(mu,nu));
        }
    }
//    S_f_inv_large->print();
//    Cat*S-1t*S_ao*SCa
    SharedMatrix Ca = wfn->Ca();
    S_f_inv_large->transform(S_ao);
    S_f_inv_large->transform(Ca);
//    S_f_inv_large->print();

    int nalpha = wfn->nalpha();
    outfile->Printf("number of electrons %d", nalpha);

    SharedMatrix SAo(new Matrix("SAo",nalpha,nalpha));
    for (int i = 0; i < nalpha; ++i){
        for (int j = 0; j < nalpha; ++j){
            SAo->set(i,j,S_f_inv_large->get(i,j));
        }
    }
    SharedMatrix Uo(new Matrix("Uo",nalpha,nalpha));
    SharedVector lo(new Vector("lo",nalpha));
    SAo->diagonalize(Uo,lo,descending);
    lo->print();

    int Aocc = 0;
    for(int o = 0; o < nalpha; o++){
       if((lo->get(o))>localcut){
           Aocc++; 
       }
    }
    nAocc_ = Aocc;
    int navir = nbf - nalpha;
    SharedMatrix SAv(new Matrix("SAv",navir,navir));
    for (int i = 0; i < navir; ++i){
        for (int j = 0; j < navir; ++j){
            SAv->set(i,j,S_f_inv_large->get(i + nalpha,j + nalpha));
        }
    }
    SharedMatrix Uv(new Matrix("Uv",navir,navir));
    SharedVector lv(new Vector("lv",navir));
    SAv->diagonalize(Uv,lv,descending);
    lv->print();

    int Bvir = 0;
    for(int v = 0; v < navir; v++){
       if((lv->get(v)) > localcut){
          Bvir++; 
       }
    }
    nAvir_ = Bvir; 
    nBocc_ = nalpha - Aocc;
    nBvir_ = navir - Bvir;

    SharedMatrix U(new Matrix("U",nbf,nbf));
    for (int i = 0; i < nalpha; ++i){
        for (int j = 0; j < nalpha; ++j){
            U->set(i,j,Uo->get(i,j));
        }
    }
    for (int a = 0; a < navir; ++a){
        for (int b = 0; b < navir; ++b){
            U->set(a + nalpha,b + nalpha,Uv->get(a,b));
        }
    }
    U->print();
    SharedMatrix Ca_new(Ca->clone());
    Ca_new->gemm(false, false, 1.0,Ca,U, 0.0);
    Ca->copy(Ca_new);


    //This checks to see if orbitals are or
    S_ao->transform(Ca);
    S_ao->print();
    Ca->print();
    boost::shared_ptr<Matrix> DA(new Matrix("DA", nbf, nbf));
    boost::shared_ptr<Matrix> DB(new Matrix("DB", nbf, nbf));
    boost::shared_ptr<Matrix> Df(new Matrix("DF", nbf, nbf));

    for(int i = 0; i < nbf; i++){
       for(int j = 0; j < nbf; j++){
          for(int nae = 0; nae < nAocc_; nae++){
             DA->add(i,j,Ca->get(i,nae)*Ca->get(j,nae));
             Df->add(i,j,Ca->get(i,nae)*Ca->get(j,nae));
          }
          for(int nbe = nAocc_; nbe < nBocc_ + nAocc_; nbe++){
             DB->add(i,j,Ca->get(i,nbe)*Ca->get(j,nbe));
             Df->add(i,j,Ca->get(i,nbe)*Ca->get(j,nbe));
          }
       }
    }
Df->print();
}

}} // End namespaces




