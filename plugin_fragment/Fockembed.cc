#include <psi4-dec.h>
#include <psifiles.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/matrix.h>
#include "Embed.h"
#include <libfock/jk.h>
#include <libpsio/psio.hpp>
//##include <libpsio/psio.h>
using namespace boost;

namespace psi{namespace plugin_fragment {

void Embed::FockEmbed()
{
     boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
     boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
     
     boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());

     boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");
     boost::shared_ptr<IntegralFactory> integral(new IntegralFactory (aoBasis, aoBasis, aoBasis, aoBasis));

     //Since I need three different J, K matrices for each system, use JK library for this
     boost::shared_ptr<JK> jkA = JK::build_JK();
     boost::shared_ptr<JK> jkAB = JK::build_JK();
     boost::shared_ptr<JK> jkB = JK::build_JK();
   
     jkA->set_memory(1000000000L);
     jkAB->set_memory(1000000000L);
     jkB->set_memory(1000000000L);
     jkA->set_cutoff(1.0E-12);
     jkAB->set_cutoff(1.0E-12);
     jkB->set_cutoff(1.0E-12);
     jkA->initialize();
     jkAB->initialize();
     jkB->initialize();
     jkA->print_header();
     jkAB->print_header();
     jkB->print_header();
    
     int nbf = nbfA_ + nbfB_;
     boost::shared_ptr<MatrixFactory> HAinB_fact(new MatrixFactory);
     HAinB_fact->init_with(nbf, nbf);
    
     SharedMatrix CF(new Matrix("CF", nbf, nbf));
     CF = wfn->Ca();
     CF->print();
   
     //The C matrix for system A, B, and AB.  
     SharedMatrix CAocc(new Matrix("CAocc", nbf, nAocc_));
     SharedMatrix CBocc(new Matrix("CBocc", nbf, nBocc_));
     SharedMatrix CABocc(new Matrix("CBocc", nbf, nAocc_+nBocc_));
     for(int i = 0; i < nbf; i++){
        for(int j = 0; j < nAocc_; j++){
           CAocc->set(i,j,CF->get(i,j)); 
           CABocc->set(i,j,CF->get(i,j));
        }
        for(int k = 0; k < nBocc_; k++){
           CBocc->set(i,k,CF->get(i,k + nAocc_));
           CABocc->set(i,k + nAocc_,CF->get(i,k + nAocc_));
        }
     }

     outfile->Printf("\nnAocc = %d\t\n nBocc = %d \n", nAocc_, nBocc_);
     CAocc->print();
     CBocc->print();
     CABocc->print();
     
     SharedMatrix DA(new Matrix("DA", nbf, nbf));
     SharedMatrix DB(new Matrix("DB", nbf, nbf));
     SharedMatrix DAB(new Matrix("DAB", nbf, nbf));
     SharedMatrix GA(new Matrix("GA", nbf, nbf));
     SharedMatrix GAB(new Matrix("GAB", nbf, nbf));
     SharedMatrix GB(new Matrix("GB", nbf, nbf));

     //A vector of matrices - used to generate J and K
     //This is how I can do do different calculations since the only thing
     //that changes is the density matrix
     std::vector<boost::shared_ptr<Matrix> >& Cl = jkA->C_left();
     std::vector<boost::shared_ptr<Matrix> >& CAB = jkAB->C_left();
     std::vector<boost::shared_ptr<Matrix> >& CB = jkB->C_left();

     Cl.clear();
     CAB.clear();
     CB.clear();
     //Creates CAocc, CABocc, CBocc
     Cl.push_back(CAocc);
     CAB.push_back(CABocc);
     CB.push_back(CBocc);

     jkA->compute();
     jkAB->compute();
     jkB->compute();
     DA = jkA->D()[0];
     DAB = jkAB->D()[0];
     DB = jkB->D()[0];

     //The code below forms G = 2J - K for A, B, and AB.  
     SharedMatrix JA = jkA->J()[0];
     SharedMatrix JAB = jkAB->J()[0];
     SharedMatrix JB = jkB->J()[0];
     SharedMatrix KA = jkA->K()[0];
     SharedMatrix KAB = jkAB->K()[0];
     SharedMatrix KB = jkB->K()[0];
     JA->scale(2.0);
     JAB->scale(2.0);
     JB->scale(2.0);
     GA->copy(JA);
     GA->subtract(KA);
     GA->print();
     GAB->copy(JAB);
     GAB->subtract(KAB);
     GB->copy(JB);
     GB->subtract(KB);
 
     double EA = 0.0;
     double EB = 0.0;
     double EAB = 0.0;
     SharedMatrix HAinB(new Matrix("PSIF_SO_H", nbf, nbf));
     //Forms HAinB by passing GAB and GA
     // HAinB = Hc + 2J_a - K_a -(2J_ab - k_ab)
     HAinB = formHAinB(GAB, GA);
     HAinB->print();
     SharedMatrix Hc(new Matrix("Hc", nbf, nbf));
     Hc = wfn->H();
     //HAinB->zero();

     double nucrep = molecule->nuclear_repulsion_energy();
     //EA = tr[DA(HAinB + GA)]
     EA = EmbedEnergy(DA,HAinB,GA);
     //EB = tr(DB(Hc + GB))
     EB = EmbedEnergy(DB,Hc,GB);
     EAB = EmbedEnergy(DAB,Hc,GAB);
     //EB = tr(DAB(Hc + GAB)) --> SCF_energy 
     //PSIF_SO_H
     SharedMatrix HAinB_write(HAinB_fact->create_matrix(PSIF_SO_H));

     HAinB_write->copy(HAinB);
     HAinB_write->save(_default_psio_lib_, PSIF_OEI);

     outfile->Printf("\n\n\tEnergy of system A = %20.12f \n\n\t", EA);
     outfile->Printf("\n\n\tEnergy of system B = %20.12f \n\n\t", EB);
     outfile->Printf("\n\n\t Embedded energy = %20.12f \n \n\t", EA + EB + nucrep);
     outfile->Printf("\n\n\tEnergy of system  = %20.12f \n\n\t", EAB + nucrep);
     outfile->Printf("\n\n\tError  = %20.12f \n\n\t", (EA + EB + nucrep) - wfn->reference_energy());
     
}
boost::shared_ptr<Matrix> Embed::formHAinB(boost::shared_ptr<Matrix>& GAB,boost::shared_ptr<Matrix>& GA){
   int nbf = nbfA_ + nbfB_;
   SharedMatrix H(new Matrix("H", nbf, nbf));
   SharedMatrix F(new Matrix("F", nbf, nbf));
   boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

   H->copy(wfn->H()); 
   H->add(GAB);
   H->subtract(GA);

   return H;
}
double Embed::EmbedEnergy(boost::shared_ptr<Matrix>& D, boost::shared_ptr<Matrix>& H, boost::shared_ptr<Matrix>& G){ 
     double EA = 0.0; 
     int nbf = nbfA_ + nbfB_;
     SharedMatrix F(new Matrix("F", nbf, nbf));
     F->copy(H);
     F->add(G);
     F->add(H);
     EA = D->vector_dot(F);

    return(EA);
}

}}
