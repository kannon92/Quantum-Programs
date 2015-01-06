#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include "Embed.h"
using namespace boost;

namespace psi{namespace plugin_fragment {

void Embed::ManbyEmbed()
{
   boost::shared_ptr<Molecule> molecule = Process::environment.molecule();  
   boost::shared_ptr<Wavefunction> wfn  = Process::environment.wavefunction();

   boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());

   boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser,molecule, "BASIS");

   boost::shared_ptr<IntegralFactory> integral(new IntegralFactory (aoBasis, aoBasis, aoBasis, aoBasis));
   SharedVector tei(new Vector("tei array", nmo*nmo*nmo*nmo));


   //Computes all the TEI for right now

   boost::shared_ptr<TwoBodyAOInt> eri(integral->eri()); 

   const double *buffer = eri->buffer();

   AOShellCombinationsIterator shellIter = integral->shells_iterator();

   int count = 0;

   for(shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
      if(eri->compute_shell(shellIter)) {
         AOIntegralsIterator intIter = shellIter.integrals_iterator();
         for(intIter.first(); intIter.is_done() == false; intIter.next()) {
            int p = intIter.i();
            int q = intIter.j();
            int r = intIter.k();
            int s = intIter.l();
            outfile->Printf("\t(%2d %2d | %2d %2d) = %20.15f\n",
                p, q, r, s, buffer[intIter.index()]);
            ++count;

            tei->set(INDEX4(p,q,r,s), buffer[intIter.index()]);
         }
      }
   }
   
   outfile->Printf("I am about to build HAinB matrix");
   int nbf = nbfA_ + nbfB_;
   //Now I will try to build the HAinB matrix
   SharedMatrix HAinB(new Matrix("HAinB", nbf, nbf));
   SharedMatrix HBinA(new Matrix("HAinB", nbf, nbf));
   //The complete density matrix (wfn->Da())
   SharedMatrix Dfull(new Matrix("Dfull", nbf, nbf));
   //The complete MO coefficients from SCF calc
   SharedMatrix CF(new Matrix("CF", nbf, nbf));
   //The complete core hamiltonian
   SharedMatrix Hc(new Matrix("Hc", nbf, nbf));
   //The density matrix for subsystem A - ca*ca
   SharedMatrix DA(new Matrix("DA", nbf, nbf));
   //Density matrix for subsystem B
   SharedMatrix DB(new Matrix("DB", nbf, nbf));

   CF = wfn->Ca();
   CF->print();

   //Grab the density matrix for system A
   outfile->Printf("NaOcc = %d    NBocc = %d", nAocc_, nBocc_);
   for(int i = 0; i < nbf; i++){
      for(int j = 0; j < nbf; j++){
         for(int nae = 0; nae < nAocc_; nae++){
            DA->add(i,j,CF->get(i,nae)*CF->get(j,nae)); 
         }
      }
   }
   boost::shared_ptr<Matrix> DANO(new Matrix("DANO", nbf, nbf));
   boost::shared_ptr<Matrix> DBNO(new Matrix("DBNO", nbf, nbf));
   boost::shared_ptr<Vector> DANOe(new Vector("DANOe",nbf));
   boost::shared_ptr<Vector> DBNOe(new Vector("DBNOe",nbf));
   DA->diagonalize(DANO, DANOe, descending);
   DANOe->print();

   //Grab the denstiy matrix for system B
   for(int i = 0; i < nbf; i++){
      for(int j = 0; j < nbf; j++){
         for(int nbe = nAocc_; nbe < nAocc_ + nBocc_; nbe++){
           DB->add(i,j, CF->get(i,nbe)*CF->get(j,nbe));
         }
      }
   }
   DA->print();
   DB->print();
   DB->diagonalize(DBNO, DBNOe, descending);
   DBNOe->print();
   
   for(int i = 0; i < nbf; i++){
      for(int  j =0; j < nbf; j++){
        for(int mA = 0; mA < nAocc_; mA++){
           Dfull->add(i,j,CF->get(i,mA)*CF->get(j,mA));
        }
        for(int mB = nAocc_; mB < nAocc_ + nBocc_; mB++){
           Dfull->add(i,j,CF->get(i,mB)*CF->get(j,mB)); 
        }
      }
   }
   Dfull->print();

   Hc = wfn->H();

   HAinB->copy(Hc);
   SharedMatrix JA(new Matrix("JA", nbf, nbf)); 
   for(int i=0; i < nbf; i++){
      for(int j = 0; j < nbf; j++){
         for(int k = 0; k < nbf; k++){
            for(int l = 0; l < nbf; l++){
               int ijkl = INDEX4(i,j,k,l);
               int ikjl = INDEX4(i,k,j,l);
               HAinB->add(i,j,Dfull->get(k,l)*(2*tei->get(ijkl) - tei->get(ikjl)));
               HBinA->add(i,j, Dfull->get(k,l)*(2*tei->get(ijkl) - tei->get(ikjl)));
               JA->add(i,j, DA->get(k,l)*(2*tei->get(ijkl)));
               //HAinB->add(i,j, DB->get(k,l)*(2*tei->get(ijkl) - tei->get(ikjl)));
            }
         }
      }  
    }
    JA->print();
    for(int i = 0; i < nbf; i++){
      for(int j = 0; j < nbf; j++){
         for(int k2 = 0; k2 < nbf; k2++){
            for(int l2 = 0; l2 < nbf; l2++){
               int ijkl2 = INDEX4(i,j,k2,l2);
               int ikjl2 = INDEX4(i,k2,j,l2);
               HAinB->add(i,j,-1.0*DA->get(k2,l2)*(2*tei->get(ijkl2)-tei->get(ikjl2)));
               HBinA->add(i,j,-1.0*DB->get(k2,l2)*(2*tei->get(ijkl2)-tei->get(ikjl2)));
            } 
         } 
      }
   }
   HAinB->print();
   boost::shared_ptr<Matrix> P(new Matrix("Proj", nbf, nbf));
   boost::shared_ptr<Matrix> S(new Matrix("overlap", nbf, nbf));
   boost::shared_ptr<Matrix> Sm(new Matrix("small overlap", nbf, nbf));
   boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap()); 
   boost::shared_ptr<Matrix> Shalf(new Matrix("S^-.5", nbf, nbf));
   sOBI->compute(S);
   Shalf->copy(S);
   Shalf->power(-1.0/2);
   P->copy(DB);
   P->transform(S);
   P ->print(); 

   double Etrace = 0.0;
   for(int i = 0; i < nbf; i++){
      for(int j = 0; j < nbf; j++){
         Etrace+=DA->get(i,j)*P->get(i,j);
      }
   }
  
   outfile->Printf("\nEtrace = %20.12f\n", Etrace);

   boost::shared_ptr<Matrix> FA(new Matrix("FA", nbf, nbf));         
   boost::shared_ptr<Matrix> FB(new Matrix("FB", nbf, nbf));         
   boost::shared_ptr<Matrix> FAB(new Matrix("FB", nbf, nbf));         
   boost::shared_ptr<Matrix> Demb(new Matrix("Demb", nbf, nbf));
   boost::shared_ptr<Matrix> Dold(new Matrix("Dold", nbf, nbf));
   boost::shared_ptr<Matrix> Co(new Matrix("Co", nbf, nbf));
   boost::shared_ptr<Matrix> Ci(new Matrix("final Ce", nbf, nbf));
   boost::shared_ptr<Vector> Eeo(new Vector("Epse", nbf));

   boost::shared_ptr<Matrix> Ftrans(new Matrix("Ftrans", nbf, nbf));
   Demb->copy(DA); 
   int iter = 0; double rms = 0.0; double Enew = 0.0;
   double Eold = 0.0;
   double EAB = 0.0;
  
   //These codes create a new Fock matrix for A, B, and the full system AB.  
   FA->copy(HAinB); 
   FB->copy(Hc); 
   FAB->copy(Hc); 
   for(int i = 0; i < nbf; i++){
        for(int j = 0; j < nbf; j++){
             for(int k = 0; k < nbf; k++){
                for(int l = 0; l < nbf; l++){
                   int ijkl = INDEX4(i,j,k,l);
                   int ikjl = INDEX4(i,k,j,l);
                   FA->add(i,j,Demb->get(k,l)*(2*tei->get(ijkl) - tei->get(ikjl)));
                   FB->add(i,j,DB->get(k,l)*(2*tei->get(ijkl) - tei->get(ikjl)));
                   FAB->add(i,j,Dfull->get(k,l)*(2*tei->get(ijkl) - tei->get(ikjl)));
                }
             }
          }
       }

       //This is the energy of the system of A.  tr(DA(HAinB + FA))
       for(int i = 0; i < nbf; i++){
          for(int j = 0; j < nbf; j++){
             Enew += Demb->get(i,j)*(HAinB->get(i,j) + FA->get(i,j));
             EAB += Dfull->get(i,j)*(Hc->get(i,j) + FAB->get(i,j));
          }
       }
   double E_emb = Enew;
   double EembB = 0.0;
   double Etr = 0.0;
   double Enonadd = 0.0;

   //Computes the energy of system B
   for(int i = 0; i < nbf; i++){
      for(int j = 0; j < nbf; j++){
         EembB+=DB->get(i,j)*(Hc->get(i, j)+FB->get(i,j));     
      }
   }
   outfile->Printf("EembB = %20.12f", EembB);
   
   double nucrep = molecule->nuclear_repulsion_energy();
   outfile->Printf("\nEmbed Energy = %20.12f\n", E_emb + EembB + Etr + Enonadd + nucrep);
   double scfenergy = wfn->reference_energy();
   outfile->Printf("\n SCF energy = %20.12f\n",scfenergy );
   outfile->Printf("\n Error = %20.12f \n", scfenergy - (E_emb + EembB + Enonadd + nucrep + Etr));
   outfile->Printf("\nEAB = %20.12f\n",((EAB + nucrep) - scfenergy)); 
   
   FockEmbed();

}

}}
