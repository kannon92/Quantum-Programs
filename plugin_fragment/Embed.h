#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/wavefunction.h>

#ifndef Embed_H
#define Embed_H

namespace psi { namespace plugin_fragment {
class Embed : public Wavefunction 
{
private:
      boost::shared_ptr<PSIO> psio_;
public:
    
     void localize_on_fragment(Options& options, boost::shared_ptr<Wavefunction> wfn);
     void localize_on_atoms(Options& options, boost::shared_ptr<Wavefunction> wfn);
     //This is here so I can actually use this class
     double compute_energy()
     { return 0.0;}
     //I am using this function to start the embed procedure
     //The constructor of a class. OH MY! 
     Embed(boost::shared_ptr<Wavefunction> wfn, Options& Options);
     //Does things to get nmo, nso, and all that jazz
     void common_init();
     // for use in common_init
     virtual ~Embed();
     long int nso, nmo, ndocc, nvirt, nfzc, nfzv, ndoccact;
     void begin(Options& options);
     //The number of basis functions with fragment A
     int nbfA_ = 0;
     //The number of basis functions with fragment B.
     int nbfB_ = 0;
     int nEA_  = 0;
     int nEB_  = 0;
     int nAocc_ = 0;
     int nBocc_ = 0;
     int nAvir_ = 0;
     int nBvir_ = 0;
     //Deconstructor.  Does nothing right now
     //The density matrix is in the full basis 
     boost::shared_ptr<Matrix> Dfull;
     //the density matrix for fragment A
     boost::shared_ptr<Matrix> D_A;
     //the density matrix for fragment B
     boost::shared_ptr<Matrix> D_B;
     //The original Hcore matrix
     boost::shared_ptr<Matrix> Hcore_;
     //HAinB
     boost::shared_ptr<Matrix> HAinB_;
     //The Manby embedding method without projectors--not used right now
     void ManbyEmbed();
     //Embedding method used libfock interface--more efficient version of ManbyEmbed
     void FockEmbed();
     /*Forms HAinB
     @GA - G matrix for system A
     &GAB - G matrix for overall system
     returns HAinB
     */
     boost::shared_ptr<Matrix> formHAinB(boost::shared_ptr<Matrix>& GA, boost::shared_ptr<Matrix>& GAB);
     /*Computes SCF energy for system A, B, and AB.
     @De - Density matrix for a system of your choosing 
     @H - A H matrix (Hcore or HAinB)
     @G - A G matrix (GA, GAB, or GB)
     returns E_k = tr(D_k(H_k + G_b))
     */
     double EmbedEnergy(boost::shared_ptr<Matrix>& De, boost::shared_ptr<Matrix>& H, boost::shared_ptr<Matrix>& G);
};
}}
#endif /* !EMBED_H */


