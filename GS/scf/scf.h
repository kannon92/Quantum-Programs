#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/wavefunction.h>
#include <vector>

#ifndef Embed_H
#define Embed_H

namespace psi { namespace scf{
class Scfclass : public Wavefunction
{
protected:
      ///Core Hamiltonian matrix
      SharedMatrix hMat_;
      ///The overlap matrix 
      SharedMatrix sMat_;
      //A standard vector that stores the tei in an array of size nmo^4
      std::vector<double> teiV_;
      std::vector<double> eps_;
      std::vector<double> fermidirac_;
     
      //General variables for use in SCF code
      int nbf_;
      int nmo_;
      //The number of occupied orbitals - assumes SCF reference
      int nocc_;
      double compute_energy();
      //If FTHF, computes the density via n_i C C^T but sum is over all ao.  
      boost::shared_ptr<Matrix> frac_occupation(SharedMatrix,int &iter, bool &t_done);
      double bisection(std::vector<double>&, double T);
      double occ_vec(std::vector<double>& bisect, double ef, double T);
      //The fermi level (n_i = N - solved using bisection method)
      double ef_ = 0.0;
      double energy = 0.0;
private:
public:
      void scf_energy();
      void scf_iteration();
      Scfclass(boost::shared_ptr<Wavefunction> wfn, Options& Options);
      void ao_ints(Options& options);
      virtual ~Scfclass();
};

}}
#endif /* !EMBED_H */

