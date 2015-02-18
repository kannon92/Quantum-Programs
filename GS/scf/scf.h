#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/wavefunction.h>

#ifndef Embed_H
#define Embed_H

namespace psi { namespace scf{
class Scfclass : public Wavefunction
{
protected:
      SharedMatrix hMat_;
      SharedMatrix sMat_;
      std::vector<double> teiV_;
      int nbf_;
      int nmo_;
      int nocc_;
      double compute_energy();
      boost::shared_ptr<Matrix> frac_occupation(SharedMatrix, SharedVector, int &iter, bool &t_done);
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

