#include <fstream>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include "Embed.h"
#include <libmints/matrix.h>
#include <libmints/wavefunction.h>

INIT_PLUGIN

using namespace boost;
namespace psi{ namespace plugin_fragment {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "PLUGIN_FRAGMENT"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_double("LocalCut",0.0000);
        options.add_int("PRINT", 1);
    }

    return true;
}
Embed::Embed(boost::shared_ptr<Wavefunction> reference_wavefunction, Options& options) : 
Wavefunction(options, _default_psio_lib_)
{
    copy(reference_wavefunction);
    common_init();
    
}

extern "C" 
PsiReturnType plugin_fragment(Options& options)
{
    //Get the number of fragments from the wfn object
   
    boost::shared_ptr<Embed> embed(
        new Embed(Process::environment.wavefunction(), options));
    embed->begin(options);
    return Success;
}

void Embed::common_init() {
    long int nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    for(int h = 0; h < nirrep_; h++){
       nfzc += frzcpi_[h];
       nfzv += frzvpi_[h];
       nso +=  nsopi_[h];
       nmo +=  nmopi_[h];
       ndocc += doccpi_[h];
    }
    SharedMatrix Hcore(new Matrix("Hcore", nmo, nmo));
    SharedMatrix Dfull(new Matrix("Dfull", nmo, nmo));
    ndoccact = ndocc - nfzc;
    nvirt = nmo - ndocc;
    Dfull = wfn->Da(); 
    Hcore = wfn->H();

}


Embed::~Embed()
{
}

void Embed::begin(Options& options)
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Molecule> mol = wfn->molecule();
    int nfragments = mol->nfragments();

    // nfragments = 1
    // localize on atoms

    if (options.get_int("PRINT") > 0){
        outfile->Printf("\n\n         ---------------------------------------------------------");
        outfile->Printf("\n                            ORBITAL LOCALIZATION");
        outfile->Printf("\n                        by Francesco A. Evangelista");
        outfile->Printf("\n         ---------------------------------------------------------");

        outfile->Printf("\n\n  ==> Details <==");
        outfile->Printf("\n\n    Number of fragments: %d",nfragments);
    }

    if (nfragments == 1){
        localize_on_atoms(options,wfn);
    }else{
        localize_on_fragment(options,wfn);
    }
    //    ManbyEmbed();
        FockEmbed();
}

}} // End namespaces



