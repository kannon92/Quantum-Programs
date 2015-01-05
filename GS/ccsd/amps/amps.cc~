/*
 *@BEGIN LICENSE
 *
 * amps by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <psifiles.h>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace amps {

class Amps : public Wavefunction
{
public:
    Amps(boost::shared_ptr<Wavefunction> reference_wavefunction, Options& options);
    virtual ~Amps();

    double compute_energy();

private:
    Dimension nvirtpi_;
    void common_init();
};

Amps::Amps(boost::shared_ptr<Wavefunction> reference_wavefunction, Options& options)
    : Wavefunction(options, _default_psio_lib_)
{
    Process::environment.set_wavefunction(reference_wavefunction);
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

Amps::~Amps()
{
}

void Amps::common_init()
{
    nso_     = reference_wavefunction_->nso();
    nirrep_  = reference_wavefunction_->nirrep();
    nmo_     = reference_wavefunction_->nmo();
    soccpi_  = reference_wavefunction_->soccpi();
    doccpi_  = reference_wavefunction_->doccpi();
    frzcpi_  = reference_wavefunction_->frzcpi();
    frzvpi_  = reference_wavefunction_->frzvpi();
    nmopi_   = reference_wavefunction_->nmopi();
    nsopi_   = reference_wavefunction_->nsopi();
    nvirtpi_ = nsopi_ - frzcpi_ - frzvpi_ - doccpi_;
}

double Amps::compute_energy()
{
    int v = nvirtpi_[0];
    int o = doccpi_[0];
    //SharedVector t2(new Vector("t2", v*v*o*o));
    //double* t2p = t2->pointer();
    double *t2 = new double[o*o*v*v];
    boost::shared_ptr<PSIO> psio (new PSIO());
    
    psio->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);
    psio->read_entry(PSIF_CC_TAMPS, "T2", (char*)&t2[0], o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_T2,1);

    return 0.0;
}

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "AMPS"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C"
PsiReturnType amps(Options& options)
{
    Amps wave(Process::environment.wavefunction(), options);
    wave.compute_energy();

    return Success;
}

}} // End namespaces

