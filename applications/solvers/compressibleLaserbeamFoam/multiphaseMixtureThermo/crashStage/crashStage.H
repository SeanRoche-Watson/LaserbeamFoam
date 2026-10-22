/*---------------------------------------------------------------------------*\
Namespace
    Foam::crashStage

Description
    Lightweight "breadcrumb" for locating crashes.

    Code marks where it is with crashStage::set("..."). If the process then
    dies on SIGFPE / SIGSEGV / SIGBUS, a handler chained in front of
    OpenFOAM's own prints, from the processor that died:

        [crashDiagnostics] processor 125 caught signal 8
            (SIGFPE: overflow (result too large))
        [crashDiagnostics] last stage reached on this processor:
            solveAlphas: phase change (metal4,metal4vapour) volume clamp

    and then hands over to OpenFOAM's handler, which prints the stack trace
    as before. Everything in the handler is async-signal-safe (write(2)).

    set() costs one bounded string copy, so it can be called freely, but
    not inside per-cell loops.

SourceFiles
    crashStage.C

\*---------------------------------------------------------------------------*/

#ifndef crashStage_H
#define crashStage_H

#include <string>

namespace Foam
{
namespace crashStage
{
    //- Record the current stage (truncated to 255 chars)
    void set(const char* stage);

    //- Record the current stage (truncated to 255 chars)
    void set(const std::string& stage);

    //- The last recorded stage
    const char* get();

    //- Chain the reporting handler in front of the existing SIGFPE,
    //  SIGSEGV and SIGBUS handlers. Call once, after OpenFOAM has set up
    //  its own (i.e. after setRootCase / createTime).
    void install(int procNo);

} // End namespace crashStage
} // End namespace Foam

#endif
