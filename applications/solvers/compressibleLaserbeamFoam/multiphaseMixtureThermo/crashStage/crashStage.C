/*---------------------------------------------------------------------------*\
    crashStage: see crashStage.H
\*---------------------------------------------------------------------------*/

#include "crashStage.H"

#include <csignal>
#include <cstring>
#include <unistd.h>

namespace
{
    char stageBuf_[256] = "startup";
    int procNo_ = 0;
    bool installed_ = false;

    struct sigaction oldFpe_;
    struct sigaction oldSegv_;
    struct sigaction oldBus_;

    void writeStr(const char* s)
    {
        const ssize_t r = ::write(2, s, std::strlen(s));
        (void)r;
    }

    void writeInt(long v)
    {
        char b[32];
        int i = 31;
        b[i] = '\0';
        const bool neg = (v < 0);
        unsigned long u = neg ? static_cast<unsigned long>(-v) : v;
        do
        {
            b[--i] = char('0' + (u % 10));
            u /= 10;
        } while (u && i > 1);
        if (neg)
        {
            b[--i] = '-';
        }
        writeStr(b + i);
    }

    const char* fpeCodeName(int code)
    {
        switch (code)
        {
            case FPE_FLTDIV: return "division by zero";
            case FPE_FLTOVF: return "overflow (result too large)";
            case FPE_FLTUND: return "underflow";
            case FPE_FLTINV:
                return "invalid operation (NaN, sqrt/log of negative, 0/0, inf-inf)";
            case FPE_FLTRES: return "inexact result";
            case FPE_INTDIV: return "integer division by zero";
            case FPE_INTOVF: return "integer overflow";
            default: return "unknown FPE code";
        }
    }

    void handler(int sig, siginfo_t* info, void* ctx)
    {
        writeStr("\n[crashDiagnostics] processor ");
        writeInt(procNo_);
        writeStr(" caught signal ");
        writeInt(sig);

        if (sig == SIGFPE)
        {
            writeStr(" (SIGFPE: ");
            writeStr(fpeCodeName(info ? info->si_code : 0));
            writeStr(")");
        }
        else if (sig == SIGSEGV)
        {
            writeStr(" (SIGSEGV: invalid memory access)");
        }
        else if (sig == SIGBUS)
        {
            writeStr(" (SIGBUS: bus error)");
        }

        writeStr("\n[crashDiagnostics] last stage reached on this processor: ");
        writeStr(stageBuf_);
        writeStr("\n[crashDiagnostics] OpenFOAM stack trace follows\n");

        struct sigaction* old =
            (sig == SIGFPE) ? &oldFpe_
          : (sig == SIGSEGV) ? &oldSegv_
          : &oldBus_;

        // Hand over to the previous (OpenFOAM) handler: it prints the stack
        // and re-raises with the default action.
        if (old->sa_flags & SA_SIGINFO)
        {
            if (old->sa_sigaction)
            {
                old->sa_sigaction(sig, info, ctx);
                return;
            }
        }
        else if (old->sa_handler != SIG_DFL && old->sa_handler != SIG_IGN)
        {
            old->sa_handler(sig);
            return;
        }

        // No previous handler: restore default and re-raise
        ::sigaction(sig, old, nullptr);
        ::raise(sig);
    }

    void installOne(int sig, struct sigaction* old)
    {
        struct sigaction act;
        std::memset(&act, 0, sizeof(act));
        act.sa_sigaction = handler;
        act.sa_flags = SA_SIGINFO;
        sigemptyset(&act.sa_mask);
        ::sigaction(sig, &act, old);
    }
}


void Foam::crashStage::set(const char* stage)
{
    std::strncpy(stageBuf_, stage, sizeof(stageBuf_) - 1);
    stageBuf_[sizeof(stageBuf_) - 1] = '\0';
}


void Foam::crashStage::set(const std::string& stage)
{
    set(stage.c_str());
}


const char* Foam::crashStage::get()
{
    return stageBuf_;
}


void Foam::crashStage::install(int procNo)
{
    procNo_ = procNo;

    if (installed_)
    {
        return;
    }
    installed_ = true;

    installOne(SIGFPE, &oldFpe_);
    installOne(SIGSEGV, &oldSegv_);
    installOne(SIGBUS, &oldBus_);
}


// ************************************************************************* //
