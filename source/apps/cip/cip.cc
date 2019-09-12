#include <csignal>
#include <iostream>
#include <sstream>
#include <vector>
#include "calq/errors.h"
#include "program-options.h"

static std::string commandLineStr(int argc, char* argv[]) {
    std::vector<std::string> args(argv, (argv + argc));
    std::stringstream commandLine;
    for (const auto& arg : args) {
        commandLine << arg << " ";
    }
    return commandLine.str();
}

static int calq_main(int argc, char* argv[]) {
    try {
        std::cout << "cip: command line: " << commandLineStr(argc, argv) << std::endl;
        cip::ProgramOptions programOptions(argc, argv);
    } catch (const calq::ErrorException& errorException) {
        std::cerr << "cip: CALQ error: " << errorException.what() << std::endl;
        return -1;
    } catch (const std::exception& stdException) {
        std::cerr << "cip: error: " << stdException.what() << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "cip: fatal: unknown error occurred" << std::endl;
        return -1;
    }

    return 0;
}

extern "C" void handleSignal(int sig) {
    // Ignore the signal
    std::signal(sig, SIG_IGN);

    // Get signal string and log it
    std::string signalString;
    switch (sig) {
        case SIGABRT:
            signalString = "SIGABRT";
            break;
        case SIGFPE:
            signalString = "SIGFPE";
            break;
        case SIGILL:
            signalString = "SIGILL";
            break;
        case SIGINT:
            signalString = "SIGINT";
            break;
        case SIGSEGV:
            signalString = "SIGSEGV";
            break;
        case SIGTERM:
            signalString = "SIGTERM";
            break;
        default:
            signalString = "unknown";
            break;
    }
    std::cout << "cip: caught signal: " << sig << " (" << signalString << ")";

    // Invoke the default signal action
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

int main(int argc, char* argv[]) {
    // Install signal handler for the following signal types:
    //   SIGABRT  abnormal termination condition, as is e.g. initiated by std::abort()
    //   SIGFPE   erroneous arithmetic operation such as divide by zero
    //   SIGILL   invalid program image, such as invalid instruction
    //   SIGINT   external interrupt, usually initiated by the user
    //   SIGSEGV  invalid memory access (segmentation fault)
    //   SIGTERM  termination request, sent to the program
    std::signal(SIGABRT, handleSignal);
    std::signal(SIGFPE, handleSignal);
    std::signal(SIGILL, handleSignal);
    std::signal(SIGINT, handleSignal);
    std::signal(SIGSEGV, handleSignal);
    std::signal(SIGTERM, handleSignal);

    // Fire up main method
    int mainRc = calq_main(argc, argv);
    if (mainRc != 0) {
        std::cerr << "cip: error: failed to run cip" << std::endl;
    }

    // The C standard makes no guarantees as to when output to stdout or stderr (standard error) is actually flushed. If
    // e.g. stdout is directed to a file and an error occurs while flushing the data (after program termination), then
    // the output may be lost. Thus we explicitly flush stdout and stderr. On failure, we notify the OS by returning
    // with EXIT_FAILURE.
    if (fflush(stdout) == EOF) {
        return EXIT_FAILURE;
    }
    if (fflush(stderr) == EOF) {
        return EXIT_FAILURE;
    }

    // Return to the caller
    int callerRc = ((mainRc == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
    std::cout << "cip: exiting with return code: " << callerRc << std::endl;
    return callerRc;
}
