#include <csignal>
#include <iostream>
#include <sstream>
#include <vector>
#include <calq/exceptions.h>
#include <util/log.h>
#include "decode.h"
#include "encode.h"
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
        util::Log::setLevel(util::Log::Level::INFO);
        cip::ProgramOptions programOptions(argc, argv);
        if (programOptions.help) {
            return 0;
        }
        util::Log::setLevel(static_cast<util::Log::Level>(programOptions.logLevel));
        LOG_INFO("command line: " + commandLineStr(argc, argv));

        if (programOptions.decompress) {
            decode(programOptions);
        } else {
            encode(programOptions);
        }
    } catch (const calq::ErrorException& errorException) {
        LOG_ERROR("CALQ error: " + errorException.whatStr());
        return -1;
    } catch (const std::exception& stdException) {
        LOG_ERROR("error: " + std::string(stdException.what()));
        return -1;
    } catch (...) {
        LOG_ERROR("fatal: unknown error occurred");
        return -1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    // Fire up main method
    int mainRc = calq_main(argc, argv);
    if (mainRc != 0) {
        LOG_ERROR("failed to run cip");
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
    LOG_DEBUG("exiting with return code: " + std::to_string(callerRc));
    return callerRc;
}
