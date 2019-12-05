#include <iostream>
#include <sstream>
#include <vector>
#include "decode.h"
#include "encode.h"
#include "exceptions.h"
#include "program-options.h"

static std::string commandLineStr(int argc, char* argv[]) {
    std::vector<std::string> args(argv, (argv + argc));
    std::stringstream commandLine;
    for (const auto& arg : args) {
        commandLine << arg << " ";
    }
    return commandLine.str();
}

static int calq_codec_main(int argc, char* argv[]) {
    try {
        calq::ProgramOptions programOptions(argc, argv);
        if (programOptions.help) {
            return 0;
        }
        std::cout << "command line: " + commandLineStr(argc, argv) << std::endl;

        if (programOptions.decompress) {
            calq::decode(programOptions);
        } else {
            calq::encode(programOptions);
        }
    } catch (const calq::ErrorException& errorException) {
        std::cerr << "error: calq: " << errorException.what() << std::endl;
        return -1;
    } catch (const std::exception& stdException) {
        std::cerr << "error: " << stdException.what() << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "error: unknown error occurred" << std::endl;
        return -1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    int rc = calq_codec_main(argc, argv);
    if (rc != 0) {
        std::cerr << "error: failed to run calq-codec" << std::endl;
    }

    return ((rc == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
}
