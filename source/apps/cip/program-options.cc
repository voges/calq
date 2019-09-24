#include "program-options.h"
#include <util/log.h>
#include <cli11@13becad/CLI11.hpp>
#include <sstream>

namespace cip {

ProgramOptions::ProgramOptions(int argc, char *argv[]) : decompress(false), force(false), help(false), logLevel(static_cast<int>(util::Log::Level::INFO)) { processCommandLine(argc, argv); }

void ProgramOptions::processCommandLine(int argc, char *argv[]) {
    CLI::App app("CALQ");

    app.add_flag("-d,--decompress", decompress, "Decompress");
    app.add_flag("-f,--force", force, "Force overwriting of output file(s)");
    app.add_option("-l,--log-level", logLevel, "Log level (0 = fatal, 1 = error, 2 = warn, 3 = info, 4 = debug, 5 = trace)");

    try {
        app.parse(argc, argv);
        validate();
    } catch (const CLI::ParseError &e) {
        if (app.count("--help")) {
            // Set our internal help flag so that the caller can act on that
            help = true;
            app.exit(e);
            return;
        } else {
            throw std::runtime_error("command line parsing failed");
            app.exit(e);
        }
    }
}

void ProgramOptions::validate() {
    validateCommon();
}

void ProgramOptions::validateCommon() {
    // decompress

    // force

    // logLevel
    if (!util::isValidLogLevel(logLevel)) {
        throw std::runtime_error("invalid log level: " + std::to_string(logLevel));
    }
}

}  // namespace cip
