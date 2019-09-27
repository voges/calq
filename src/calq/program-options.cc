#include "program-options.h"
#include <cli11@13becad/CLI11.hpp>
#include <sstream>

namespace calq {

ProgramOptions::ProgramOptions(int argc, char *argv[]) : decompress(false), force(false), help(false) {
    processCommandLine(argc, argv);
}

void ProgramOptions::processCommandLine(int argc, char *argv[]) {
    CLI::App app("CALQ");

    app.add_flag("-d,--decompress", decompress, "Decompress");
    app.add_flag("-f,--force", force, "Force overwriting of output file(s)");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        if (app.count("--help")) {
            // Set our internal help flag so that the caller can act on that
            help = true;
            app.exit(e);
            return;
        } else {
            app.exit(e);
            throw std::runtime_error("command line parsing failed");
        }
    }
}

}  // namespace calq
