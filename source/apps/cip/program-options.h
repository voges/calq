#ifndef CIP_PROGRAM_OPTIONS_H_
#define CIP_PROGRAM_OPTIONS_H_

#include <string>

namespace cip {

class ProgramOptions {
   public:
    ProgramOptions(int argc, char *argv[]);

    bool decompress;
    bool force;
    bool help;
    int logLevel;

   private:
    void processCommandLine(int argc, char *argv[]);
    void validate();
    void validateCommon();
};

}  // namespace cip

#endif  // CIP_PROGRAM_OPTIONS_H_
