#ifndef CALQ_PROGRAM_OPTIONS_H_
#define CALQ_PROGRAM_OPTIONS_H_

#include <string>

namespace calq {

class ProgramOptions {
   public:
    ProgramOptions() = delete;
    ProgramOptions(int argc, char *argv[]);
    ProgramOptions(const ProgramOptions &) = delete;
    ProgramOptions &operator=(const ProgramOptions &) = delete;
    ProgramOptions(ProgramOptions &&) = delete;
    ProgramOptions &operator=(ProgramOptions &&) = delete;
    ~ProgramOptions() = default;

    bool decompress;
    bool force;
    bool help;

   private:
    void processCommandLine(int argc, char *argv[]);
};

}  // namespace calq

#endif  // CALQ_PROGRAM_OPTIONS_H_
