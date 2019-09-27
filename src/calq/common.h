#include <chrono>

namespace cip {

template <typename Diff>
void logProgress(Diff d) {
    std::cout << "took " << std::chrono::duration_cast<std::chrono::milliseconds>(d).count()
              << " ms ~= " << std::chrono::duration_cast<std::chrono::seconds>(d).count()
              << " s ~= " << std::chrono::duration_cast<std::chrono::minutes>(d).count()
              << " min ~= " << std::chrono::duration_cast<std::chrono::hours>(d).count() << " h" << std::endl;
}

}  // namespace cip