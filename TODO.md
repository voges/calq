- [ ] Use log::prefix() for logging
- [ ] Clean up exceptions (use only simple derived class - w/o all the special things)
- [ ] Re-structure: cmake-build-*/, doc/, thirdparty/, src/, test/, data/
- [ ] u(32) issue in GABAC
- [ ] Make CIP application (source/apps/cip/) work again (including GABAC and rANS)
- [ ] Think about a logging solution (singleton class from Genie's ureads-encoder?)
- [ ] Format log and error messages uniformly
- [ ] Check whether both SamRecord and MinSamRecord classes are needed
- [ ] Check comments
  - [ ] Sentences should start with an uppercase letter and end with a full stop.
  - [ ] Statements should start with an uppercase letter and not end with a full stop.
- [ ] Resolve FIXME's and TODO's
- [ ] Resolve clang-tidy warnings

- [ ] Check for correct usage of assertions versus exceptions
- [ ] Follow the C++ Core Guidelines (https://isocpp.github.io/CppCoreGuidelines), in particular C.21



# Final checks

- [ ] Update README.md once everything is working
- [ ] Check everything with Valgrind
