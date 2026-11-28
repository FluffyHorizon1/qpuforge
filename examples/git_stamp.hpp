// Source stamp of the review and study programs: the commit the binary was built from, and whether tracked files
// differed from that commit when it was built. review_io.hpp writes both into the provenance block of every output
// ("source_commit", "source_tree_dirty").
//
// A CMake build takes the stamp on every build, not when the build directory is configured. The always-run target
// qpuforge_git_stamp_update calls cmake/git_stamp.cmake, which writes the two definitions into git_stamp.cpp in the
// build tree (only when they change). That file is the one source of the static library qpuforge_git_stamp, which the
// programs link and which defines QPUFORGE_GIT_STAMP for them. No other translation unit depends on the commit, so a
// new commit costs one small compile and a relink of each program.
//
// Compiled by hand (g++ -std=c++17 -Iinclude examples/surgery_study.cpp) there is no such library, and the stamp is
// "unknown" / false. It is the same in a CMake build outside a git checkout or without git.
#pragma once

namespace review {

#ifdef QPUFORGE_GIT_STAMP
extern const char* const git_stamp_commit;   // full hash of HEAD, or "unknown"
extern const bool git_stamp_dirty;           // tracked files modified (untracked files do not count)
#else
inline constexpr const char* git_stamp_commit = "unknown";
inline constexpr bool git_stamp_dirty = false;
#endif

} // namespace review
