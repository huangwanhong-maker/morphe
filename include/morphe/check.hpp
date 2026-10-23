// Checks of exact results, reported line by line.
//
// morphe depends on no other library of the project (CLAUDE.md, invariant 1), so its checks
// cannot use calx's verification; this header prints the same report, so that the runner of
// The Science Infrastructures adds morphe's checks to the others (decision 0004, §2.6):
//
//   ok    name: detail
//   FAIL  name: detail
//   skip  what
//
//   41 of 41 checks passed.
//
// Most of morphe's checks are exact, and are stated with require(); check() compares a
// floating-point deviation with a limit, for the few structures built from real numbers.

#pragma once

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>

namespace morphe {

class Suite {
 public:
  bool require(std::string_view name, bool condition, std::string_view detail = "") {
    count_++;
    if (!condition) failures_++;
    std::cout << (condition ? "ok  " : "FAIL") << "  " << name << (detail.empty() ? "" : ": ") << detail << "\n";
    return condition;
  }
  // Two values that must be equal, written in the report as the computed and the expected.
  template <class A, class B>
  bool equal(std::string_view name, const A& got, const B& want) {
    const bool ok = got == want;
    return require(name, ok, ok ? text(got) : text(got) + ", expected " + text(want));
  }
  // value: a deviation; it passes when finite and |value| <= limit.
  bool check(std::string_view name, double value, double limit) {
    char buf[96];
    std::snprintf(buf, sizeof buf, "%.3e (limit %g)", value, limit);
    return require(name, std::isfinite(value) && std::fabs(value) <= limit, buf);
  }
  void skip(std::string_view what) { std::cout << "skip  " << what << "\n"; }
  // Prints the count; the exit status of the program.
  int done() {
    std::cout << "\n" << (count_ - failures_) << " of " << count_ << " checks passed.\n";
    std::cout.flush();
    return failures_ ? 1 : 0;
  }
  int failures() const { return failures_; }

 private:
  template <class T>
  static std::string text(const T& x) {
    if constexpr (std::is_same_v<T, bool>) return x ? "true" : "false";
    else if constexpr (requires { std::string(x); }) return std::string(x);
    else if constexpr (requires { x.str(); }) return x.str();
    else if constexpr (requires { std::to_string(x); }) return std::to_string(x);
    else if constexpr (requires { x.begin(); x.end(); }) {
      std::string s = "(";
      for (const auto& e : x) s += (s.size() > 1 ? ", " : "") + text(e);
      return s + ")";
    } else return "(value)";
  }
  int count_ = 0, failures_ = 0;
};

}  // namespace morphe
